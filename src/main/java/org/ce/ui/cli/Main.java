package org.ce.ui.cli;

import static org.ce.domain.cluster.AllClusterData.ClusterData;

import org.ce.domain.cluster.*;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.engine.cvm.CVMFreeEnergy;
import org.ce.domain.engine.mcs.MCSEngine;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.Workspace.SystemId;
import org.ce.storage.Workspace;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.InputLoader;
import org.ce.workflow.CalculationService;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicRequest;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.util.List;
import java.io.PrintStream;

/**
 * Main entry point for the CE Workbench CLI.
 *
 * Usage:
 *   ./gradlew run                                  -- full pipeline, built-in defaults
 *   ./gradlew run --args="type1a"                  -- cluster identification only
 *   ./gradlew run --args="type1b"                  -- Hamiltonian scaffold only
 *   ./gradlew run --args="type2"                   -- thermodynamic calculation only
 *   ./gradlew run --args="all  Nb-Ti BCC_A2 T"     -- explicit system, all modes
 *   ./gradlew run --args="all  Nb-Ti BCC_A2 T_CVCF" -- explicit system with CVCF basis
 *   ./gradlew run --args="scaffold_cvcf Nb-Ti BCC_A2 T_CVCF" -- scaffold CEC using CVCF basis
 *   ./gradlew run --args="type2 Nb-Ti BCC_A2 T"    -- explicit system, type-2 only
 *   ./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 0.5"  -- calculate with minimization
 *   ./gradlew run --args="calc_min Nb-Ti BCC_A2 T_CVCF 1000 0.5"  -- CVCF calc_min
 *   ./gradlew run --args="calc_fixed Nb-Ti BCC_A2 T 1000 0.5 0.1 0.2 0.3 0.4"  -- calculate with fixed CFs
 *   ./gradlew run --args="calc_fixed Nb-Ti BCC_A2 T_CVCF 1000 0.5 0.1 0.2 0.3 0.4"  -- CVCF calc_fixed
 *
 * Arguments:  mode  [elements  structure  model]
 *   mode       : type1a | type1b | type2 | all | calc_min | calc_fixed  (default: all)
 *   elements   : element pair, e.g. Nb-Ti  (default: A-B)
 *   structure  : structure ID,  e.g. BCC_A2 (default: BCC_B2)
 *   model      : model ID,      e.g. T      (default: T)
 *   for calc_min: temp comp
 *   for calc_fixed: temp comp cf1 cf2 ...
 */
public class Main {

    /**
     * Registry of supported CVCF transformations (T model of CVM).
     *
     * <p>Maps (structure, numComponents) pairs to their corresponding CVCF basis transformations
     * as defined in the T model of CVM. All transformation matrices are sourced from the
     * {@link BccA2TModelCvCfTransformations} class, which encodes the orthogonal-to-CVCF basis
     * transformations derived from first-principles theory for the T model.
     *
     * <p>Supported systems (BCC_A2 structure, T model only):
     * <ul>
     *   <li>Binary BCC_A2_T (2-component): 6 orthogonal CFs → 6 CVCF CFs (4 non-point + 2 point)</li>
     *   <li>Ternary BCC_A2_T (3-component): 21 orthogonal CFs → 21 CVCF CFs (18 non-point + 3 point)</li>
     *   <li>Quaternary BCC_A2_T (4-component): 55 orthogonal CFs → 55 CVCF CFs (51 non-point + 4 point)</li>
     * </ul>
     *
     * <p>Format: key = "{STRUCTURE}|{NCOMP}" e.g. "BCC_A2|2"
     * <p>Source matrices: BccA2TModelCvCfTransformations.{BINARY_T, TERNARY_T, QUATERNARY_T}
     */

    public static void main(String[] args) {


        // Calculate with minimization
        if (args.length > 0 && args[0].equals("calc_min")) {
            if (args.length < 6) {
                System.err.println("Usage: calc_min <elements> <structure> <model> <temp> <comp1> [<comp2> <comp3> ...]");
                System.exit(1);
            }
            String elements = args[1];
            String structure = args[2];
            String model = args[3];
            double temp = Double.parseDouble(args[4]);
            
            // Handle variable-length composition array
            double[] composition;
            if (args.length == 6) {
                // Binary: single composition value given
                double comp = Double.parseDouble(args[5]);
                composition = new double[]{1 - comp, comp};
            } else {
                // Ternary/Quaternary: all composition values given explicitly
                composition = new double[args.length - 5];
                for (int i = 5; i < args.length; i++) {
                    composition[i - 5] = Double.parseDouble(args[i]);
                }
            }
            runCalcMin(elements, structure, model, temp, composition);
            return;
        }

        // Calculate with fixed CFs
        if (args.length > 0 && args[0].equals("calc_fixed")) {
            String elements = args.length > 1 ? args[1] : "";
            int K = elements.isBlank() ? 2 : elements.split("-").length;
            if (args.length < 6 + K) {
                System.err.println("Usage: calc_fixed <elements> <structure> <model> <temp> <comp1> [<comp2> ...compK] <cf1> <cf2> ...");
                System.err.println("  For a " + K + "-component system provide " + K + " composition values followed by CF values.");
                System.exit(1);
            }
            String structure = args[2];
            String model = args[3];
            double temp = Double.parseDouble(args[4]);
            double[] composition = new double[K];
            for (int i = 0; i < K; i++) {
                composition[i] = Double.parseDouble(args[5 + i]);
            }
            double[] cfs = new double[args.length - 5 - K];
            for (int i = 0; i < cfs.length; i++) {
                cfs[i] = Double.parseDouble(args[5 + K + i]);
            }
            runCalcFixed(elements, structure, model, temp, composition, cfs);
            return;
        }

        // Parse positional arguments — same three values as the GUI's system context
        String mode      = args.length > 0 ? args[0].toLowerCase() : "all";
        String elements  = args.length > 1 ? args[1] : "A-B";
        String structure = args.length > 2 ? args[2] : "BCC_B2";
        String model     = args.length > 3 ? args[3] : "T";

        if (!mode.equals("type1a") && !mode.equals("type1b")
                && !mode.equals("type2") && !mode.equals("all")
                && !mode.equals("calc_min") && !mode.equals("calc_fixed")
                && !mode.equals("scaffold_cvcf")) {
            System.err.println("Unknown mode: " + mode);
            System.err.println("Usage: [mode] [elements] [structure] [model]");
            System.err.println("  mode: type1a | type1b | type2 | all | calc_min | calc_fixed | scaffold_cvcf");
            System.exit(1);
        }

        if (mode.equals("scaffold_cvcf")) {
            scaffoldCecCvcf(elements, structure, model);
            return;
        }

        SystemId system       = new SystemId(elements, structure, model);
        String   CLUSTER_ID   = system.clusterId();
        String   HAMILTONIAN_ID = system.hamiltonianId();
        int numComponents = elements.split("-").length;

        System.out.println("================================================================================");
        System.out.println("                    CE THERMODYNAMICS WORKBENCH");
        System.out.println("================================================================================");
        System.out.println("Mode      : " + mode);
        System.out.println("System    : " + elements + "  /  " + structure + "  /  " + model);
        System.out.println("Cluster ID: " + CLUSTER_ID + "   Hamiltonian ID: " + HAMILTONIAN_ID);

        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            ClusterDataStore clusterStore = appCtx.getClusterStore();
            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();
            CalculationService service = appCtx.getCalculationService();

            // ------------------------------------------------------------------
            // TYPE-1a: Cluster Identification + Save
            // ------------------------------------------------------------------
            if (mode.equals("type1a") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1a: Cluster Identification ===\n");

                String disClusFile = "clus/BCC_A2-T.txt";
                String disSymGroup = "BCC_A2-SG";
                String ordClusFile = "clus/BCC_B2-T.txt";
                String ordSymGroup = "BCC_B2-SG";
                if ("BCC_A2".equalsIgnoreCase(structure)) {
                    ordClusFile = "clus/BCC_A2-T.txt";
                    ordSymGroup = "BCC_A2-SG";
                }

                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .disorderedClusterFile(disClusFile)
                        .orderedClusterFile(ordClusFile)
                        .disorderedSymmetryGroup(disSymGroup)
                        .orderedSymmetryGroup(ordSymGroup)
                        .numComponents(numComponents)
                        .build();

                AllClusterData clusterData = appCtx.runType1a(CLUSTER_ID, config, System.out::println);
                System.out.println("Identification complete: " + clusterData.getSummary());
                System.out.println("Saved: " + appCtx.getWorkspace().clusterDataFile(CLUSTER_ID));
            }

            // ------------------------------------------------------------------
            // TYPE-1b: Scaffold empty Hamiltonian from saved cluster data
            // ------------------------------------------------------------------
            if (mode.equals("type1b") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1b: Hamiltonian Scaffold ===\n");

                if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                    System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                    System.out.println("  Delete " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID) + " to re-scaffold.");
                } else {
                    appCtx.runType1b(HAMILTONIAN_ID, elements, structure, model);
                    System.out.println("Saved: " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID));
                    System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
                }
            }

            // ------------------------------------------------------------------
            // TYPE-2: Thermodynamic Calculation (CVM temperature scan)
            // ------------------------------------------------------------------
            if (mode.equals("type2") || mode.equals("all")) {

                System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");

                double[] composition = {0.5, 0.5};
                if (numComponents > 2) {
                    composition = new double[numComponents];
                    double x = 1.0 / numComponents;
                    for (int i = 0; i < numComponents; i++) {
                        composition[i] = x;
                    }
                }
                double tStart = 1000.0;
                double tEnd   = 1000.0;
                double tStep  = 100.0;

                System.out.println("System      : " + CLUSTER_ID + " / " + HAMILTONIAN_ID);
                System.out.println("Composition : " + java.util.Arrays.toString(composition));
                System.out.println("T range     : " + tStart + " K to " + tEnd + " K, step " + tStep + " K\n");

                List<ThermodynamicResult> results = service.runLineScanTemperature(
                        CLUSTER_ID, HAMILTONIAN_ID, composition, tStart, tEnd, tStep, "CVM"
                );

                System.out.println(String.format("  %-8s  %-16s  %-16s",
                        "T (K)", "G (J/mol)", "H (J/mol)"));
                System.out.println("  " + "-".repeat(48));
                for (ThermodynamicResult r : results) {
                    System.out.println(String.format("  %-8.1f  %-16.4f  %-16.4f",
                            r.temperature, r.gibbsEnergy, r.enthalpy));
                }
            }

        } catch (Exception e) {
            System.err.println("\nError: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("\n================================================================================");
    }

    private static void printSupportedCvcfStructures(PrintStream out) {
        out.println(CvCfBasis.Registry.INSTANCE.supportedSummary());
        out.println();
    }

    private static void scaffoldCecCvcf(String elements, String structure, String model) {
        int numComponents = elements.split("-").length;

        if (!CvCfBasis.Registry.INSTANCE.isSupported(structure, numComponents)) {
            System.err.println("CVCF basis not registered for structure '" + structure + "' with " + numComponents + " components.");
            System.err.println();
            printSupportedCvcfStructures(System.err);
            System.exit(1);
        }

        var basis = CvCfBasis.Registry.INSTANCE.get(structure, numComponents);
        System.out.println("Using CVCF basis for " + structure + " (" + numComponents + " components)");
        System.out.println("Clustering to CVCF transformation matrix (orthogonal->CVCF):");
        System.out.println(basis);

        SystemId system = new SystemId(elements, structure, model);
        String CLUSTER_ID = system.clusterId();
        String HAMILTONIAN_ID = system.hamiltonianId();

        System.out.println("Scaffolding CEC for CVCF basis:");
        System.out.println("  Elements      : " + elements);
        System.out.println("  Structure     : " + structure);
        System.out.println("  Model         : " + model);
        System.out.println("  Cluster ID    : " + CLUSTER_ID);
        System.out.println("  Hamiltonian ID: " + HAMILTONIAN_ID);

        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();

            if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                System.out.println("  Delete " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID) + " to re-scaffold.");
            } else {
                appCtx.runType1b(HAMILTONIAN_ID, elements, structure, model);
                System.out.println("Saved: " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID));
                System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
            }
        } catch (Exception e) {
            System.err.println("Error scaffolding CVCF: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static double[] extractEciFromEntry(CECEntry entry, double temperature) {
        double[] eci = new double[entry.ncf];
        for (int i = 0; i < entry.ncf && i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i] = term.a + term.b * temperature;
        }
        return eci;
    }

    private static void runCalcMin(String elements, String structure, String model, double temp, double[] composition) {
        try {
            // Set up context (same as GUI)
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            CalculationService service = appCtx.getCalculationService();

            // Derive IDs from system
            SystemId system = new SystemId(elements, structure, model);
            String clusterId = system.clusterId();
            String hamiltonianId = system.hamiltonianId();

            // Run calculation through CalculationService (same path as GUI)
            ThermodynamicResult result = service.runSinglePoint(
                    new ThermodynamicRequest(clusterId, hamiltonianId, temp, composition, "CVM",
                            System.out::println, null));

            // Format and display results
            System.out.println("System: " + hamiltonianId);
            System.out.println("Temperature: " + temp + " K");
            System.out.print("Composition: [");
            for (int i = 0; i < composition.length; i++) {
                System.out.printf("%.4f", composition[i]);
                if (i < composition.length - 1) System.out.print(", ");
            }
            System.out.println("]");
            System.out.println();

            System.out.println("EQUILIBRIUM PROPERTIES (after minimization)");
            System.out.println("-".repeat(80));
            System.out.printf("Gibbs Energy (G):     %20.10f J/mol%n", result.gibbsEnergy);
            System.out.printf("Enthalpy (H):         %20.10f J/mol%n", result.enthalpy);
            System.out.printf("Entropy (S):          %20.10f J/(mol·K)%n",
                    (result.enthalpy - result.gibbsEnergy) / temp);
            System.out.println();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void runCalcFixed(String elements, String structure, String model, double temp, double[] composition, double[] cfs) {
        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            SystemId system = new SystemId(elements, structure, model);
            String CLUSTER_ID = system.clusterId();
            String rawId = system.hamiltonianId();
            String preferredId = rawId + "_CVCF";
            
            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();
            ClusterDataStore clusterStore = appCtx.getClusterStore();
            
            String HAMILTONIAN_ID = hamiltonianStore.exists(preferredId) ? preferredId : rawId;
            
            AllClusterData allData = clusterStore.load(CLUSTER_ID);
            CECEntry entry = hamiltonianStore.load(HAMILTONIAN_ID);
            double[] eci = extractEciFromEntry(entry, temp);

            // Extract CVM parameters from AllClusterData
            ClusterIdentificationResult disResult = allData.getDisorderedClusterResult();
            CFIdentificationResult cfResult = allData.getDisorderedCFResult();
            CMatrix.Result cmatResult = allData.getCMatrixResult();

            List<Double> mhdis = disResult.getDisClusterData().getMultiplicities();
            double[] kb = disResult.getKbCoefficients();
            int tcdis = disResult.getTcdis();
            int[] lc = disResult.getLc();
            double[][] mh = disResult.getMh();
            int[][] lcv = cmatResult.getLcv();
            List<List<double[][]>> cmat = cmatResult.getCmat();
            List<List<int[]>> wcv = cmatResult.getWcv();
            
            int ncf = cfResult.getNcf();
            int numComponents = composition.length;
            int tcf = ncf + numComponents;

            System.out.println("System: " + HAMILTONIAN_ID);
            System.out.println("Temperature: " + temp + " K");
            System.out.println("Composition: " + java.util.Arrays.toString(composition));
            System.out.println("Input CFs (" + cfs.length + " values):");
            for (int i = 0; i < cfs.length; i++) {
                System.out.printf("  v[%d] = %.10f%n", i, cfs[i]);
            }
            System.out.println();

            System.out.println("ECI VALUES (at " + temp + " K):");
            System.out.println("-".repeat(80));
            for (int i = 0; i < eci.length; i++) {
                System.out.printf("  ECI[%d] = %20.10f J/mol%n", i, eci[i]);
            }
            System.out.println();

            // Evaluate free energy at fixed CFs
            CVMFreeEnergy.EvalResult result = CVMFreeEnergy.evaluate(
                    cfs, composition, temp, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

            System.out.println("THERMODYNAMIC PROPERTIES:");
            System.out.println("-".repeat(80));
            System.out.printf("Gibbs Energy (G):     %20.10f J/mol%n", result.G);
            System.out.printf("Enthalpy (H):         %20.10f J/mol%n", result.H);
            System.out.printf("Entropy (S):          %20.10f J/(mol·K)%n", result.S);
            System.out.println();

            System.out.println("FIRST DERIVATIVES (Gradient):");
            System.out.println("-".repeat(80));
            System.out.println("dG/dv (Gcu):");
            for (int i = 0; i < ncf; i++) {
                System.out.printf("  Gcu[%d] = %20.10e%n", i, result.Gcu[i]);
            }
            System.out.println();
            System.out.println("dH/dv (Hcu):");
            for (int i = 0; i < ncf; i++) {
                System.out.printf("  Hcu[%d] = %20.10e%n", i, result.Hcu[i]);
            }
            System.out.println();
            System.out.println("dS/dv (Scu):");
            for (int i = 0; i < ncf; i++) {
                System.out.printf("  Scu[%d] = %20.10e%n", i, result.Scu[i]);
            }
            System.out.println();

            System.out.println("SECOND DERIVATIVES (Hessian):");
            System.out.println("-".repeat(80));
            System.out.println("d²G/dv² (Gcuu) [non-zero values only]:");
            boolean hasNonZero = false;
            for (int i = 0; i < ncf; i++) {
                for (int j = i; j < ncf; j++) {
                    if (Math.abs(result.Gcuu[i][j]) > 1e-15) {
                        System.out.printf("  Gcuu[%d][%d] = %20.10e%n", i, j, result.Gcuu[i][j]);
                        hasNonZero = true;
                    }
                }
            }
            if (!hasNonZero) System.out.println("  (all zero)");
            System.out.println();

            System.out.println("d²H/dv² (Hcuu) [non-zero values only]:");
            hasNonZero = false;
            for (int i = 0; i < ncf; i++) {
                for (int j = i; j < ncf; j++) {
                    if (Math.abs(result.Hcuu[i][j]) > 1e-15) {
                        System.out.printf("  Hcuu[%d][%d] = %20.10e%n", i, j, result.Hcuu[i][j]);
                        hasNonZero = true;
                    }
                }
            }
            if (!hasNonZero) System.out.println("  (all zero)");
            System.out.println();

            System.out.println("d²S/dv² (Scuu) [non-zero values only]:");
            hasNonZero = false;
            for (int i = 0; i < ncf; i++) {
                for (int j = i; j < ncf; j++) {
                    if (Math.abs(result.Scuu[i][j]) > 1e-15) {
                        System.out.printf("  Scuu[%d][%d] = %20.10e%n", i, j, result.Scuu[i][j]);
                        hasNonZero = true;
                    }
                }
            }
            if (!hasNonZero) System.out.println("  (all zero)");
            System.out.println();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
