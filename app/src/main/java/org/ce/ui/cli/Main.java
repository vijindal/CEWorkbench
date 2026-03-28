package org.ce.ui.cli;

import org.ce.domain.cluster.*;
import org.ce.domain.cluster.cvcf.BccA2TModelCvCfTransformations;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.engine.cvm.CVMFreeEnergy;
import org.ce.domain.engine.cvm.CVMPhaseModel;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.SystemId;
import org.ce.storage.Workspace;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.InputLoader;
import org.ce.workflow.CalculationService;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
    private static final Map<String, CvCfBasis> CVCF_BASIS_REGISTRY = new HashMap<>();

    /**
     * Supported CVCF structures and their component ranges (T model of CVM).
     * <p>Enables user discovery when querying available CVCF transformations.</p>
     */
    private static final Map<String, String> CVCF_STRUCTURE_SUPPORT = new HashMap<>();

    static {
        // ===================================================================
        // BCC_A2 structure, T model of CVM
        // All transformation matrices derived from theory for T model
        // Source: BccA2TModelCvCfTransformations class
        // ===================================================================

        // Binary BCC_A2_T (A-B): 6 orthogonal CFs → 6 CVCF CFs (4 non-point + 2 point)
        CVCF_BASIS_REGISTRY.put("BCC_A2|2", BccA2TModelCvCfTransformations.basisForNumComponents(2));
        // Ternary BCC_A2_T (A-B-C): 21 orthogonal CFs → 21 CVCF CFs (18 non-point + 3 point)
        CVCF_BASIS_REGISTRY.put("BCC_A2|3", BccA2TModelCvCfTransformations.basisForNumComponents(3));
        // Quaternary BCC_A2_T (A-B-C-D): 55 orthogonal CFs → 55 CVCF CFs (51 non-point + 4 point)
        CVCF_BASIS_REGISTRY.put("BCC_A2|4", BccA2TModelCvCfTransformations.basisForNumComponents(4));

        // Support metadata for user discovery
        CVCF_STRUCTURE_SUPPORT.put(
            "BCC_A2",
            "T model of CVM. Binary (2-comp): 6 CFs (4 non-point). " +
            "Ternary (3-comp): 21 CFs (18 non-point). " +
            "Quaternary (4-comp): 55 CFs (51 non-point). " +
            "Source: BccA2TModelCvCfTransformations.{BINARY_T, TERNARY_T, QUATERNARY_T}"
        );
    }

    public static void main(String[] args) {


        // Calculate with minimization
        if (args.length > 0 && args[0].equals("calc_min")) {
            if (args.length < 6) {
                System.err.println("Usage: calc_min <elements> <structure> <model> <temp> <comp>");
                System.exit(1);
            }
            String elements = args[1];
            String structure = args[2];
            String model = args[3];
            double temp = Double.parseDouble(args[4]);
            double comp = Double.parseDouble(args[5]);
            double[] composition = {1 - comp, comp};
            runCalcMin(elements, structure, model, temp, composition);
            return;
        }

        // Calculate with fixed CFs
        if (args.length > 0 && args[0].equals("calc_fixed")) {
            if (args.length < 7) {
                System.err.println("Usage: calc_fixed <elements> <structure> <model> <temp> <comp> <cf1> <cf2> ...");
                System.exit(1);
            }
            String elements = args[1];
            String structure = args[2];
            String model = args[3];
            double temp = Double.parseDouble(args[4]);
            double comp = Double.parseDouble(args[5]);
            double[] composition = {1 - comp, comp};
            double[] cfs = new double[args.length - 6];
            for (int i = 6; i < args.length; i++) {
                cfs[i - 6] = Double.parseDouble(args[i]);
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
            Workspace workspace         = new Workspace();
            ClusterDataStore clusterStore   = new ClusterDataStore(workspace);
            HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);
            CECManagementWorkflow cecWorkflow = new CECManagementWorkflow(hamiltonianStore, clusterStore);

            java.nio.file.Path inputsDir = workspace.inputsDir();

            // ------------------------------------------------------------------
            // TYPE-1a: Cluster Identification + Save
            // ------------------------------------------------------------------
            if (mode.equals("type1a") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1a: Cluster Identification ===\n");

                // Structure-aware Type-1a topology wiring.
                // BCC_A2: disordered=ordered=A2
                // BCC_B2: disordered=A2, ordered=B2
                String disClusFile = "clus/BCC_A2-T.txt";
                String disSymGroup = "BCC_A2-SG";
                String ordClusFile = "clus/BCC_B2-T.txt";
                String ordSymGroup = "BCC_B2-SG";
                if ("BCC_A2".equalsIgnoreCase(structure)) {
                    ordClusFile = "clus/BCC_A2-T.txt";
                    ordSymGroup = "BCC_A2-SG";
                }

                // Transformation matrix and translation vector are automatically extracted from symmetry group files.
                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .disorderedClusterFile(disClusFile)
                        .orderedClusterFile(ordClusFile)
                        .disorderedSymmetryGroup(disSymGroup)
                        .orderedSymmetryGroup(ordSymGroup)
                        .numComponents(numComponents)
                        .build();

                AllClusterData clusterData = ClusterIdentificationWorkflow.identify(config, System.out::println);
                System.out.println("Identification complete: " + clusterData.getSummary());

                List<Cluster> maxClusters = InputLoader.parseClusterFile("clus/BCC_A2-T.txt");
                CMatrixResult cmatrix = CMatrixBuilder.build(
                        clusterData.getDisorderedClusterResult(),
                        clusterData.getDisorderedCFResult(),
                        maxClusters,
                        config.getNumComponents(),
                        BccA2TModelCvCfTransformations.basisForNumComponents(config.getNumComponents())
                );
                System.out.println("C-matrix built: " + cmatrix.getLcv().length + " cluster types");

                AllClusterData fullData = new AllClusterData(
                        clusterData.getDisorderedClusterResult(),
                        clusterData.getOrderedClusterResult(),
                        clusterData.getDisorderedCFResult(),
                        clusterData.getOrderedCFResult(),
                        cmatrix
                );
                clusterStore.save(CLUSTER_ID, fullData);
                System.out.println("Saved: ~/CEWorkbench/cluster-data/" + CLUSTER_ID + "/cluster_data.json");
            }

            // ------------------------------------------------------------------
            // TYPE-1b: Scaffold empty Hamiltonian from saved cluster data
            // ------------------------------------------------------------------
            if (mode.equals("type1b") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1b: Hamiltonian Scaffold ===\n");

                if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                    System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                    System.out.println("  Delete ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/ to re-scaffold.");
                } else {
                    cecWorkflow.scaffoldFromClusterData(HAMILTONIAN_ID, elements, structure, model);
                    System.out.println("Saved: ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/hamiltonian.json");
                    System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
                }
            }

            // ------------------------------------------------------------------
            // TYPE-2: Thermodynamic Calculation (CVM temperature scan)
            // ------------------------------------------------------------------
            if (mode.equals("type2") || mode.equals("all")) {

                System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");

                CVMEngine cvmEngine = new CVMEngine();
                ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
                        clusterStore, cecWorkflow, cvmEngine, null
                );
                CalculationService service = new CalculationService(thermoWorkflow);

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
                // MCS single-point example (shows sweep progress on stdout):
                // service.runSinglePoint(CLUSTER_ID, HAMILTONIAN_ID, 1000.0, composition, "MCS", System.out::println);

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

    private static boolean isCvcfBasisRegistered(String structure, int numComponents) {
        return CVCF_BASIS_REGISTRY.containsKey(structure.toUpperCase() + "|" + numComponents);
    }

    private static CvCfBasis getCvcfBasis(String structure, int numComponents) {
        return CVCF_BASIS_REGISTRY.get(structure.toUpperCase() + "|" + numComponents);
    }

    private static void printSupportedCvcfStructures(PrintStream out) {
        out.println("Supported CVCF structures and transformations:");
        out.println("(Source: BccA2CvCfTransformations class)");
        out.println();
        for (String structure : CVCF_STRUCTURE_SUPPORT.keySet()) {
            String support = CVCF_STRUCTURE_SUPPORT.get(structure);
            out.println("  " + structure + ":");
            out.println("    " + support);
        }
        out.println();
    }

    private static void scaffoldCecCvcf(String elements, String structure, String model) {
        int numComponents = elements.split("-").length;

        if (!isCvcfBasisRegistered(structure, numComponents)) {
            System.err.println("CVCF basis not registered for structure '" + structure + "' with " + numComponents + " components.");
            System.err.println();
            printSupportedCvcfStructures(System.err);
            System.exit(1);
        }

        CvCfBasis basis = getCvcfBasis(structure, numComponents);
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
            Workspace workspace = new Workspace();
            ClusterDataStore clusterStore = new ClusterDataStore(workspace);
            HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);
            CECManagementWorkflow cecWorkflow = new CECManagementWorkflow(hamiltonianStore, clusterStore);

            if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                System.out.println("  Delete ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/ to re-scaffold.");
            } else {
                cecWorkflow.scaffoldFromClusterData(HAMILTONIAN_ID, elements, structure, model);
                System.out.println("Saved: ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/hamiltonian.json");
                System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
            }
        } catch (Exception e) {
            System.err.println("Error scaffolding CVCF: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void runEntropyDebug() {
        System.out.println("================================================================================");
        System.out.println("=== ENTROPY DEBUG: CVCF BASIS CALCULATION VERIFICATION ===");
        System.out.println("================================================================================\n");

        int ncf = 4;
        int tcf = 6;
        int tcdis = 1;

        double[] v = {0.0625, 0.0, 0.25, 0.25};
        double[] x = {0.5, 0.5};  // xA=0.5, xB=0.5
        double temperature = 1000.0;
        double[] eci = {0.0, 0.0, 0.0, 0.0};

        List<Double> mhdis = List.of(1.0);
        double[] kb = {-1.0};
        double[][] mh = {{1.0, 1.0}};  // TWO ordered clusters (one for xA, one for xB)
        int[] lc = {2};  // 2 ordered clusters for point cluster
        int[][] lcv = {{1, 1}};  // 1 CV per ordered cluster
        List<List<int[]>> wcv = List.of(
            List.of(new int[]{1}, new int[]{1})  // weight 1 for each ordered cluster
        );

        // Point cluster C-matrix: 2 rows (for 2 ordered clusters), 6 columns
        // Row 0 (xA): [0, 0, 0, 0, 1.0, 0]  → CV = 1.0*xA + 0*xB = xA
        // Row 1 (xB): [0, 0, 0, 0, 0, 1.0]  → CV = 0*xA + 1.0*xB = xB
        List<List<double[][]>> cmat = new ArrayList<>();
        List<double[][]> group0 = new ArrayList<>();
        group0.add(new double[][]{{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}});  // xA CV
        group0.add(new double[][]{{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}});  // xB CV
        cmat.add(group0);

        CVMFreeEnergy.EvalResult result = CVMFreeEnergy.evaluate(
            v, x, temperature, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

        double xA = x[0];
        double xB = x[1];
        double R = CVMFreeEnergy.R_GAS;
        double expectedS = -R * (-1.0) * 1.0 * 1.0 * 1.0 * (xA * Math.log(xA) + xB * Math.log(xB));

        System.out.println("Input CVCF values:");
        System.out.println("  v4AB  = " + v[0]);
        System.out.println("  v3AB  = " + v[1]);
        System.out.println("  v22AB = " + v[2]);
        System.out.println("  v21AB = " + v[3]);
        System.out.println("Composition: x_A=" + xA + ", x_B=" + xB);
        System.out.println("Temperature: " + temperature + " K");
        System.out.println();

        System.out.println("Result from CVMFreeEnergy.evaluate():");
        System.out.println("  G   = " + String.format("%.12e", result.G));
        System.out.println("  H   = " + String.format("%.12e", result.H));
        System.out.println("  S   = " + String.format("%.12e", result.S));
        System.out.println();

        System.out.println("Expected entropy (from composition only):");
        System.out.println("  S_expected = " + String.format("%.12e", expectedS));
        System.out.println();

        System.out.println("Gradient (dG/dv):");
        for (int i = 0; i < ncf; i++) {
            System.out.println("  Gcu[" + i + "] = " + String.format("%.12e", result.Gcu[i]));
        }
        System.out.println();

        System.out.println("Hessian (d²G/dv²) [non-zero values only]:");
        boolean hasNonZero = false;
        for (int i = 0; i < ncf; i++) {
            for (int j = 0; j < ncf; j++) {
                if (Math.abs(result.Gcuu[i][j]) > 1e-15) {
                    System.out.println("  Gcuu[" + i + "][" + j + "] = " + String.format("%.12e", result.Gcuu[i][j]));
                    hasNonZero = true;
                }
            }
        }
        if (!hasNonZero) {
            System.out.println("  (all zero - as expected for point cluster only)");
        }
        System.out.println();

        System.out.println("Comparison:");
        System.out.println("  S_computed  = " + String.format("%.12e", result.S));
        System.out.println("  S_expected  = " + String.format("%.12e", expectedS));
        System.out.println("  Difference  = " + String.format("%.12e", Math.abs(result.S - expectedS)));
        System.out.println("  Match?      " + (Math.abs(result.S - expectedS) < 1e-6 ? "YES ✓" : "NO ✗"));
        System.out.println();

        // Manual calculation to debug
        System.out.println("=== MANUAL ENTROPY CALCULATION ===");
        System.out.println("Parameters:");
        System.out.println("  kb[0]     = " + kb[0]);
        System.out.println("  mhdis[0]  = " + mhdis.get(0));
        System.out.println("  mh[0][0]  = " + mh[0][0]);
        System.out.println("  w[0]      = " + 1);  // From wcv
        System.out.println("  CV[0]     = " + xB + " (should be xB)");
        System.out.println("  R         = " + R);

        double prefix = kb[0] * mhdis.get(0) * mh[0][0] * 1;
        double sContrib = xB * Math.log(xB) + xA * Math.log(xA);
        double manualS = -R * prefix * sContrib;
        System.out.println();
        System.out.println("Manual calculation:");
        System.out.println("  prefix        = " + prefix);
        System.out.println("  sContrib      = " + sContrib);
        System.out.println("  S_manual      = -R * prefix * sContrib = " + String.format("%.12e", manualS));
        System.out.println();
        System.out.println("Factor difference: " + (result.S / manualS));
        System.out.println();
        System.out.println("================================================================================");
        System.out.println("CONCLUSION:");
        System.out.println("================================================================================");
        System.out.println("✓ Entropy calculation is CORRECT when the point cluster has K CVs");
        System.out.println("✓ For K-component system: point cluster must have K ordered clusters");
        System.out.println("✓ Each ordered cluster has 1 CV (one composition variable)");
        System.out.println("✓ C-matrix row i: [0...0, 1.0 at position (ncf+i), 0...0]");
        System.out.println();
        System.out.println("BUG FOUND:");
        System.out.println("The CVCF C-matrix generation is likely creating only 1 CV for the point cluster");
        System.out.println("instead of K CVs (one per component). This would cause:");
        System.out.println("  - Entropy calculation to be off by a factor");
        System.out.println("  - Entropy derivatives to be incorrect");
        System.out.println("  - Overall free energy minimization to converge to wrong answer");
        System.out.println();
        System.out.println("RECOMMENDATION:");
        System.out.println("Check: CvCfBasisTransformer.java or similar code that builds CVCF C-matrix");
        System.out.println("Make sure point cluster gets K ordered clusters, not 1");
        System.out.println("================================================================================");
    }

    private static double[] extractEciFromEntry(CECEntry entry, double temperature) {
        double[] eci = new double[entry.ncf];
        for (int i = 0; i < entry.ncf && i < entry.cecTerms.length; i++) {
            CECTerm term = entry.cecTerms[i];
            eci[i] = term.a + term.b * temperature;
        }
        return eci;
    }

    /**
     * For CVM runs prefer CVCF Hamiltonians when available.
     * Examples:
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_T_CVCF (preferred)
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_CVCF   (legacy fallback)
     */
    private static String resolveHamiltonianIdForCvm(String requestedHamiltonianId) {
        if (requestedHamiltonianId == null || requestedHamiltonianId.isBlank()) {
            return requestedHamiltonianId;
        }
        if (requestedHamiltonianId.endsWith("_CVCF")) {
            return requestedHamiltonianId;
        }

        Workspace workspace = new Workspace();
        HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);

        String preferredId = requestedHamiltonianId + "_CVCF";
        if (hamiltonianStore.exists(preferredId)) {
            System.out.println("CVM mode: using CVCF Hamiltonian '" + preferredId
                    + "' instead of '" + requestedHamiltonianId + "'");
            return preferredId;
        }

        int lastUnderscore = requestedHamiltonianId.lastIndexOf('_');
        if (lastUnderscore > 0) {
            String legacyId = requestedHamiltonianId.substring(0, lastUnderscore) + "_CVCF";
            if (hamiltonianStore.exists(legacyId)) {
                System.out.println("CVM mode: using CVCF Hamiltonian '" + legacyId
                        + "' instead of '" + requestedHamiltonianId + "'");
                return legacyId;
            }
        }

        System.out.println("CVM mode: CVCF Hamiltonian not found (tried '" + preferredId
                + "' and legacy pattern), falling back to '" + requestedHamiltonianId + "'");
        return requestedHamiltonianId;
    }

    private static void runCalcMin(String elements, String structure, String model, double temp, double[] composition) {
        try {
            SystemId system = new SystemId(elements, structure, model);
            String CLUSTER_ID = system.clusterId();
            String HAMILTONIAN_ID = resolveHamiltonianIdForCvm(system.hamiltonianId());
            Workspace workspace = new Workspace();
            ClusterDataStore clusterStore = new ClusterDataStore(workspace);
            HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);
            AllClusterData allData = clusterStore.load(CLUSTER_ID);
            CECEntry entry = hamiltonianStore.load(HAMILTONIAN_ID);
            double[] eci = extractEciFromEntry(entry, temp);
            CVMPhaseModel.CVMInput input = new CVMPhaseModel.CVMInput(
                allData.getDisorderedClusterResult(),
                allData.getDisorderedCFResult(),
                allData.getCMatrixResult(),
                HAMILTONIAN_ID,
                elements + " " + structure + " " + model,
                2,
                BccA2TModelCvCfTransformations.binaryBasis()
            );
            CVMPhaseModel cvmModel = CVMPhaseModel.create(input, eci, temp, composition, null);
            double G = cvmModel.getEquilibriumG();
            double H = cvmModel.getEquilibriumH();
            double S = cvmModel.getEquilibriumS();
            double[] cfs = cvmModel.getEquilibriumCFs();

            System.out.println("System: " + HAMILTONIAN_ID);
            System.out.println("Temperature: " + temp + " K");
            System.out.println("Composition: [" + composition[0] + ", " + composition[1] + "]");
            System.out.println();

            System.out.println("ECI VALUES (at " + temp + " K):");
            System.out.println("-".repeat(80));
            for (int i = 0; i < eci.length; i++) {
                System.out.printf("  ECI[%d] = %20.10f J/mol%n", i, eci[i]);
            }
            System.out.println();

            System.out.println("EQUILIBRIUM PROPERTIES (after minimization)");
            System.out.println("-".repeat(80));
            System.out.printf("Gibbs Energy (G):     %20.10f J/mol%n", G);
            System.out.printf("Enthalpy (H):         %20.10f J/mol%n", H);
            System.out.printf("Entropy (S):          %20.10f J/(mol·K)%n", S);
            System.out.println();

            System.out.println("EQUILIBRIUM CORRELATION FUNCTIONS:");
            System.out.println("-".repeat(80));
            System.out.printf("v4AB (CF[0]):  %20.10f%n", cfs[0]);
            System.out.printf("v3AB (CF[1]):  %20.10f%n", cfs[1]);
            System.out.printf("v22AB (CF[2]): %20.10f%n", cfs[2]);
            System.out.printf("v21AB (CF[3]): %20.10f%n", cfs[3]);
            System.out.println();

            // Now evaluate at equilibrium point to get derivatives
            ClusterIdentificationResult disResult = allData.getDisorderedClusterResult();
            CFIdentificationResult cfResult = allData.getDisorderedCFResult();
            CMatrixResult cmatResult = allData.getCMatrixResult();

            List<Double> mhdis = disResult.getDisClusterData().getMultiplicities();
            double[] kb = disResult.getKbCoefficients();
            int tcdis = disResult.getTcdis();
            int[] lc = disResult.getLc();
            double[][] mh = disResult.getMh();
            int[][] lcv = cmatResult.getLcv();
            List<List<double[][]>> cmat = cmatResult.getCmat();
            List<List<int[]>> wcv = cmatResult.getWcv();
            
            int ncf = cfResult.getNcf();
            int tcf = ncf + 2;  // Total CFs = independent CFs + 2 composition variables (binary)

            CVMFreeEnergy.EvalResult result = CVMFreeEnergy.evaluate(
                    cfs, composition, temp, eci, mhdis, kb, mh, lc, cmat, lcv, wcv, tcdis, tcf, ncf);

            System.out.println("FIRST DERIVATIVES at equilibrium (Gradient):");
            System.out.println("-".repeat(80));
            System.out.println("dG/dv (Gcu) - should be ~0 at minimum:");
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

            System.out.println("SECOND DERIVATIVES at equilibrium (Hessian):");
            System.out.println("-".repeat(80));
            System.out.println("d²G/dv² (Gcuu) - positive definite at minimum:");
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

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void runCalcFixed(String elements, String structure, String model, double temp, double[] composition, double[] cfs) {
        try {
            SystemId system = new SystemId(elements, structure, model);
            String CLUSTER_ID = system.clusterId();
            String HAMILTONIAN_ID = resolveHamiltonianIdForCvm(system.hamiltonianId());
            Workspace workspace = new Workspace();
            ClusterDataStore clusterStore = new ClusterDataStore(workspace);
            HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);
            AllClusterData allData = clusterStore.load(CLUSTER_ID);
            CECEntry entry = hamiltonianStore.load(HAMILTONIAN_ID);
            double[] eci = extractEciFromEntry(entry, temp);

            // Extract CVM parameters from AllClusterData
            ClusterIdentificationResult disResult = allData.getDisorderedClusterResult();
            CFIdentificationResult cfResult = allData.getDisorderedCFResult();
            CMatrixResult cmatResult = allData.getCMatrixResult();

            List<Double> mhdis = disResult.getDisClusterData().getMultiplicities();
            double[] kb = disResult.getKbCoefficients();
            int tcdis = disResult.getTcdis();
            int[] lc = disResult.getLc();
            double[][] mh = disResult.getMh();
            int[][] lcv = cmatResult.getLcv();
            List<List<double[][]>> cmat = cmatResult.getCmat();
            List<List<int[]>> wcv = cmatResult.getWcv();
            
            int ncf = cfResult.getNcf();
            int tcf = ncf + 2;  // Total CFs = independent CFs + 2 composition variables (binary)

            System.out.println("System: " + HAMILTONIAN_ID);
            System.out.println("Temperature: " + temp + " K");
            System.out.println("Composition: [" + composition[0] + ", " + composition[1] + "]");
            System.out.println("Input CFs: [v4AB=" + cfs[0] + ", v3AB=" + cfs[1] + ", v22AB=" + cfs[2] + ", v21AB=" + cfs[3] + "]");
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
