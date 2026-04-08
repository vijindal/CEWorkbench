package org.ce.ui.cli;

import org.ce.domain.cluster.*;
import org.ce.domain.engine.cvm.CVMGibbsModel;
import org.ce.domain.engine.cvm.CVMGibbsModel.ModelResult;
import org.ce.storage.Workspace.SystemId;
import org.ce.storage.HamiltonianStore;
import org.ce.workflow.CalculationService;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.workflow.thermo.ThermodynamicRequest;

import java.util.List;
import java.util.function.Consumer;

/**
 * Main entry point for the CE Workbench CLI.
 */
public class Main {

    private static boolean verbose = false;

    private static void setupLogging() {
        if (!verbose) {
            java.util.logging.LogManager.getLogManager().reset();
            java.util.logging.Logger rootLogger = java.util.logging.Logger.getLogger("");
            rootLogger.setLevel(java.util.logging.Level.OFF);
        }
    }

    public static void main(String[] args) {
        // Extract flags
        java.util.List<String> argList = new java.util.ArrayList<>(java.util.Arrays.asList(args));
        verbose = argList.remove("--verbose") || argList.remove("-v");
        args = argList.toArray(new String[0]);

        setupLogging();

        // Calculate with minimization
        if (args.length > 0 && args[0].equals("calc_min")) {
            if (args.length < 6) {
                System.err
                        .println("Usage: calc_min <elements> <structure> <model> <temp> <comp1> [<comp2> <comp3> ...] [--verbose]");
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
                composition = new double[] { 1 - comp, comp };
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
                System.err.println(
                        "Usage: calc_fixed <elements> <structure> <model> <temp> <comp1> [<comp2> ...compK] <cf1> <cf2> ... [--verbose]");
                System.err.println("  For a " + K + "-component system provide " + K
                        + " composition values followed by CF values.");
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
        String mode = args.length > 0 ? args[0].toLowerCase() : "all";
        String elements = args.length > 1 ? args[1] : "Nb-Ti";
        String structure = args.length > 2 ? args[2] : "BCC_B2";
        String model = args.length > 3 ? args[3] : "T";

        if (!mode.equals("type1a") && !mode.equals("type1b")
                && !mode.equals("type2") && !mode.equals("all")
                && !mode.equals("calc_min") && !mode.equals("calc_fixed")
                && !mode.equals("view")) {
            System.err.println("Unknown mode: " + mode);
            System.err.println("Usage: <mode> [elements] [structure] [model] [--verbose]");
            System.err.println("  mode: type1a | type1b | type2 | all | calc_min | calc_fixed | view");
            System.exit(1);
        }

        if (mode.equals("view")) {
            viewHamiltonian(elements, structure, model);
            return;
        }

        SystemId system = new SystemId(elements, structure, model);
        String HAMILTONIAN_ID = system.hamiltonianId();
        int numComponents = elements.split("-").length;

        System.out.println("================================================================================");
        System.out.println("                    CE THERMODYNAMICS WORKBENCH");
        System.out.println("================================================================================");
        System.out.println("Mode      : " + mode);
        System.out.println("System    : " + elements + "  /  " + structure + "  /  " + model);
        System.out.println("Hamiltonian ID: " + HAMILTONIAN_ID);

        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();
            CalculationService service = appCtx.getCalculationService();

            Consumer<String> sink = verbose ? System.out::println : null;

            // ------------------------------------------------------------------
            // TYPE-1a: Cluster Identification
            // ------------------------------------------------------------------
            if (mode.equals("type1a") || mode.equals("all")) {

                if (verbose) {
                    System.out.println("\n=== TYPE-1a: Cluster Identification ===\n");
                }

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
                        .model(model)
                        .build();

                appCtx.identifyClusters(config, sink);
                
                if (verbose) {
                    System.out.println("Identification complete (on-the-fly, not saved to disk).");
                }
            }

            // ------------------------------------------------------------------
            // TYPE-1b: Scaffold empty Hamiltonian
            // ------------------------------------------------------------------
            if (mode.equals("type1b") || mode.equals("all")) {

                if (verbose) {
                    System.out.println("\n=== TYPE-1b: Hamiltonian Scaffold ===\n");
                }

                if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                    if (verbose) {
                        System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                        System.out.println(
                                "  Delete " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID) + " to re-scaffold.");
                    }
                } else {
                    appCtx.getCecWorkflow().scaffoldFromClusterData(HAMILTONIAN_ID, elements, structure, model);
                    if (verbose) {
                        System.out.println("Saved: " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID));
                        System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
                    } else {
                        System.out.println("Hamiltonian scaffolded.");
                    }
                }
            }

            // ------------------------------------------------------------------
            // TYPE-2: Thermodynamic Calculation (CVM temperature scan)
            // ------------------------------------------------------------------
            if (mode.equals("type2") || mode.equals("all")) {

                if (verbose) {
                    System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");
                }

                double[] composition = { 0.5, 0.5 };
                if (numComponents > 2) {
                    composition = new double[numComponents];
                    double x = 1.0 / numComponents;
                    for (int i = 0; i < numComponents; i++) {
                        composition[i] = x;
                    }
                }
                double tStart = 1000.0;
                double tEnd = 1000.0;
                double tStep = 100.0;

                if (verbose) {
                    System.out.println("System      : " + HAMILTONIAN_ID);
                    System.out.println("Composition : " + java.util.Arrays.toString(composition));
                    System.out.println("T range     : " + tStart + " K to " + tEnd + " K, step " + tStep + " K\n");
                }

                List<ThermodynamicResult> results = service.runLineScanTemperature(
                        null, HAMILTONIAN_ID, composition, tStart, tEnd, tStep, "CVM");

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
            if (verbose) {
                e.printStackTrace();
            }
            System.exit(1);
        }

        System.out.println("\n================================================================================");
    }

    private static void viewHamiltonian(String elements, String structure, String model) {
        SystemId system = new SystemId(elements, structure, model);
        String hamiltonianId = system.hamiltonianId();
        
        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            CECEntry entry = appCtx.getCecWorkflow().loadAndValidateCEC(null, hamiltonianId);
            
            System.out.println("\n=== HAMILTONIAN: " + hamiltonianId + " ===");
            System.out.println("Elements: " + entry.elements + " | Structure: " + entry.structurePhase + " | Model: " + entry.model);
            System.out.println(String.format("\n  %-4s  %-12s  %-14s  %-14s", "Idx", "Name", "a (J/mol)", "b (J/mol/K)"));
            System.out.println("  " + "-".repeat(50));
            for (int i = 0; i < entry.cecTerms.length; i++) {
                CECTerm term = entry.cecTerms[i];
                System.out.println(String.format("  [%02d]  %-12s  %14.6f  %14.6f", i, term.name, term.a, term.b));
            }
            System.out.println();
        } catch (Exception e) {
            System.err.println("Error loading Hamiltonian: " + e.getMessage());
            if (verbose) e.printStackTrace();
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
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            CalculationService service = appCtx.getCalculationService();

            SystemId system = new SystemId(elements, structure, model);
            String hamiltonianId = system.hamiltonianId();

            ThermodynamicResult result = service.runSinglePoint(
                    new ThermodynamicRequest(null, hamiltonianId, temp, composition, "CVM",
                            verbose ? System.out::println : null, null));

            System.out.println("System: " + hamiltonianId);
            System.out.println("Temperature: " + temp + " K");
            System.out.print("Composition: [");
            for (int i = 0; i < composition.length; i++) {
                System.out.printf("%.4f", composition[i]);
                if (i < composition.length - 1)
                    System.out.print(", ");
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
            if (verbose) {
                e.printStackTrace();
            } else {
                System.err.println("Error: " + e.getMessage());
            }
        }
    }

    private static void runCalcFixed(String elements, String structure, String model, double temp, double[] composition,
            double[] cfs) {
        try {
            org.ce.CEWorkbenchContext appCtx = new org.ce.CEWorkbenchContext();
            SystemId system = new SystemId(elements, structure, model);
            String rawId = system.hamiltonianId();
            String preferredId = rawId + "_CVCF";

            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();
            String HAMILTONIAN_ID = hamiltonianStore.exists(preferredId) ? preferredId : rawId;

            // Run fresh identification
            ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                    .numComponents(composition.length)
                    .structurePhase(structure)
                    .model(model)
                    .build();
            AllClusterData allData = appCtx.identifyClusters(config, null);
            
            CECEntry entry = hamiltonianStore.load(HAMILTONIAN_ID);
            double[] eci = extractEciFromEntry(entry, temp);

            ClusterIdentificationResult disResult = allData.getDisorderedClusterResult();
            CFIdentificationResult cfResult = allData.getDisorderedCFResult();
            CMatrix.Result cmatResult = allData.getCMatrixResult();

            java.util.List<Double> mhdis = disResult.getDisClusterData().getMultiplicities();
            double[] kb = disResult.getKbCoefficients();
            int tcdis = disResult.getTcdis();
            int[] lc = disResult.getLc();
            double[][] mh = disResult.getMh();
            int[][] lcv = cmatResult.getLcv();
            java.util.List<java.util.List<double[][]>> cmat = cmatResult.getCmat();
            java.util.List<java.util.List<int[]>> wcv = cmatResult.getWcv();

            int ncf = cfResult.getNcf();
            int numComponents = composition.length;
            int tcf = ncf + numComponents;

            System.out.println("System: " + HAMILTONIAN_ID);
            System.out.println("Temperature: " + temp + " K");
            System.out.println("Composition: " + java.util.Arrays.toString(composition));
            
            ModelResult result = CVMGibbsModel.evaluate(
                    cfs, composition, temp, eci, tcdis, tcf, ncf, mhdis, kb, mh, lc, cmat, lcv, wcv);

            System.out.println("THERMODYNAMIC PROPERTIES:");
            System.out.println("-".repeat(80));
            System.out.printf("Gibbs Energy (G):     %20.10f J/mol%n", result.G);
            System.out.printf("Enthalpy (H):         %20.10f J/mol%n", result.H);
            System.out.printf("Entropy (S):          %20.10f J/(mol·K)%n", result.S);
            System.out.println();

        } catch (Exception e) {
            if (verbose) {
                e.printStackTrace();
            } else {
                System.err.println("Error: " + e.getMessage());
            }
        }
    }
}
