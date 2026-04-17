package org.ce.ui.cli;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.model.storage.DataStore.HamiltonianStore;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.ThermodynamicResult;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationResult;
import org.ce.calculation.CalculationSpecifications;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.ResultFormatter;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.cluster.ClusterIdentificationRequest;

import java.util.List;
import java.util.function.Consumer;

/**
 * Main entry point for the CE Workbench CLI.
 *
 * <p>Result formatting delegates to {@link ResultFormatter} (shared with GUI).
 * Quantity availability is determined via {@link QuantityDescriptor} (shared with GUI).</p>
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
        java.util.List<String> argList = new java.util.ArrayList<>(java.util.Arrays.asList(args));
        verbose = argList.remove("--verbose") || argList.remove("-v");
        args = argList.toArray(new String[0]);

        setupLogging();

        // Single shared context for all sub-commands
        CEWorkbenchContext appCtx;
        try {
            appCtx = new CEWorkbenchContext();
        } catch (Exception e) {
            System.err.println("Failed to initialise context: " + e.getMessage());
            System.exit(1);
            return;
        }

        // ── calc_min ─────────────────────────────────────────────────────────
        if (args.length > 0 && args[0].equals("calc_min")) {
            if (args.length < 6) {
                System.err.println(
                        "Usage: calc_min <elements> <structure> <model> <temp> <comp1> [<comp2> ...] [prop (G|H|S)] [--verbose]");
                System.exit(1);
            }
            String elements  = args[1];
            String structure = args[2];
            String model     = args[3];
            double temp      = Double.parseDouble(args[4]);

            String lastArg = args[args.length - 1];
            Property requestedProp = Property.GIBBS_ENERGY;
            int compEndIndex = args.length;

            if (lastArg.equalsIgnoreCase("G") || lastArg.equalsIgnoreCase("H") || lastArg.equalsIgnoreCase("S")) {
                requestedProp = switch (lastArg.toUpperCase()) {
                    case "H" -> Property.ENTHALPY;
                    case "S" -> Property.ENTROPY;
                    default  -> Property.GIBBS_ENERGY;
                };
                compEndIndex--;
            }

            double[] composition;
            if (compEndIndex == 6) {
                double comp = Double.parseDouble(args[5]);
                composition = new double[]{1 - comp, comp};
            } else {
                composition = new double[compEndIndex - 5];
                for (int i = 5; i < compEndIndex; i++) composition[i - 5] = Double.parseDouble(args[i]);
            }
            runCalcMin(appCtx, elements, structure, model, temp, composition, requestedProp);
            return;
        }

        // ── positional-arg modes ─────────────────────────────────────────────
        String mode      = args.length > 0 ? args[0].toLowerCase() : "all";
        String elements  = args.length > 1 ? args[1] : "Nb-Ti";
        String structure = args.length > 2 ? args[2] : "BCC_B2";
        String model     = args.length > 3 ? args[3] : "T";

        if (!mode.equals("type1a") && !mode.equals("type1b")
                && !mode.equals("type2") && !mode.equals("all")
                && !mode.equals("view")) {
            System.err.println("Unknown mode: " + mode);
            System.err.println("Usage: <mode> [elements] [structure] [model] [--verbose]");
            System.err.println("  mode: type1a | type1b | type2 | all | calc_min | view");
            System.exit(1);
        }

        if (mode.equals("view")) {
            viewHamiltonian(appCtx, elements, structure, model);
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
            HamiltonianStore hamiltonianStore = appCtx.getHamiltonianStore();
            CalculationService service = appCtx.getCalculationService();
            Consumer<String> sink = verbose ? System.out::println : null;

            // ── TYPE-1a: Cluster Identification ──────────────────────────────
            if (mode.equals("type1a") || mode.equals("all")) {
                if (verbose) System.out.println("\n=== TYPE-1a: Cluster Identification ===\n");

                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .numComponents(numComponents)
                        .structurePhase(structure)
                        .model(model)
                        .build();

                org.ce.model.cluster.ClusterCFIdentificationPipeline.runFullWorkflow(config, sink);
                if (verbose) System.out.println("Identification complete.");
            }

            // ── TYPE-1b: Scaffold empty Hamiltonian ───────────────────────────
            if (mode.equals("type1b") || mode.equals("all")) {
                if (verbose) System.out.println("\n=== TYPE-1b: Hamiltonian Scaffold ===\n");

                if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                    if (verbose) {
                        System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                        System.out.println("  Delete " + appCtx.getWorkspace().hamiltonianFile(HAMILTONIAN_ID) + " to re-scaffold.");
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

            // ── TYPE-2: Thermodynamic Calculation (CVM temperature scan) ──────
            if (mode.equals("type2") || mode.equals("all")) {
                if (verbose) System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");

                double[] composition = new double[numComponents];
                double x = 1.0 / numComponents;
                for (int i = 0; i < numComponents; i++) composition[i] = x;

                double tStart = 1000.0, tEnd = 1000.0, tStep = 100.0;

                ModelSpecifications modelSpecs = new ModelSpecifications(elements, structure, model, EngineConfig.CVM);
                CalculationSpecifications calcSpecs = new CalculationSpecifications(Property.GIBBS_ENERGY, Mode.ANALYSIS);
                calcSpecs.set(Parameter.COMPOSITION, composition);
                calcSpecs.set(Parameter.T_START, tStart);
                calcSpecs.set(Parameter.T_END, tEnd);
                calcSpecs.set(Parameter.T_STEP, tStep);

                if (verbose) {
                    System.out.println("System      : " + modelSpecs);
                    System.out.println("Composition : " + java.util.Arrays.toString(composition));
                    System.out.println("T range     : " + tStart + " K to " + tEnd + " K, step " + tStep + " K\n");
                }

                printResult(service.execute(modelSpecs, calcSpecs, sink, null));
            }

        } catch (Exception e) {
            System.err.println("\nError: " + e.getMessage());
            if (verbose) e.printStackTrace();
            System.exit(1);
        }

        System.out.println("\n================================================================================");
    }

    // =========================================================================
    // Sub-commands
    // =========================================================================

    private static void viewHamiltonian(CEWorkbenchContext appCtx, String elements, String structure, String model) {
        SystemId system = new SystemId(elements, structure, model);
        String hamiltonianId = system.hamiltonianId();

        try {
            CECEntry entry = appCtx.getCecWorkflow().loadAndValidateCEC(null, hamiltonianId);

            System.out.println("\n=== HAMILTONIAN: " + hamiltonianId + " ===");
            System.out.println("Elements: " + entry.elements
                    + " | Structure: " + entry.structurePhase
                    + " | Model: " + entry.model);
            System.out.println(String.format("\n  %-4s  %-12s  %-14s  %-14s",
                    "Idx", "Name", "a (J/mol)", "b (J/mol/K)"));
            System.out.println("  " + "-".repeat(50));
            for (int i = 0; i < entry.cecTerms.length; i++) {
                CECEntry.CECTerm term = entry.cecTerms[i];
                System.out.println(String.format("  [%02d]  %-12s  %14.6f  %14.6f",
                        i, term.name, term.a, term.b));
            }
            System.out.println();
        } catch (Exception e) {
            System.err.println("Error loading Hamiltonian: " + e.getMessage());
            if (verbose) e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Single-point calculation at given temperature and composition.
     */
    private static void runCalcMin(CEWorkbenchContext appCtx,
                                   String elements, String structure, String model,
                                   double temp, double[] composition, Property requestedProp) {
        try {
            CalculationService service = appCtx.getCalculationService();
            Consumer<String> sink = verbose ? System.out::println : null;

            ModelSpecifications modelSpecs = new ModelSpecifications(elements, structure, model, EngineConfig.CVM);
            CalculationSpecifications calcSpecs = new CalculationSpecifications(requestedProp, Mode.ANALYSIS);
            calcSpecs.set(Parameter.T_START, temp);
            calcSpecs.set(Parameter.T_END,   temp);
            calcSpecs.set(Parameter.T_STEP,  0.0);

            // X_STARTS holds the independent fractions (x2, x3, ...) — x1 is derived as 1-sum
            double[] xIndep = java.util.Arrays.copyOfRange(composition, 1, composition.length);
            calcSpecs.set(Parameter.X_STARTS, xIndep);
            calcSpecs.set(Parameter.X_ENDS,   xIndep);
            calcSpecs.set(Parameter.X_STEPS,  new double[xIndep.length]);

            System.out.println("System: " + modelSpecs);
            System.out.println();
            printResult(service.execute(modelSpecs, calcSpecs, sink, null));

        } catch (Exception e) {
            if (verbose) e.printStackTrace();
            else System.err.println("Error: " + e.getMessage());
        }
    }

    // =========================================================================
    // Output helpers
    // =========================================================================

    private static void printResult(CalculationResult calcResult) {
        switch (calcResult) {
            case CalculationResult.Single s ->
                System.out.print(ResultFormatter.fullBlock(s.value()));
            case CalculationResult.Grid g -> {
                List<List<ThermodynamicResult>> grid = g.values();
                List<ThermodynamicResult> flat = grid.stream().flatMap(List::stream).toList();
                System.out.print(flat.size() == 1
                        ? ResultFormatter.fullBlock(flat.get(0))
                        : ResultFormatter.table(flat));
            }
        }
    }
}
