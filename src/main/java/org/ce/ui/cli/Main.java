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
import org.ce.calculation.Conditions;
import org.ce.calculation.ConditionsScan;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.Range;
import org.ce.calculation.ResultFormatter;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.cluster.ClusterIdentificationRequest;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
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

        // ── api ──────────────────────────────────────────────────────────────
        // JSON request on stdin -> JSON response on stdout. Diagnostics go to
        // stderr so the payload stays machine-parseable.
        if (args.length > 0 && args[0].equals("api")) {
            System.exit(ApiCommand.run(appCtx, System.in));
            return;
        }

        // ── calc_min ─────────────────────────────────────────────────────────
        if (args.length > 0 && args[0].equals("calc_min")) {
            if (args.length < 6) {
                System.err.println(
                        "Usage: calc_min <elements> <structure> <model> <temp> <El>=<x> [<El>=<x> ...] [G|H|S] [--verbose]\n" +
                        "  e.g.  calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5\n" +
                        "  e.g.  calc_min Nb-Ti-V-Zr BCC_A2 T 1273 Ti=0.25 V=0.25 Zr=0.25 S\n" +
                        "  Any element may be omitted; its fraction is derived as 1 - sum(given).");
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

            Map<String, Double> comp = new LinkedHashMap<>();
            for (int i = 5; i < compEndIndex; i++) {
                Map.Entry<String, Range> tok = parseCompToken(args[i]);
                if (tok == null) {
                    System.err.println(
                            "Error: positional mole fractions are no longer supported.\n" +
                            "  Old: calc_min Nb-Ti-V BCC_A2 T 1000 0.3 0.4\n" +
                            "  New: calc_min Nb-Ti-V BCC_A2 T 1000 Ti=0.3 V=0.4");
                    System.exit(1);
                }
                comp.put(tok.getKey(), tok.getValue().start());
            }
            runCalcMin(appCtx, elements, structure, model, temp, comp, requestedProp);
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

                // Equiatomic composition: 1/K for elements 2..K, element 1 derived.
                List<String> elemList = system.elementList();
                Map<String, Double> comp = new LinkedHashMap<>();
                double x = 1.0 / numComponents;
                for (int i = 1; i < elemList.size(); i++) comp.put(elemList.get(i), x);

                double tStart = 1000.0, tEnd = 1000.0, tStep = 100.0;

                ModelSpecifications modelSpecs = new ModelSpecifications(elements, structure, model, EngineConfig.CVM);
                ModelSession session = service.getOrBuildSession(modelSpecs, sink);
                ConditionsScan scan = new ConditionsScan(new Range(tStart, tEnd, tStep), toRanges(comp));

                if (verbose) {
                    System.out.println("System      : " + modelSpecs);
                    System.out.println("Composition : " + comp);
                    System.out.println("T range     : " + tStart + " K to " + tEnd + " K, step " + tStep + " K\n");
                }

                List<ThermodynamicResult> results =
                        service.calculateScan(session, scan, Property.GIBBS_ENERGY, sink, null);
                printResult(new CalculationResult.Grid(List.of(results)));
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
                                   double temp, Map<String, Double> composition, Property requestedProp) {
        try {
            CalculationService service = appCtx.getCalculationService();
            Consumer<String> sink = verbose ? System.out::println : null;

            ModelSpecifications modelSpecs = new ModelSpecifications(elements, structure, model, EngineConfig.CVM);
            ModelSession session = service.getOrBuildSession(modelSpecs, sink);
            Conditions conditions = new Conditions(temp, composition);

            System.out.println("System: " + modelSpecs);
            System.out.println();
            Consumer<org.ce.model.ProgressEvent> eventSink = verbose ? ev -> {
                if (ev instanceof org.ce.model.ProgressEvent.CvmIteration it) {
                    System.out.printf("  [ITER %d] G=%.6f H=%.6f S=%.6f gradNorm=%.6e cfs=%s%n",
                            it.iteration, it.gibbsEnergy, it.enthalpy, it.entropy, it.gradientNorm,
                            java.util.Arrays.toString(it.cfs));
                }
            } : null;
            ThermodynamicResult result = service.calculate(session, conditions, requestedProp, sink, eventSink);
            printResult(new CalculationResult.Single(result));

        } catch (Exception e) {
            if (verbose) e.printStackTrace();
            else System.err.println("Error: " + e.getMessage());
        }
    }

    /** Parses "Ti=0.3" or "Ti=0.1:0.9:0.1". Returns null if not a composition token. */
    private static Map.Entry<String, Range> parseCompToken(String tok) {
        int eq = tok.indexOf('=');
        if (eq <= 0) return null;
        String sym = tok.substring(0, eq).trim();
        String val = tok.substring(eq + 1).trim();
        String[] parts = val.split(":");
        try {
            return switch (parts.length) {
                case 1 -> Map.entry(sym, Range.fixed(Double.parseDouble(parts[0])));
                case 3 -> Map.entry(sym, new Range(Double.parseDouble(parts[0]),
                                                    Double.parseDouble(parts[1]),
                                                    Double.parseDouble(parts[2])));
                default -> throw new IllegalArgumentException(
                        "Bad composition token '" + tok + "'. Use El=x or El=start:end:step");
            };
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Non-numeric mole fraction in '" + tok + "'");
        }
    }

    private static Map<String, Range> toRanges(Map<String, Double> fixed) {
        Map<String, Range> out = new LinkedHashMap<>();
        for (var e : fixed.entrySet()) out.put(e.getKey(), Range.fixed(e.getValue()));
        return out;
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
