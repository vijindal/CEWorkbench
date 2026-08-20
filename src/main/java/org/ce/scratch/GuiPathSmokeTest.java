package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.JobSpecifications;
import org.ce.calculation.CalculationDescriptor.Mode;
import org.ce.calculation.CalculationDescriptor.ModelSpecifications;
import org.ce.calculation.CalculationDescriptor.Parameter;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.CalculationResult;
import org.ce.calculation.ConditionsScan;
import org.ce.calculation.Range;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.ModelSession;
import org.ce.model.ThermodynamicResult;
import org.ce.model.storage.Workspace;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Exercises the exact call the GUI makes, without opening a window.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.GuiPathSmokeTest
 * </pre>
 *
 * <p>{@code DynamicCalculationPanel.startExecution} builds a
 * {@link ModelSpecifications} and a {@link JobSpecifications}, sets a
 * {@link ConditionsScan}, and calls {@link CalculationService#execute} on a
 * {@code SwingWorker}. That {@code execute} overload is the one calculation
 * path the CLI and the JSON API never touch -- they call {@code calculate} /
 * {@code calculateScan} directly -- so a break in it would not show up in any
 * other gate.</p>
 *
 * <p>Everything below the {@code SwingWorker} is plain Java, so running it on
 * the main thread reproduces the GUI's behaviour exactly. What this does
 * <em>not</em> cover is Swing itself: widget layout, the Nimbus dark theme, and
 * the event-dispatch threading rules.</p>
 */
public final class GuiPathSmokeTest {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(72));
        System.out.println("  GUI calculation path: CalculationService.execute");
        System.out.println("=".repeat(72));

        CEWorkbenchContext context = new CEWorkbenchContext(new Workspace());
        CalculationService service = context.getCalculationService();

        // The GUI's own defaults, as DynamicCalculationPanel pre-fills them.
        ModelSpecifications modelSpecs = new ModelSpecifications(
                "Nb-Ti", "BCC_A2", "T", ModelSession.EngineConfig.CVM);

        JobSpecifications specs = new JobSpecifications(Property.GIBBS_ENERGY, Mode.ANALYSIS);

        // A short temperature scan at fixed composition -- the shape the panel
        // builds whenever composition spinners are present.
        Map<String, Range> ranges = new LinkedHashMap<>();
        ranges.put("Ti", Range.fixed(0.5));
        specs.set(Parameter.CONDITIONS_SCAN,
                new ConditionsScan(new Range(900, 1100, 100), ranges));

        StringBuilder log = new StringBuilder();
        CalculationResult result = service.execute(
                modelSpecs, specs, line -> log.append(line).append('\n'), null);

        check("result is not null", result != null);
        if (result == null) {
            finish();
            return;
        }

        List<ThermodynamicResult> points = flatten(result);
        System.out.printf("%n  scan returned %d point(s)%n", points.size());
        check("scan produced 3 points (900, 1000, 1100)", points.size() == 3);

        boolean sawBaseline = false;
        for (ThermodynamicResult p : points) {
            System.out.printf("    T=%.0f  G=%.10f  converged=%s%n",
                    p.temperature, p.gibbsEnergy, p.converged);
            check(String.format("T=%.0f converged", p.temperature),
                    Boolean.TRUE.equals(p.converged));

            // The CLAUDE.md baseline sits inside this scan. If the GUI path
            // reaches a different number than the CLI does, they have diverged.
            if (Math.abs(p.temperature - 1000.0) < 1e-9) {
                sawBaseline = true;
                check("T=1000 matches the CLI baseline -3480.5209063902",
                        Math.abs(p.gibbsEnergy - (-3480.5209063902)) < 1e-9);
            }
        }
        check("scan included T=1000", sawBaseline);

        // The panel streams these lines into OutputPanel; empty would mean the
        // GUI shows a blank log even on a successful run.
        check("progress sink received output", log.length() > 0);
        System.out.printf("  progress sink: %d character(s)%n", log.length());

        finish();
    }

    /** Both CalculationResult variants, flattened to a point list. */
    private static List<ThermodynamicResult> flatten(CalculationResult r) {
        List<ThermodynamicResult> out = new ArrayList<>();
        if (r instanceof CalculationResult.Single s) {
            out.add(s.value());
        } else if (r instanceof CalculationResult.Grid g) {
            for (List<ThermodynamicResult> row : g.values()) {
                out.addAll(row);
            }
        }
        return out;
    }

    private static void finish() {
        System.out.println("\n" + "=".repeat(72));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(72));
        if (failures > 0) {
            throw new AssertionError(failures + " GUI-path checks failed");
        }
    }

    private static void check(String what, boolean ok) {
        if (!ok) {
            failures++;
            System.out.println("    [!] FAIL  " + what);
        }
    }

    private GuiPathSmokeTest() {
    }
}
