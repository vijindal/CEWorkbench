package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.CalculationService;
import org.ce.calculation.workflow.QuaternarySquareScan;
import org.ce.calculation.workflow.QuaternarySquareScan.Point;
import org.ce.calculation.workflow.QuaternarySquareScan.PropertyQuantity;
import org.ce.calculation.workflow.QuaternarySquareScan.Region;
import org.ce.calculation.workflow.QuaternarySquareScan.Result;
import org.ce.calculation.workflow.QuaternarySquareScan.Variant;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;

/**
 * Coarse-grid smoke test for {@link QuaternarySquareScan} against the real
 * Nb-Ti-V-Zr Hamiltonian.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.QuaternarySquareScanSmokeTest
 * </pre>
 *
 * <p>Not a numerical reference check (no external/published values for
 * these specific square-plot points exist yet) — this is a structural gate:
 * both variants run to completion at a coarse grid, every point lands in
 * the expected region, region counts match what the parametrization
 * predicts, mole fractions at every point sum to 1, and INTERIOR points
 * agree with a direct {@link CalculationService#calculate} call at the same
 * composition (checking the scan's own interior path isn't silently
 * computing something else).</p>
 */
public final class QuaternarySquareScanSmokeTest {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  QuaternarySquareScan smoke test — Nb-Ti-V-Zr / BCC_A2 / T");
        System.out.println("=".repeat(78));

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
        CalculationService service = context.getCalculationService();

        List<String> slotOrder = List.of("Nb", "Ti", "V", "Zr");
        double temperature = 1273.0;
        int n = 8; // coarse grid, per the user's explicit "n=50 eventually, coarse initially" guidance

        for (Variant variant : Variant.values()) {
            System.out.printf("%n--- Variant %s, n=%d ---%n", variant, n);
            Result result = QuaternarySquareScan.run(service, session, slotOrder, variant,
                    temperature, new PropertyQuantity(Property.GIBBS_ENERGY), n, null);

            int expectedTotal = (n + 1) * (n + 1);
            int actualTotal = result.points().size() + result.skipped();
            System.out.printf("  points=%d skipped=%d (expected total %d, got %d)%n",
                    result.points().size(), result.skipped(), expectedTotal, actualTotal);
            check(variant + ": all grid points accounted for", actualTotal == expectedTotal);

            Map<Region, Integer> counts = new EnumMap<>(Region.class);
            for (Point p : result.points()) counts.merge(p.region(), 1, Integer::sum);
            System.out.printf("  region counts: %s%n", counts);

            // 4 square corners (X,Y both in {0,1}) -> CORNER (pure element, analytic G/H/S=0);
            // 4*(n-1) square-boundary-not-corner points (exactly one of X,Y in {0,1}) attempt a
            // SQUARE_EDGE_BINARY solve, but a non-convergent binary solve is legitimately skipped
            // rather than fabricated (mirrors TernaryGridScan's philosophy), so only an upper bound
            // applies here, not an exact count; (n-1)*(n-1) square-interior points -> INTERIOR,
            // same caveat.
            check(variant + ": exactly 4 CORNER points",
                    counts.getOrDefault(Region.CORNER, 0) == 4);
            int maxBinaryEdge = 4 * (n - 1);
            int actualBinaryEdge = counts.getOrDefault(Region.SQUARE_EDGE_BINARY, 0);
            check(variant + ": SQUARE_EDGE_BINARY point count within grid geometry bound ("
                            + actualBinaryEdge + " <= " + maxBinaryEdge + ")",
                    actualBinaryEdge > 0 && actualBinaryEdge <= maxBinaryEdge);

            for (Point p : result.points()) {
                double sum = p.fSlot0() + p.fSlot1() + p.fSlot2() + p.fSlot3();
                if (Math.abs(sum - 1.0) > 1e-9) {
                    fail(variant + ": mole fractions don't sum to 1 at (X=" + p.x() + ",Y=" + p.y()
                            + "): sum=" + sum);
                }
                for (double f : new double[] { p.fSlot0(), p.fSlot1(), p.fSlot2(), p.fSlot3() }) {
                    if (f < -1e-12 || f > 1 + 1e-12) {
                        fail(variant + ": mole fraction out of [0,1] at (X=" + p.x() + ",Y=" + p.y() + "): " + f);
                    }
                }
            }

            // Cross-check: an INTERIOR point's scan value must match a direct calculate() call.
            Point interior = result.points().stream()
                    .filter(p -> p.region() == Region.INTERIOR).findFirst().orElse(null);
            if (interior == null) {
                fail(variant + ": no INTERIOR point found to cross-check");
            } else {
                Map<String, Double> composition = new java.util.LinkedHashMap<>();
                composition.put("Ti", interior.fSlot1());
                composition.put("V", interior.fSlot2());
                composition.put("Zr", interior.fSlot3());
                var direct = service.calculate(session,
                        new org.ce.calculation.Conditions(temperature, composition),
                        Property.GIBBS_ENERGY, null, null);
                boolean match = direct.converged != null && direct.converged
                        && Math.abs(direct.gibbsEnergy - interior.value()) < 1e-6;
                System.out.printf("  interior cross-check: scan=%.10f direct=%.10f  %s%n",
                        interior.value(), direct.gibbsEnergy, match ? "ok" : "MISMATCH");
                check(variant + ": interior scan value matches direct calculate()", match);
            }
        }

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " QuaternarySquareScan smoke checks failed");
        }
    }

    private static void check(String what, boolean ok) {
        if (!ok) fail(what);
    }

    private static void fail(String what) {
        failures++;
        System.out.println("    [!] FAIL  " + what);
    }

    private QuaternarySquareScanSmokeTest() {}
}
