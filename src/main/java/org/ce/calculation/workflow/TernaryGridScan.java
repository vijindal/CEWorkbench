package org.ce.calculation.workflow;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.Conditions;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.SroCalculator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Composition-grid sweep over a ternary system at fixed temperature, used to
 * render isothermal sections (e.g. Figs. 15-20 of Jindal &amp; Lele 2025,
 * CALPHAD 89, 102825).
 *
 * <p>{@link org.ce.calculation.ConditionsScan} deliberately supports only one
 * varying axis at a time (see its class doc), so a true 2-D composition sweep
 * is driven here directly against {@link CalculationService#calculate}
 * point-by-point, reusing the session cache rather than routing through
 * {@code calculateScan}.</p>
 *
 * <p><b>Near-edge non-convergence:</b> the ternary CVM Newton-Raphson solver
 * can fail to converge in a thin composition band adjacent to a binary edge,
 * even though it converges exactly on the edge (one component = 0) and
 * further into the interior. This is a solver limitation, not a scan bug —
 * see CLAUDE.md discussion. Rather than leave a gap, such points are bridged
 * by linear interpolation (in barycentric distance from the edge) between the
 * exact edge value and a converged interior point on the same composition
 * ray.</p>
 */
public final class TernaryGridScan {

    private TernaryGridScan() {}

    /**
     * What to plot at each grid point. Either a standard thermodynamic
     * {@link Property} (G/H/S) or the 1st-neighbour Cowley-Warren pair SRO
     * parameter for one unlike species pair.
     *
     * <p>Every {@link Property} calculation already computes both SRO shells
     * as a side effect (see {@code ThermodynamicWorkflow.computeSro}), so
     * {@code PairSroQuantity} reuses the same single-point {@code calculate}
     * call and just extracts a different field from the result — no separate
     * SRO-only calculation path is needed.</p>
     *
     * <p>Only 1st-nearest-neighbour ("1NN") pair SRO is exposed for now.
     * Triangle/tetrahedron multi-site SRO and the 2NN shell are deferred —
     * see the CLAUDE.md discussion on why some CVCF correlation functions
     * (e.g. the ternary binary-triangle CFs) are antisymmetric differences of
     * two probabilities rather than a single probability, and so do not have
     * a directly meaningful Cowley-Warren-style alpha.</p>
     */
    public sealed interface Quantity {
        String label();
    }

    /** A standard G/H/S calculation. */
    public record PropertyQuantity(Property property) implements Quantity {
        @Override public String label() { return property.name(); }
    }

    /** 1st-neighbour Cowley-Warren SRO alpha for one unlike species pair, e.g. "Nb-Ti". */
    public record PairSroQuantity(String elementA, String elementB) implements Quantity {
        @Override public String label() { return "SRO_" + elementA + "-" + elementB + "_1NN"; }
    }

    /** One computed (or interpolated) grid point. */
    public record Point(double fa, double fb, double fc, double value, boolean interpolated) {}

    /** Full scan result: element order (a, b, c = triangle corners), temperature, quantity, points. */
    public record Result(List<String> elements, double temperature, Quantity quantity, List<Point> points,
                          int skipped) {}

    /**
     * Runs the scan. {@code elements} must be exactly the 3 elements of the
     * session's system, in the order corresponding to the desired triangle
     * corners (a, b, c).
     *
     * @param n           grid subdivisions per triangle edge (point count is
     *                    (n+1)(n+2)/2)
     * @param eventSink   receives {@link ProgressEvent.ScanPoint} after each
     *                    grid point completes (converged, interpolated, or
     *                    skipped); may be null
     */
    public static Result run(CalculationService service, ModelSession session, List<String> elements,
            double temperature, Quantity quantity, int n, Consumer<ProgressEvent> eventSink) throws Exception {

        if (elements.size() != 3)
            throw new IllegalArgumentException("Ternary grid scan requires exactly 3 elements, got " + elements);
        // CVM always computes GIBBS_ENERGY internally regardless of the requested
        // property (the equilibrium state is the same); requesting it explicitly
        // here keeps SRO extraction on the same well-tested code path as a normal
        // property scan.
        Property property = (quantity instanceof PropertyQuantity pq) ? pq.property() : Property.GIBBS_ENERGY;

        List<double[]> grid = barycentricGrid(n);
        int total = grid.size();

        // Cache single-point results within this scan so edge/interior lookups used
        // for interpolation don't repeat identical calculations.
        Map<String, Double> cache = new LinkedHashMap<>();

        List<Point> points = new ArrayList<>(total);
        int skipped = 0;
        int index = 0;
        for (double[] f : grid) {
            index++;
            Double value = query(service, session, elements, temperature, property, quantity, f[0], f[1], f[2], cache);
            boolean interpolated = false;
            if (value == null) {
                double[][] edgeAndInterior = nearestEdgeAndInterior(f[0], f[1], f[2]);
                if (edgeAndInterior != null) {
                    double[] edge = edgeAndInterior[0];
                    double[] interior = edgeAndInterior[1];
                    double t = edgeAndInterior[2][0];
                    Double edgeVal = query(service, session, elements, temperature, property, quantity,
                            edge[0], edge[1], edge[2], cache);
                    Double interiorVal = query(service, session, elements, temperature, property, quantity,
                            interior[0], interior[1], interior[2], cache);
                    if (edgeVal != null && interiorVal != null) {
                        value = (1 - t) * edgeVal + t * interiorVal;
                        interpolated = true;
                    }
                }
            }
            if (value != null) {
                points.add(new Point(f[0], f[1], f[2], value, interpolated));
            } else {
                skipped++;
            }
            if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
        }

        return new Result(elements, temperature, quantity, points, skipped);
    }

    private static Double query(CalculationService service, ModelSession session, List<String> elements,
            double temperature, Property property, Quantity quantity, double fa, double fb, double fc,
            Map<String, Double> cache) throws Exception {
        String key = fa + "," + fb + "," + fc;
        if (cache.containsKey(key)) return cache.get(key);

        Map<String, Double> composition = new LinkedHashMap<>();
        composition.put(elements.get(1), fb);
        composition.put(elements.get(2), fc);
        Conditions conditions = new Conditions(temperature, composition);

        Double result;
        try {
            ThermodynamicResult r = service.calculate(session, conditions, property, null, null);
            result = (Boolean.FALSE.equals(r.converged)) ? null : extractQuantity(r, quantity, elements);
        } catch (Exception e) {
            result = null;
        }
        cache.put(key, result);
        return result;
    }

    private static double extractQuantity(ThermodynamicResult r, Quantity quantity, List<String> elements) {
        if (quantity instanceof PropertyQuantity pq) {
            return switch (pq.property()) {
                case GIBBS_ENERGY -> r.gibbsEnergy;
                case ENTHALPY -> r.enthalpy;
                case ENTROPY -> r.entropy;
                default -> Double.NaN;
            };
        }
        if (quantity instanceof PairSroQuantity sq) {
            if (r.sro == null) return Double.NaN;
            List<SroCalculator.PairSro> shell = r.sro.get("1NN");
            if (shell == null) return Double.NaN;
            int ia = elements.indexOf(sq.elementA());
            int ib = elements.indexOf(sq.elementB());
            if (ia < 0 || ib < 0) return Double.NaN;
            for (SroCalculator.PairSro pair : shell) {
                if ((pair.i == ia && pair.j == ib) || (pair.i == ib && pair.j == ia)) {
                    return pair.alpha;
                }
            }
            return Double.NaN;
        }
        return Double.NaN;
    }

    /** Barycentric grid over the 2-simplex with n subdivisions per edge. Each row is {fa, fb, fc}. */
    static List<double[]> barycentricGrid(int n) {
        List<double[]> points = new ArrayList<>();
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n - i; j++) {
                int k = n - i - j;
                points.add(new double[] { (double) i / n, (double) j / n, (double) k / n });
            }
        }
        return points;
    }

    /**
     * For composition (fa, fb, fc), finds the nearest edge (component -> 0) and a
     * point further along the same ray into the interior.
     *
     * @return {edgePoint, interiorPoint, {t}} where t is the barycentric fraction
     *         of the original point between edgePoint (t=0) and interiorPoint
     *         (t=1); or null if no valid ray could be constructed.
     */
    private static double[][] nearestEdgeAndInterior(double fa, double fb, double fc) {
        double[] comps = { fa, fb, fc };
        int axis = 0;
        for (int i = 1; i < 3; i++) if (comps[i] < comps[axis]) axis = i;
        double tEdge = comps[axis];

        double[] edge = scale(comps, axis, 0.0);
        double interiorTarget = Math.min(tEdge + 0.08, 1.0);
        if (edge == null || interiorTarget <= tEdge) return null;
        double[] interior = scale(comps, axis, interiorTarget);
        if (interior == null) return null;

        double t = tEdge / interiorTarget;
        return new double[][] { edge, interior, { t } };
    }

    /**
     * Rescales composition so {@code comps[axis] == targetAxisValue}, keeping the
     * ratio of the other two components fixed (moves along the ray toward/away
     * from the edge where axis == 0).
     */
    private static double[] scale(double[] comps, int axis, double targetAxisValue) {
        int o1 = (axis + 1) % 3, o2 = (axis + 2) % 3;
        double remaining = 1.0 - targetAxisValue;
        double otherSum = comps[o1] + comps[o2];
        if (otherSum <= 0) return null;
        double[] scaled = comps.clone();
        scaled[axis] = targetAxisValue;
        scaled[o1] = comps[o1] / otherSum * remaining;
        scaled[o2] = comps[o2] / otherSum * remaining;
        return scaled;
    }
}
