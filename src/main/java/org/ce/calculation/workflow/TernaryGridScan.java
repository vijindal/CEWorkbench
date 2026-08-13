package org.ce.calculation.workflow;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.Conditions;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.SroCalculator;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace.SystemId;

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
 * <p><b>Four composition regions, each handled differently:</b></p>
 * <ul>
 *   <li><b>Interior</b> (all three mole fractions &gt; 0) — a genuine ternary
 *       CVM solve. If it fails to converge, the point is simply skipped; no
 *       interpolation or synthetic data is substituted.</li>
 *   <li><b>Edge</b> (exactly one mole fraction = 0) — never evaluated on the
 *       ternary Hamiltonian. Swept separately as a 1-D scan on a genuine
 *       <em>binary</em> CVM session, extracted from the ternary Hamiltonian's
 *       inherited pair CECs via {@link BinarySubsystemExtractor}
 *       (binary-cluster CECs are unchanged in a higher-order system per
 *       Eq. 30 of the paper — no transformation, only extraction). The
 *       ternary solver is known to be numerically fragile in the composition
 *       band adjacent to a binary edge (see CLAUDE.md); the true binary
 *       problem has a much smaller configuration space and solves far more
 *       reliably.</li>
 *   <li><b>Corner</b> (two mole fractions = 0, i.e. a pure element) — for
 *       G/H/S this is trivially 0 by definition (nothing to mix with a
 *       single component); no calculation is run. For pair SRO the Eq. 40
 *       reference {@code x_P*x_R} is 0 there, making alpha mathematically
 *       undefined, so no point is emitted at all.</li>
 * </ul>
 *
 * <p>This mirrors how the reference dataset behind the SRO feature (a
 * companion spreadsheet for Nb-Ti-V) was itself generated: its composition
 * sweep never samples an exact zero-composition point either.</p>
 */
public final class TernaryGridScan {

    private TernaryGridScan() {}

    /**
     * Minimum mole fraction, for either element of a pair-SRO request, below
     * which the point is excluded rather than computed. Cowley-Warren alpha
     * (Eq. 40) divides by {@code x_P*x_R}: as either mole fraction shrinks
     * toward zero — whether at the pair's own edge or just approaching one of
     * the *other* two edges (where the third element dominates and one of the
     * pair's elements is incidentally small too) — that reference shrinks
     * with it, amplifying ordinary numerical noise in the converged cluster
     * probability into large, non-physical swings in alpha. Points below this
     * threshold are mathematically defined but not representative; excluding
     * them keeps the plotted/reported SRO values meaningful.
     */
    private static final double SRO_MIN_MOLE_FRACTION = 0.03;

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

    /** Which composition region a point came from. */
    public enum Region { INTERIOR, EDGE, CORNER }

    /** One computed grid point. */
    public record Point(double fa, double fb, double fc, double value, Region region) {}

    /** Full scan result: element order (a, b, c = triangle corners), temperature, quantity, points. */
    public record Result(List<String> elements, double temperature, Quantity quantity, List<Point> points,
                          int skipped) {}

    /**
     * Runs the scan: corners (trivial or omitted), each binary edge (1-D
     * binary CVM sweep), and the ternary interior (2-D ternary CVM sweep,
     * failures skipped — no interpolation).
     *
     * @param n           grid subdivisions per triangle edge
     * @param eventSink   receives {@link ProgressEvent.ScanPoint} after each
     *                    grid point completes; may be null
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

        Map<String, Double> cache = new LinkedHashMap<>();
        Map<String, ModelSession> binarySessions = new LinkedHashMap<>();

        List<Point> points = new ArrayList<>(total);
        int skipped = 0;
        int index = 0;

        for (double[] f : grid) {
            index++;
            int zeroCount = (f[0] == 0.0 ? 1 : 0) + (f[1] == 0.0 ? 1 : 0) + (f[2] == 0.0 ? 1 : 0);

            // For pair SRO, exclude any point where either paired element's
            // mole fraction is below SRO_MIN_MOLE_FRACTION — not just at the
            // pair's own edge, but anywhere in the triangle (e.g. approaching
            // one of the OTHER two edges, where the third element dominates
            // and one of the pair's elements is incidentally near-zero too).
            // See SRO_MIN_MOLE_FRACTION's doc for why this region is noisy
            // rather than simply small.
            if (quantity instanceof PairSroQuantity sq) {
                int ia = elements.indexOf(sq.elementA()), ib = elements.indexOf(sq.elementB());
                if (ia >= 0 && ib >= 0
                        && (f[ia] < SRO_MIN_MOLE_FRACTION || f[ib] < SRO_MIN_MOLE_FRACTION)
                        && f[ia] > 0.0 && f[ib] > 0.0) {
                    // Exact-zero cases (corners, non-pair edges) are already
                    // handled below with their own exclusion reasons; this
                    // only trims the noisy near-zero interior/edge band.
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }
            }

            Double value;
            Region region;
            if (zeroCount == 2) {
                // Corner (pure element): G/H/S = 0 by definition; SRO is
                // mathematically undefined (reference mole fraction is 0).
                if (quantity instanceof PairSroQuantity) {
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }
                value = 0.0;
                region = Region.CORNER;
            } else if (zeroCount == 1) {
                int axis = zeroAxis(f);
                int o1 = (axis + 1) % 3, o2 = (axis + 2) % 3;
                String elA = elements.get(o1), elB = elements.get(o2);
                double xB = f[o2]; // f[o1] + f[o2] == 1 here

                // For pair SRO, only that pair's own edge is physically
                // meaningful — the other two edges don't involve both of the
                // requested pair's elements at all (one is the excluded
                // axis). Not attempted, not counted as skipped, exactly like
                // a corner.
                if (quantity instanceof PairSroQuantity sq && !isPairEdge(sq, elA, elB)) {
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }

                ModelSession binarySession;
                try {
                    binarySession = getOrBuildBinarySession(session, elA, elB, binarySessions);
                } catch (Exception e) {
                    skipped++;
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue; // no matching CVCF terms for this pair
                }
                value = binaryEdgePointValue(service, binarySession, elA, elB, xB, temperature, quantity);
                region = Region.EDGE;
            } else {
                value = query(service, session, elements, temperature, property, quantity, f[0], f[1], f[2], cache);
                region = Region.INTERIOR;
            }

            if (value != null) {
                points.add(new Point(f[0], f[1], f[2], value, region));
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

    /** Solves the binary system at composition xB (fraction of elB) and extracts the requested quantity. */
    private static Double binaryEdgePointValue(CalculationService service, ModelSession binarySession,
            String elA, String elB, double xB, double temperature, Quantity quantity) throws Exception {

        Property property = (quantity instanceof PropertyQuantity pq) ? pq.property() : Property.GIBBS_ENERGY;
        Map<String, Double> composition = new LinkedHashMap<>();
        composition.put(elB, xB);
        Conditions conditions = new Conditions(temperature, composition);

        try {
            ThermodynamicResult r = service.calculate(binarySession, conditions, property, null, null);
            if (Boolean.FALSE.equals(r.converged)) return null;
            double value = extractQuantity(r, quantity, List.of(elA, elB));
            return Double.isNaN(value) ? null : value;
        } catch (Exception e) {
            return null;
        }
    }

    private static ModelSession getOrBuildBinarySession(ModelSession ternarySession, String elA, String elB,
            Map<String, ModelSession> binarySessions) throws Exception {
        String pairKey = elA + "-" + elB;
        ModelSession binarySession = binarySessions.get(pairKey);
        if (binarySession == null) {
            CECEntry binaryEntry = BinarySubsystemExtractor.extractBinary(
                    ternarySession.cecEntry, ternarySession.elements(), elA, elB);
            SystemId binarySystemId = new SystemId(binaryEntry.elements,
                    ternarySession.systemId.structure(), ternarySession.systemId.model());
            binarySession = new ModelSession.Builder(null)
                    .build(binarySystemId, ternarySession.engineConfig, binaryEntry, null);
            binarySessions.put(pairKey, binarySession);
        }
        return binarySession;
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

    /** Index (0,1,2) of the single zero component, for a point with exactly one zero. */
    private static int zeroAxis(double[] f) {
        if (f[0] == 0.0) return 0;
        if (f[1] == 0.0) return 1;
        return 2;
    }

    /** True if the edge's two elements are exactly {@code sq}'s pair (in either order). */
    private static boolean isPairEdge(PairSroQuantity sq, String elA, String elB) {
        return (elA.equalsIgnoreCase(sq.elementA()) && elB.equalsIgnoreCase(sq.elementB()))
                || (elA.equalsIgnoreCase(sq.elementB()) && elB.equalsIgnoreCase(sq.elementA()));
    }
}
