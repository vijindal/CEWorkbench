package org.ce.calculation.workflow;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.Conditions;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * (X,Y)-square composition sweep over a quaternary system at fixed
 * temperature, used to render the "square plots" of Fig. 20, Jindal &amp;
 * Lele 2025, CALPHAD 89, 102825.
 *
 * <p>The paper maps the unit square to the quaternary composition simplex
 * (a 3-D object) via four role slots, illustrated here with the paper's own
 * Nb-Ti-V-Zr example (slot0=Nb, slot1=Ti, slot2=V, slot3=Zr):</p>
 * <pre>
 *   x(slot0) = Y(1-X)
 *   x(slot1) = (1-X)(1-Y)
 *   x(slot2) = X(1-Y)
 *   x(slot3) = XY
 * </pre>
 * <p>at fixed X,Y &isin; [0,1]. Which system element plays which slot role is
 * caller policy (see {@link #run}'s {@code slotOrder} parameter) so this
 * class is not hardwired to any one quaternary system. This covers only a
 * 2-D slice of the full 3-D quaternary simplex, not all of it: a single
 * square reaches all four corners (X,Y both in {0,1}, a pure element) and
 * all four boundaries (exactly one of X,Y in {0,1}, a true binary edge of
 * the composition tetrahedron — see {@link Region}), but the interior of
 * the square is one particular family of paths through the tetrahedron's
 * interior, not every possible quaternary composition. The paper's Fig. 20
 * uses two such squares (a {@link Variant} each) specifically so that,
 * between the two, every one of the tetrahedron's six binary edges is
 * reached at least once.</p>
 *
 * <p>Mirrors {@link TernaryGridScan}'s general approach (composition-region
 * classification, true lower-order sessions at a boundary rather than
 * evaluating the parent Hamiltonian at zero composition, skipping rather
 * than patching non-convergent points) but the region taxonomy itself
 * differs — see {@link Region}'s doc for why the square has no ternary-edge
 * analogue of the triangle's {@code EDGE} region.</p>
 */
public final class QuaternarySquareScan {

    private QuaternarySquareScan() {}

    /**
     * Which square parametrization to use. The paper uses two, arranged so
     * that between them all six binary edges of the quaternary tetrahedron
     * are covered. Using the paper's own Nb-Ti-V-Zr slot assignment as a
     * concrete example: {@code STANDARD}'s four square edges are Ti-Zr
     * (X=0), Nb-V (X=1), Ti-V (Y=0), Nb-Zr (Y=1); the interior sweeps all
     * four elements together. {@code V_ZR_SWAPPED} swaps slots 2 and 3's
     * formula terms (V and Zr's roles, in that example), which reaches the
     * two edges {@code STANDARD} cannot reach at all: Nb-Ti and V-Zr become
     * the X=0/X=1 edges instead. For a different {@code slotOrder}, "V" and
     * "Zr" above simply mean "whatever occupies slot 2 / slot 3".
     */
    public enum Variant { STANDARD, V_ZR_SWAPPED }

    /** What to plot at each grid point — reuses the same quantity model as {@link TernaryGridScan}. */
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

    /**
     * Which composition region a point came from. Unlike the ternary
     * triangle, the (X,Y) square's parametrization never produces a point
     * with exactly one mole fraction at zero — working through the formula
     * at a boundary (say X=0): slot2=X(1-Y)=0 and slot3=XY=0 simultaneously,
     * so a square boundary always zeroes exactly <em>two</em> slots at once,
     * and a square corner (X,Y both in {0,1}) zeroes exactly
     * <em>three</em>. There is no {@code EDGE}/ternary-session region here
     * (contrast {@link TernaryGridScan.Region}, whose triangle boundary
     * zeroes exactly one component); only two regions are physically
     * reachable:
     * <ul>
     *   <li>{@code INTERIOR} — X,Y both strictly inside (0,1); all four mole
     *       fractions &gt; 0; a genuine quaternary CVM solve.</li>
     *   <li>{@code SQUARE_EDGE_BINARY} — exactly one of X,Y is 0 or 1 (a
     *       square boundary, not corner); exactly two mole fractions are 0;
     *       solved as a true binary CVM session (two-of-four elements).</li>
     *   <li>{@code CORNER} — both X and Y are 0 or 1 (a square corner);
     *       exactly three mole fractions are 0, i.e. a pure element. G/H/S
     *       are trivially 0 by definition; no calculation is run.</li>
     * </ul>
     */
    public enum Region { INTERIOR, SQUARE_EDGE_BINARY, CORNER }

    /**
     * One computed grid point, in (X,Y) square coordinates plus the resolved
     * mole fractions in slot order — {@code fSlot0} is the mole fraction of
     * {@code Result.elements().get(0)}, etc. (see {@code run}'s
     * {@code slotOrder} parameter).
     */
    public record Point(double x, double y, double fSlot0, double fSlot1, double fSlot2, double fSlot3,
                         double value, Region region) {}

    /** Full scan result. {@code elements} is the slot order passed to {@code run} (slot0..slot3). */
    public record Result(List<String> elements, Variant variant, double temperature, Quantity quantity,
                          List<Point> points, int skipped) {}

    /**
     * Runs the scan over an n&times;n grid of the unit (X,Y) square.
     *
     * @param session     a quaternary CVM session (4 elements)
     * @param slotOrder   the session's 4 elements, reordered so that
     *                    {@code slotOrder.get(0)} plays "Nb"'s role in the
     *                    paper's formula, {@code slotOrder.get(1)} plays
     *                    "Ti"'s role, {@code slotOrder.get(2)} plays "V"'s
     *                    role, and {@code slotOrder.get(3)} plays "Zr"'s
     *                    role (for {@code V_ZR_SWAPPED}, slots 2 and 3 swap
     *                    which formula term feeds them, but the slot-to-
     *                    element assignment itself is unchanged — pass the
     *                    same {@code slotOrder} for both variants). This is
     *                    independent of {@code session}'s own canonical
     *                    element order; slot assignment is caller policy,
     *                    typically "reproduce the paper's Fig. 20 for system
     *                    X" by naming which of X's four elements should be
     *                    treated as Nb/Ti/V/Zr.
     * @param n           grid subdivisions per axis (n+1 points per axis)
     * @param eventSink   receives {@link ProgressEvent.ScanPoint} after each grid point; may be null
     */
    public static Result run(CalculationService service, ModelSession session, List<String> slotOrder,
            Variant variant, double temperature, Quantity quantity, int n, Consumer<ProgressEvent> eventSink)
            throws Exception {

        if (slotOrder.size() != 4)
            throw new IllegalArgumentException("Quaternary square scan requires exactly 4 elements, got " + slotOrder);
        if (session.elements().size() != 4)
            throw new IllegalArgumentException("Quaternary square scan requires a 4-element session, got "
                    + session.elements());
        for (String el : slotOrder) {
            if (session.elements().stream().noneMatch(el::equalsIgnoreCase)) {
                throw new IllegalArgumentException("slotOrder element '" + el
                        + "' is not one of the session's elements " + session.elements());
            }
        }
        Property property = (quantity instanceof PropertyQuantity pq) ? pq.property() : Property.GIBBS_ENERGY;

        List<double[]> grid = squareGrid(n);
        int total = grid.size();

        Map<String, Double> cache = new LinkedHashMap<>();
        Map<String, ModelSession> binarySessions = new LinkedHashMap<>();

        List<Point> points = new ArrayList<>(total);
        int skipped = 0;
        int index = 0;

        for (double[] xy : grid) {
            index++;
            double X = xy[0], Y = xy[1];
            // frac[i] is the mole fraction of slotOrder.get(i) — slot order fixed by
            // the parameter, not by toMoleFractions (V_ZR_SWAPPED only changes which
            // formula expression feeds slots 2/3, never which element owns which slot).
            double[] frac = toMoleFractions(X, Y, variant);
            List<String> order = slotOrder;

            int zeroCount = 0;
            for (double f : frac) if (f == 0.0) zeroCount++;

            Double value;
            Region region;

            if (zeroCount == 3) {
                // Square corner (X,Y both in {0,1}): a pure element. G/H/S trivially 0;
                // SRO's reference mole fraction is 0, making alpha undefined (mirrors
                // TernaryGridScan's corner handling).
                if (quantity instanceof PairSroQuantity) {
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }
                value = 0.0;
                region = Region.CORNER;
            } else if (zeroCount == 2) {
                // Square boundary (exactly one of X,Y in {0,1}): a true binary edge.
                String[] pair = nonZeroPair(order, frac);
                if (pair == null) {
                    skipped++;
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }
                String elA = pair[0], elB = pair[1];
                double xB = fractionOf(order, frac, elB); // fractions of elA/elB sum to 1 here

                if (quantity instanceof PairSroQuantity sq && !isPair(sq, elA, elB)) {
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }

                ModelSession binarySession;
                try {
                    binarySession = getOrBuildBinarySession(session, elA, elB, binarySessions);
                } catch (Exception e) {
                    skipped++;
                    if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
                    continue;
                }
                value = twoElementPointValue(service, binarySession, List.of(elA, elB),
                        Map.of(elB, xB), temperature, quantity);
                region = Region.SQUARE_EDGE_BINARY;
            } else {
                Map<String, Double> composition = new LinkedHashMap<>();
                composition.put(order.get(1), frac[1]);
                composition.put(order.get(2), frac[2]);
                composition.put(order.get(3), frac[3]);
                value = query(service, session, order, temperature, property, quantity, frac, cache);
                region = Region.INTERIOR;
            }

            if (value != null) {
                // frac slots are fixed as [slotOrder(0)="Nb", 1="Ti", 2="V", 3="Zr"]
                // by toMoleFractions regardless of variant; see run()'s slotOrder doc.
                points.add(new Point(X, Y, frac[0], frac[1], frac[2], frac[3], value, region));
            } else {
                skipped++;
            }
            if (eventSink != null) eventSink.accept(new ProgressEvent.ScanPoint(index, total));
        }

        return new Result(slotOrder, variant, temperature, quantity, points, skipped);
    }

    /**
     * Maps (X,Y) to mole fractions in the caller's [Nb,Ti,V,Zr] slot order
     * (slot 0 = "Nb"'s role, ..., slot 3 = "Zr"'s role) per the paper's
     * formula. {@code V_ZR_SWAPPED} swaps which formula expression feeds
     * slots 2 and 3 (not which element owns which slot — that's fixed by the
     * caller's {@code slotOrder}) so the square's edges trace a different
     * pair of the tetrahedron's six binary edges — see {@link Variant}'s doc.
     */
    private static double[] toMoleFractions(double X, double Y, Variant variant) {
        double slot0 = Y * (1 - X);
        double slot1 = (1 - X) * (1 - Y);
        double slot2 = X * (1 - Y);
        double slot3 = X * Y;
        if (variant == Variant.V_ZR_SWAPPED) {
            double tmp = slot2;
            slot2 = slot3;
            slot3 = tmp;
        }
        return new double[] { slot0, slot1, slot2, slot3 };
    }

    private static double fractionOf(List<String> order, double[] frac, String element) {
        int i = order.indexOf(element);
        if (i < 0) throw new IllegalArgumentException("Element '" + element + "' not in order " + order);
        return frac[i];
    }

    private static String[] nonZeroPair(List<String> order, double[] frac) {
        List<String> nonZero = new ArrayList<>();
        for (int i = 0; i < frac.length; i++) if (frac[i] > 0.0) nonZero.add(order.get(i));
        if (nonZero.size() != 2) return null;
        return new String[] { nonZero.get(0), nonZero.get(1) };
    }

    private static boolean isPair(PairSroQuantity sq, String elA, String elB) {
        return (elA.equalsIgnoreCase(sq.elementA()) && elB.equalsIgnoreCase(sq.elementB()))
                || (elA.equalsIgnoreCase(sq.elementB()) && elB.equalsIgnoreCase(sq.elementA()));
    }

    private static Double query(CalculationService service, ModelSession session, List<String> order,
            double temperature, Property property, Quantity quantity, double[] frac,
            Map<String, Double> cache) throws Exception {
        String key = frac[0] + "," + frac[1] + "," + frac[2] + "," + frac[3];
        if (cache.containsKey(key)) return cache.get(key);

        Map<String, Double> composition = new LinkedHashMap<>();
        composition.put(order.get(1), frac[1]);
        composition.put(order.get(2), frac[2]);
        composition.put(order.get(3), frac[3]);
        Conditions conditions = new Conditions(temperature, composition);

        Double result;
        try {
            ThermodynamicResult r = service.calculate(session, conditions, property, null, null);
            result = (Boolean.FALSE.equals(r.converged)) ? null : extractQuantity(r, quantity, order);
        } catch (Exception e) {
            result = null;
        }
        cache.put(key, result);
        return result;
    }

    /** Solves a binary session at the given composition and extracts the requested quantity. */
    private static Double twoElementPointValue(CalculationService service, ModelSession lowerSession,
            List<String> lowerElements, Map<String, Double> composition, double temperature, Quantity quantity)
            throws Exception {

        Property property = (quantity instanceof PropertyQuantity pq) ? pq.property() : Property.GIBBS_ENERGY;
        Conditions conditions = new Conditions(temperature, composition);

        try {
            ThermodynamicResult r = service.calculate(lowerSession, conditions, property, null, null);
            if (Boolean.FALSE.equals(r.converged)) return null;
            double value = extractQuantity(r, quantity, lowerElements);
            return Double.isNaN(value) ? null : value;
        } catch (Exception e) {
            return null;
        }
    }

    private static ModelSession getOrBuildBinarySession(ModelSession parentSession, String elA, String elB,
            Map<String, ModelSession> binarySessions) throws Exception {
        String pairKey = elA + "-" + elB;
        ModelSession binarySession = binarySessions.get(pairKey);
        if (binarySession == null) {
            CECEntry binaryEntry = BinarySubsystemExtractor.extractBinary(
                    parentSession.cecEntry, parentSession.elements(), elA, elB);
            SystemId binarySystemId = new SystemId(binaryEntry.elements,
                    parentSession.systemId.structure(), parentSession.systemId.model());
            binarySession = new ModelSession.Builder(null)
                    .build(binarySystemId, parentSession.engineConfig, binaryEntry, null);
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
            List<org.ce.model.cvm.CVMGibbsModel.PairSro> shell = r.sro.get("1NN");
            if (shell == null) return Double.NaN;
            int ia = elements.indexOf(sq.elementA());
            int ib = elements.indexOf(sq.elementB());
            if (ia < 0 || ib < 0) return Double.NaN;
            for (org.ce.model.cvm.CVMGibbsModel.PairSro pair : shell) {
                if ((pair.i == ia && pair.j == ib) || (pair.i == ib && pair.j == ia)) {
                    return pair.alpha;
                }
            }
            return Double.NaN;
        }
        return Double.NaN;
    }

    /** Full n&times;n grid over the unit square [0,1]&times;[0,1], n+1 points per axis. */
    static List<double[]> squareGrid(int n) {
        List<double[]> points = new ArrayList<>();
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                points.add(new double[] { (double) i / n, (double) j / n });
            }
        }
        return points;
    }
}
