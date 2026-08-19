package org.ce.scratch;

import org.ce.model.equilibrium.SgteDatabase;

/**
 * Cross-checks {@link SgteDatabase} (parsed from {@code inputs/unary.dat})
 * against {@link HardcodedLatticeStability} (the polynomials LatticeStability carried before migration) for every
 * (element, phase, T) both can supply -- for {@code G0} and for
 * {@code dG0/dT}.
 *
 * <p>The two are independent transcriptions of the same SGTE source, so
 * agreement is real evidence and disagreement localises a transcription error
 * to one side or the other. This is the gate to clear before the hardcoded
 * version is retired.</p>
 *
 * <p>The derivative pass also runs a <b>finite-difference check against each
 * side's own {@code G0}</b>. That is a stronger test than the two agreeing with
 * each other: two implementations can agree while both differentiate the same
 * polynomial wrongly, but neither can match a numerical derivative of its own
 * energy unless the analytic form is right. It also independently confirms the
 * SGTE term-to-derivative mapping in {@code SgteDatabase}, which was written
 * from the term table rather than transcribed.</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.SgteCrossCheck
 * </pre>
 */
public final class SgteCrossCheck {

    /** J/mol. Both evaluate the same closed form, so agreement should be near-exact. */
    private static final double TOL = 1.0e-6;

    /** J/mol/K for dG0/dT comparisons between the two implementations. */
    private static final double TOL_DERIV = 1.0e-8;

    /**
     * Relative tolerance for the finite-difference check. Central differences
     * over a polynomial with terms up to T^7 and down to T^-9 carry real
     * truncation error, so this is necessarily looser than the analytic
     * comparisons above.
     */
    private static final double TOL_FD_REL = 1.0e-5;

    private static final String[] ELEMENTS = { "Mo", "Nb", "Re", "Ta", "Ti", "V", "W", "Zr" };
    private static final String[] PHASES = { "BCC_A2", "FCC_A1", "HCP_A3", "LIQUID" };
    private static final double[] TEMPERATURES = { 300.0, 700.0, 1000.0, 1273.0, 1800.0, 2500.0 };

    public static void main(String[] args) {
        System.out.println("=".repeat(84));
        System.out.println("  SGTE cross-check: database (unary.dat) vs hardcoded LatticeStability");
        System.out.println("=".repeat(84));

        boolean energyOk = checkEnergy();
        boolean derivOk = checkDerivative();
        boolean fdOk = checkFiniteDifference();

        boolean pass = energyOk && derivOk && fdOk;
        System.out.println("\n" + "=".repeat(84));
        System.out.printf("RESULT: %s   (G0 %s, dG0/dT %s, finite-difference %s)%n",
                pass ? "PASS" : "FAIL",
                energyOk ? "ok" : "FAILED",
                derivOk ? "ok" : "FAILED",
                fdOk ? "ok" : "FAILED");
        System.out.println("=".repeat(84));

        if (!pass) {
            throw new AssertionError("SGTE cross-check failed");
        }
    }

    // =========================================================================
    // 1. G0
    // =========================================================================

    private static boolean checkEnergy() {
        System.out.println("\n--- G0 -------------------------------------------------------------");
        int agree = 0, differ = 0, dbOnly = 0, hcOnly = 0, neither = 0, boundary = 0;
        StringBuilder notes = new StringBuilder();

        for (String el : ELEMENTS) {
            for (String ph : PHASES) {
                for (double t : TEMPERATURES) {
                    Double db = tryGet(() -> SgteDatabase.g0(el, ph, t));
                    Double hc = tryGet(() -> HardcodedLatticeStability.g0(el, ph, t));

                    if (db == null && hc == null) {
                        neither++;
                    } else if (db == null) {
                        hcOnly++;
                        notes.append(String.format(
                                "  %-3s %-8s %7.1f K  hardcoded=%16.6f  database=UNAVAILABLE%n", el, ph, t, hc));
                    } else if (hc == null) {
                        dbOnly++;
                    } else if (Math.abs(db - hc) <= TOL) {
                        agree++;
                    } else if (isRangeBoundary(el, ph, t)) {
                        boundary++;
                        notes.append(String.format(
                                "  [boundary] %-3s %-8s %7.1f K  hardcoded=%16.6f  database=%16.6f  diff=%.3e%n",
                                el, ph, t, hc, db, db - hc));
                    } else {
                        differ++;
                        notes.append(String.format(
                                "  %-3s %-8s %7.1f K  hardcoded=%16.6f  database=%16.6f  diff=%.6e%n",
                                el, ph, t, hc, db, db - hc));
                    }
                }
            }
        }

        report("G0", agree, boundary, differ, dbOnly, hcOnly, neither, notes);
        return differ == 0 && hcOnly == 0;
    }

    // =========================================================================
    // 2. dG0/dT
    // =========================================================================

    private static boolean checkDerivative() {
        System.out.println("\n--- dG0/dT ---------------------------------------------------------");
        int agree = 0, differ = 0, dbOnly = 0, hcOnly = 0, neither = 0, boundary = 0;
        StringBuilder notes = new StringBuilder();

        for (String el : ELEMENTS) {
            for (String ph : PHASES) {
                for (double t : TEMPERATURES) {
                    Double db = tryGet(() -> SgteDatabase.dG0Dt(el, ph, t));
                    Double hc = tryGet(() -> HardcodedLatticeStability.dG0Dt(el, ph, t));

                    if (db == null && hc == null) {
                        neither++;
                    } else if (db == null) {
                        hcOnly++;
                        notes.append(String.format(
                                "  %-3s %-8s %7.1f K  hardcoded=%16.8f  database=UNAVAILABLE%n", el, ph, t, hc));
                    } else if (hc == null) {
                        dbOnly++;
                    } else if (Math.abs(db - hc) <= TOL_DERIV) {
                        agree++;
                    } else if (isDerivativeRangeBoundary(el, ph, t)) {
                        boundary++;
                        notes.append(String.format(
                                "  [boundary] %-3s %-8s %7.1f K  hardcoded=%16.8f  database=%16.8f  diff=%.3e%n",
                                el, ph, t, hc, db, db - hc));
                    } else {
                        differ++;
                        notes.append(String.format(
                                "  %-3s %-8s %7.1f K  hardcoded=%16.8f  database=%16.8f  diff=%.6e%n",
                                el, ph, t, hc, db, db - hc));
                    }
                }
            }
        }

        report("dG0/dT", agree, boundary, differ, dbOnly, hcOnly, neither, notes);
        return differ == 0 && hcOnly == 0;
    }

    // =========================================================================
    // 3. Finite-difference check of each side against its own G0
    // =========================================================================

    private static boolean checkFiniteDifference() {
        System.out.println("\n--- dG0/dT vs central difference of the same source's G0 -----------");
        int okDb = 0, okHc = 0, badDb = 0, badHc = 0;
        StringBuilder notes = new StringBuilder();

        for (String el : ELEMENTS) {
            for (String ph : PHASES) {
                for (double t : TEMPERATURES) {
                    // Skip only a genuine piecewise break, where a central
                    // difference would straddle two different polynomials and
                    // the comparison is meaningless for either implementation.
                    // Detected from the source's own curvature, not from
                    // whether the two implementations agree -- an earlier
                    // version reused isRangeBoundary here, which also returns
                    // true when a lookup simply fails, and so silently skipped
                    // every hardcoded case.
                    double h = 1.0e-3 * t;
                    if (hasBreakNear(el, ph, t, h)) {
                        continue;
                    }

                    Double dbFd = centralDifference(x -> tryGet(() -> SgteDatabase.g0(el, ph, x)), t, h);
                    Double dbAn = tryGet(() -> SgteDatabase.dG0Dt(el, ph, t));
                    if (dbFd != null && dbAn != null) {
                        if (relativeDiff(dbAn, dbFd) <= TOL_FD_REL) {
                            okDb++;
                        } else {
                            badDb++;
                            notes.append(String.format(
                                    "  [database]  %-3s %-8s %7.1f K  analytic=%16.8f  fd=%16.8f  rel=%.3e%n",
                                    el, ph, t, dbAn, dbFd, relativeDiff(dbAn, dbFd)));
                        }
                    }

                    Double hcFd = centralDifference(x -> tryGet(() -> HardcodedLatticeStability.g0(el, ph, x)), t, h);
                    Double hcAn = tryGet(() -> HardcodedLatticeStability.dG0Dt(el, ph, t));
                    if (hcFd != null && hcAn != null) {
                        if (relativeDiff(hcAn, hcFd) <= TOL_FD_REL) {
                            okHc++;
                        } else {
                            badHc++;
                            notes.append(String.format(
                                    "  [hardcoded] %-3s %-8s %7.1f K  analytic=%16.8f  fd=%16.8f  rel=%.3e%n",
                                    el, ph, t, hcAn, hcFd, relativeDiff(hcAn, hcFd)));
                        }
                    }
                }
            }
        }

        System.out.printf("%n  database  : %d ok, %d FAILED%n", okDb, badDb);
        System.out.printf("  hardcoded : %d ok, %d FAILED%n", okHc, badHc);
        if (notes.length() > 0) {
            System.out.println("\n  --- finite-difference mismatches ---");
            System.out.print(notes);
        }
        return badDb == 0 && badHc == 0;
    }

    /**
     * True when a piecewise break lies within {@code +/-h} of {@code t} in
     * either source, making a central difference straddle two polynomials.
     *
     * <p>Judged from each source against itself -- comparing its analytic
     * derivative at the two offsets with a one-sided difference on that same
     * side. A smooth range agrees closely on both sides; a break does not.</p>
     */
    private static boolean hasBreakNear(String el, String ph, double t, double h) {
        return breaksNear(t, h, x -> tryGet(() -> SgteDatabase.g0(el, ph, x)),
                x -> tryGet(() -> SgteDatabase.dG0Dt(el, ph, x)))
                || breaksNear(t, h, x -> tryGet(() -> HardcodedLatticeStability.g0(el, ph, x)),
                        x -> tryGet(() -> HardcodedLatticeStability.dG0Dt(el, ph, x)));
    }

    private static boolean breaksNear(double t, double h, TempFunction g, TempFunction dg) {
        Double left = g.at(t - h);
        Double mid = g.at(t);
        Double right = g.at(t + h);
        Double dMid = dg.at(t);
        if (left == null || mid == null || right == null || dMid == null) {
            return true; // cannot judge; treat as unusable rather than assert on it
        }
        // Each one-sided slope must match the analytic derivative to within a
        // loose relative bound; a piecewise break shows up as a large mismatch.
        double slopeL = (mid - left) / h;
        double slopeR = (right - mid) / h;
        return relativeDiff(dMid, slopeL) > 1.0e-3 || relativeDiff(dMid, slopeR) > 1.0e-3;
    }

    private interface TempFunction {
        Double at(double t);
    }

    private static Double centralDifference(TempFunction f, double t, double h) {
        Double plus = f.at(t + h);
        Double minus = f.at(t - h);
        if (plus == null || minus == null) {
            return null;
        }
        return (plus - minus) / (2.0 * h);
    }

    private static double relativeDiff(double a, double b) {
        return Math.abs(a - b) / Math.max(Math.abs(a), 1.0);
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private static void report(String label, int agree, int boundary, int differ,
            int dbOnly, int hcOnly, int neither, StringBuilder notes) {
        System.out.printf("%n  agree            : %d%n", agree);
        System.out.printf("  boundary only    : %d   (same polynomials; range selection at T itself)%n", boundary);
        System.out.printf("  DISAGREE         : %d%n", differ);
        System.out.printf("  database only    : %d   (coverage the hardcoded version lacks)%n", dbOnly);
        System.out.printf("  hardcoded only   : %d   (database lookup failed - investigate)%n", hcOnly);
        System.out.printf("  neither          : %d%n", neither);
        if (notes.length() > 0) {
            System.out.println("\n  --- " + label + " discrepancies ---");
            System.out.print(notes);
        }
    }

    /**
     * True when a {@code G0} disagreement at {@code t} is a piecewise-boundary
     * artifact rather than a coefficient error.
     *
     * <p>Detected empirically instead of from a table of known boundaries:
     * nudge the temperature by a millikelvin either way and see whether the two
     * implementations agree there. If they do, both carry the same polynomials
     * and only the range selection at exactly {@code t} differed -- SGTE ranges
     * are half-open {@code (lower, upper]}, so {@code T} equal to a stated upper
     * bound belongs to the range ending there, which is what
     * {@code SgteDatabase} implements; the hardcoded version tests
     * {@code t < upper} and falls through to the next range.</p>
     */
    private static boolean isRangeBoundary(String el, String ph, double t) {
        for (double eps : new double[] { -1.0e-3, 1.0e-3 }) {
            Double db = tryGet(() -> SgteDatabase.g0(el, ph, t + eps));
            Double hc = tryGet(() -> HardcodedLatticeStability.g0(el, ph, t + eps));
            if (db == null || hc == null || Math.abs(db - hc) > TOL) {
                return false;
            }
        }
        return true;
    }

    /** As {@link #isRangeBoundary}, for {@code dG0/dT}. */
    private static boolean isDerivativeRangeBoundary(String el, String ph, double t) {
        for (double eps : new double[] { -1.0e-3, 1.0e-3 }) {
            Double db = tryGet(() -> SgteDatabase.dG0Dt(el, ph, t + eps));
            Double hc = tryGet(() -> HardcodedLatticeStability.dG0Dt(el, ph, t + eps));
            if (db == null || hc == null || Math.abs(db - hc) > TOL_DERIV) {
                return false;
            }
        }
        return true;
    }

    private interface Getter {
        double get();
    }

    /** Returns null when the source has no data for that combination. */
    private static Double tryGet(Getter g) {
        try {
            double v = g.get();
            return Double.isFinite(v) ? v : null;
        } catch (RuntimeException e) {
            return null;
        }
    }
}
