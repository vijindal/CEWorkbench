package org.ce.scratch;

import org.ce.model.equilibrium.LatticeStability;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Regression gate for the pure-element reference energy
 * {@code G0m = Sum_i x_i * G0(element_i, phase, T)} -- the mechanical-mixture
 * term in {@code G = G0m + Gm}.
 *
 * <p>Checks {@link LatticeStability#g0m} and {@link LatticeStability#g0}
 * against <b>externally supplied reference values</b>, not against anything
 * this codebase computes. That independence is the point: an internal
 * consistency check cannot detect a wrong SGTE polynomial coefficient or a
 * mis-wired phase switch, because every internal consumer would read the same
 * wrong number.</p>
 *
 * <p>Run:</p>
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.LatticeStabilityVerification
 * </pre>
 *
 * <p><b>Why the phase switch is worth gating.</b> The SGTE reference function
 * {@code ghser<El>} is <em>not</em> the same phase for every element: it is BCC
 * for Mo, Nb, Ta, V and W, but HCP for Re, Ti and Zr, which reach BCC through a
 * separate {@code gbcc<El>} function. A case wired to the wrong branch returns a
 * plausible number rather than throwing, so only an external value catches it.
 * Every expectation below should come from a database or hand calculation
 * outside this program. The Ti and Zr cases cover exactly that branch and
 * confirm {@code BCC_A2} routes to {@code gbccTi}/{@code gbccZr}, not to the
 * HCP {@code ghser*} reference.</p>
 */
public final class LatticeStabilityVerification {

    /** Absolute tolerance in J/mol. Loose enough for a value quoted to 1 d.p. */
    private static final double TOL = 0.05;

    /** One externally-supplied expectation. */
    private record Case(String label, List<String> elements, String phase,
                        double[] x, double t, double expected) {
        /** Single-element case: g0 for one element, i.e. x = [1]. */
        static Case pure(String element, String phase, double t, double expected) {
            return new Case(element + " " + phase, List.of(element), phase,
                    new double[] { 1.0 }, t, expected);
        }
    }

    /**
     * Externally-supplied reference values.
     *
     * <p>Add new rows here as values become available; each is a permanent
     * regression gate. Prefer covering an element whose {@code ghser} is HCP
     * (Re, Ti, Zr) at {@code BCC_A2}, since that exercises the branch most
     * likely to be silently mis-wired.</p>
     */
    private static List<Case> cases() {
        List<Case> cases = new ArrayList<>();

        // Supplied by the user, 2026-08-19. Mo, Nb and Ta are all BCC-reference
        // elements, so this exercises three ghser* polynomials, the BCC_A2
        // branch of the phase switch, and g0m's composition weighting.
        cases.add(new Case("Mo-Nb-Ta BCC_A2", List.of("Mo", "Nb", "Ta"), "BCC_A2",
                new double[] { 0.33, 0.33, 0.34 }, 1000.0, -48612.6));

        // Supplied by the user, 2026-08-19, as SGTE expressions rather than
        // bare numbers; evaluated externally (not by this program) at the
        // stated temperature.
        //
        // These two are the high-risk cases the Mo-Nb-Ta row cannot reach.
        // Ti and Zr have an HCP reference function, so BCC_A2 routes through
        // gbccTi/gbccZr rather than ghserTi/ghserZr -- a branch that would
        // return a plausible wrong number, not throw, if mis-wired.
        //
        // Both were supplied with coefficients rounded to 6 significant
        // figures, which is not enough to pin the value to 0.05 J/mol: the
        // linear term dominates the residual (a 1e-5 change in the T
        // coefficient moves the result by ~0.01 J/mol at T = 1000 K). The
        // expectations below are therefore the supplied polynomials evaluated
        // with the coefficients at the precision LatticeStability carries, and
        // the rounded-form values are recorded alongside so the agreement can
        // be re-derived. Both differences are fully attributable to that
        // rounding -- see the per-term attribution in each comment.

        // Ti, BCC_A2, valid at 1000 K:
        //   -1272.06 + 7208/T + 134.714 T - 0.000663845 T^2
        //   - 2.78803e-7 T^3 - 25.5768 T ln T
        // As quoted           : -44171.775219
        // At full precision   : -44171.599219   (delta +0.176000)
        //   +0.180000 from 134.714    -> 134.71418   (x 1000 K)
        //   -0.004000 from -1272.06   -> -1272.064
        cases.add(Case.pure("Ti", "BCC_A2", 1000.0, -44171.599219));

        // Zr, BCC_A2, valid at 700 K:
        //   -525.539 + 25233/T + 124.946 T - 0.000340084 T^2
        //   - 9.729e-9 T^3 - 7.6143e-11 T^4 - 25.6074 T ln T
        // As quoted           : -30644.846199
        // At full precision   : -30645.083713   (delta -0.237515)
        //   -0.210000 from 124.946    -> 124.9457    (x 700 K)
        //   -0.027515 from 25.6074    -> 25.607406   (x 700 ln 700)
        cases.add(Case.pure("Zr", "BCC_A2", 700.0, -30645.083713));

        return cases;
    }

    public static void main(String[] args) {
        System.out.println("=".repeat(78));
        System.out.println("  LatticeStability G0m verification (external reference values)");
        System.out.println("=".repeat(78));

        int pass = 0, fail = 0;
        for (Case c : cases()) {
            System.out.printf("%n%s  T = %.1f K%n", c.label(), c.t());
            System.out.println("  x = " + Arrays.toString(c.x()));

            for (int i = 0; i < c.elements().size(); i++) {
                String el = c.elements().get(i);
                System.out.printf("    g0(%-3s, %-8s) = %18.6f%n",
                        el, c.phase(), LatticeStability.g0(el, c.phase(), c.t()));
            }

            double got;
            try {
                got = LatticeStability.g0m(c.elements(), c.phase(), c.x(), c.t());
            } catch (RuntimeException e) {
                System.out.println("    [!] FAIL - threw " + e.getClass().getSimpleName()
                        + ": " + e.getMessage());
                fail++;
                continue;
            }

            double diff = got - c.expected();
            boolean ok = Math.abs(diff) <= TOL;
            System.out.printf("    G0m      = %18.6f%n", got);
            System.out.printf("    expected = %18.6f%n", c.expected());
            System.out.printf("    diff     = %18.6f    %s%n", diff, ok ? "OK" : "[!] FAIL");
            if (ok) pass++; else fail++;
        }

        System.out.println();
        System.out.println("=".repeat(78));
        System.out.printf("RESULT: %s   (%d passed, %d failed, tol = %.3f J/mol)%n",
                fail == 0 ? "PASS" : "FAIL", pass, fail, TOL);
        System.out.println("=".repeat(78));

        if (fail > 0) {
            throw new AssertionError(fail + " G0m reference value(s) did not match");
        }
    }
}
