package org.ce.scratch;

import org.ce.model.equilibrium.HillertSolver.MassBalanceReport;

import java.util.Arrays;

/**
 * Regression gate for {@link MassBalanceReport}'s defensive-copy contract
 * (STEP 5, PART 1).
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.MassBalanceReportImmutability
 * </pre>
 *
 * <p>{@code MassBalanceReport} is a record with two {@code double[]} components
 * ({@code targetOverall}, {@code calculatedOverall}). Records auto-generate an
 * accessor that returns the field by reference, and the canonical constructor
 * stores its argument by reference -- so before this fix a caller could mutate
 * the report's inventory vectors either through the array it passed in or
 * through an array the accessor handed back. The record now clones on the way
 * in (compact canonical constructor) and on the way out (overridden
 * accessors).</p>
 *
 * <p>Checks: (1) mutating the arrays passed to the constructor does not change
 * the report; (2) mutating an array returned by an accessor does not change the
 * report; (3) two accessor calls return independent arrays; (4) a {@code null}
 * component is preserved as {@code null} (a rejected-problem report can carry
 * one); (5) the scalar components are untouched.</p>
 */
public final class MassBalanceReportImmutability {

    private static int failures = 0;

    public static void main(String[] args) {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(78));
        System.out.println("  MassBalanceReport defensive-copy contract");
        System.out.println("=".repeat(78));

        checkConstructorInputIsolation();
        checkAccessorOutputIsolation();
        checkAccessorsReturnIndependentArrays();
        checkNullPreserved();
        checkScalarsIntact();

        System.out.println("\n" + "=".repeat(78));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(78));
        if (failures > 0) {
            throw new AssertionError(failures + " immutability checks failed");
        }
    }

    private static void checkConstructorInputIsolation() {
        System.out.println("\n--- Constructor input isolation ---");

        double[] target = { 0.25, 0.75 };
        double[] calc = { 0.30, 0.70 };
        MassBalanceReport r = new MassBalanceReport(target, calc, 1.5e-3, 2.0e-6, 4.0e-3);

        double[] targetSnapshot = target.clone();
        double[] calcSnapshot = calc.clone();

        // Mutate the caller's original arrays after construction.
        Arrays.fill(target, -999.0);
        Arrays.fill(calc, Double.NaN);

        check("targetOverall() unaffected by mutating the constructor input",
                Arrays.equals(r.targetOverall(), targetSnapshot),
                Arrays.toString(r.targetOverall()));
        check("calculatedOverall() unaffected by mutating the constructor input",
                Arrays.equals(r.calculatedOverall(), calcSnapshot),
                Arrays.toString(r.calculatedOverall()));
    }

    private static void checkAccessorOutputIsolation() {
        System.out.println("\n--- Accessor output isolation ---");

        MassBalanceReport r = new MassBalanceReport(
                new double[] { 0.4, 0.6 }, new double[] { 0.5, 0.5 }, 1.0e-3, 1.0e-6, 2.0e-3);

        double[] pristineTarget = r.targetOverall().clone();
        double[] pristineCalc = r.calculatedOverall().clone();

        // Scribble on what the accessors return.
        double[] handed1 = r.targetOverall();
        Arrays.fill(handed1, 12345.0);
        double[] handed2 = r.calculatedOverall();
        Arrays.fill(handed2, -1.0);

        check("targetOverall() unchanged after scribbling on a prior return value",
                Arrays.equals(r.targetOverall(), pristineTarget),
                Arrays.toString(r.targetOverall()));
        check("calculatedOverall() unchanged after scribbling on a prior return value",
                Arrays.equals(r.calculatedOverall(), pristineCalc),
                Arrays.toString(r.calculatedOverall()));
    }

    private static void checkAccessorsReturnIndependentArrays() {
        System.out.println("\n--- Accessors return independent arrays ---");

        MassBalanceReport r = new MassBalanceReport(
                new double[] { 0.1, 0.2, 0.7 }, new double[] { 0.2, 0.2, 0.6 },
                3.0e-3, 3.0e-6, 5.0e-3);

        check("two targetOverall() calls are distinct objects",
                r.targetOverall() != r.targetOverall(), "same reference");
        check("two calculatedOverall() calls are distinct objects",
                r.calculatedOverall() != r.calculatedOverall(), "same reference");
    }

    private static void checkNullPreserved() {
        System.out.println("\n--- Null component preserved ---");

        MassBalanceReport r = new MassBalanceReport(
                null, null, Double.NaN, Double.NaN, Double.NaN);
        check("null targetOverall stays null", r.targetOverall() == null,
                String.valueOf((Object) r.targetOverall()));
        check("null calculatedOverall stays null", r.calculatedOverall() == null,
                String.valueOf((Object) r.calculatedOverall()));
    }

    private static void checkScalarsIntact() {
        System.out.println("\n--- Scalar components intact ---");

        MassBalanceReport r = new MassBalanceReport(
                new double[] { 0.5, 0.5 }, new double[] { 0.5, 0.5 },
                1.25e-4, 6.25e-7, 9.0e-4);
        check("maxAbsResidual preserved", r.maxAbsResidual() == 1.25e-4,
                String.valueOf(r.maxAbsResidual()));
        check("maxRelResidual preserved", r.maxRelResidual() == 6.25e-7,
                String.valueOf(r.maxRelResidual()));
        check("residualBeforeLastStep preserved", r.residualBeforeLastStep() == 9.0e-4,
                String.valueOf(r.residualBeforeLastStep()));
    }

    private static void check(String label, boolean ok, String detail) {
        if (ok) {
            System.out.printf("    %-64s OK%n", label);
        } else {
            failures++;
            System.out.printf("    %-64s [!] FAIL  %s%n", label, detail);
        }
    }

    private MassBalanceReportImmutability() {
    }
}
