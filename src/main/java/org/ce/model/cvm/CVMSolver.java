package org.ce.model.cvm;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.ProgressEvent;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMGibbsModel.ModelResult;
import org.ce.model.cvm.CVMEquilibriumState;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CancellationException;

/**
 * Numerical solver for CVM equilibrium calculations using Newton-Raphson.
 */
public class CVMSolver {

    private static final int MAX_ITER = 20;
    private static final double TOLX = 1.0e-12;
    private static final double R_GAS = 8.3144598; // J/(mol·K)

    /**
     * Result of CVM equilibrium calculation.
     *
     * <p>
     * Separates thermodynamic contract (state) from solver metadata (convergence
     * info, trace).
     * The calculation layer should check {@link #converged} before using
     * {@link #state}.
     * </p>
     */
    public static final class EquilibriumResult {

        /** The thermodynamic equilibrium state. Use this in the calculation layer. */
        public final CVMEquilibriumState state;

        /** Solver metadata: convergence flag. Check this before using state. */
        public final boolean converged;

        /** Solver metadata: iteration count. */
        public final int iterations;

        /** Solver metadata: final gradient norm ||∇G||. */
        public final double finalGradientNorm;

        /** Debug trace: iteration snapshots. Only populated for diagnostics. */
        public final List<IterationSnapshot> trace;

        /**
         * @param state             thermodynamic equilibrium state (u, G/H/S)
         * @param converged         convergence flag
         * @param iterations        iteration count at convergence/failure
         * @param finalGradientNorm final ||∇G|| norm
         * @param trace             iteration snapshots for debugging (may be empty)
         */
        public EquilibriumResult(CVMEquilibriumState state, boolean converged, int iterations,
                double finalGradientNorm, List<IterationSnapshot> trace) {
            this.state = state;
            this.converged = converged;
            this.iterations = iterations;
            this.finalGradientNorm = finalGradientNorm;
            this.trace = trace;
        }
    }

    public static final class IterationSnapshot {
        public final int iteration;
        public final double G, H, S;
        public final double gradientNorm;
        public final double[] u;

        public IterationSnapshot(int iteration, double G, double H, double S, double gradientNorm, double[] u) {
            this.iteration = iteration;
            this.G = G;
            this.H = H;
            this.S = S;
            this.gradientNorm = gradientNorm;
            this.u = u;
        }
    }

    /**
     * Minimizes the Gibbs free energy for the given model and conditions.
     */
    public EquilibriumResult minimize(CVMGibbsModel model, double[] moleFractions,
            double temperature, double tolerance,
            java.util.function.Consumer<String> progressSink,
            java.util.function.Consumer<ProgressEvent> eventSink) {

        if (eventSink != null) {
            eventSink.accept(new ProgressEvent.EngineStart("CVM", 0));
        }
        int n = model.getNcf();
        double[] u = model.computeRandomCFs(moleFractions);
        List<IterationSnapshot> trace = new ArrayList<>();

        double errf = 0;
        ModelResult current = null;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException();

            // Evaluate physics from model
            current = model.evaluate(u, moleFractions, temperature);

            // Calculate gradient norm (L1 norm for simplicity and consistency with legacy)
            errf = 0;
            for (double g : current.Gu)
                errf += Math.abs(g);

            trace.add(new IterationSnapshot(its, current.G, current.H, current.S, errf, u.clone()));

            // Progress reporting
            if (progressSink != null) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("    iter %3d  |∇G| = %.3e  G = %11.4f  H = %11.4f  S = %9.6f",
                        its, errf, current.G, current.H, current.S));
                if (current.cfs != null && current.cfs.length > 0) {
                    sb.append("  CFs: [");
                    for (int i = 0; i < current.cfs.length; i++) {
                        if (i > 0) sb.append(", ");
                        sb.append(String.format("%.4f", current.cfs[i]));
                    }
                    sb.append("]");
                }
                progressSink.accept(sb.toString());
            }
            if (eventSink != null) {
                eventSink.accept(new ProgressEvent.CvmIteration(its, current.G, errf, current.H, current.S, u));
            }

            // Convergence check
            if (errf <= tolerance) {
                CVMEquilibriumState state = new CVMEquilibriumState(u, current, temperature, R_GAS);
                return new EquilibriumResult(state, true, its, errf, trace);
            }

            try {
                // Newton step: Guu * p = -Gu
                double[] negGu = new double[n];
                for (int i = 0; i < n; i++)
                    negGu[i] = -current.Gu[i];

                double[] p = LinearAlgebra.solve(current.Guu, negGu);

                // Step size limiting (physical bounds)
                double alpha = model.calculateStepLimit(u, p, moleFractions);

                double errx = 0;
                for (int i = 0; i < n; i++) {
                    double delta = alpha * p[i];
                    u[i] += delta;
                    errx += Math.abs(delta);
                }

                // X-convergence check
                if (errx <= TOLX) {
                    current = model.evaluate(u, moleFractions, temperature);
                    CVMEquilibriumState state = new CVMEquilibriumState(u, current, temperature, R_GAS);
                    return new EquilibriumResult(state, true, its, errf, trace);
                }

            } catch (Exception e) {
                CVMEquilibriumState state = (current != null) ? new CVMEquilibriumState(u, current, temperature, R_GAS)
                        : null;
                return new EquilibriumResult(state, false, its, errf, trace);
            }
        }

        CVMEquilibriumState state = (current != null) ? new CVMEquilibriumState(u, current, temperature, R_GAS) : null;
        return new EquilibriumResult(state, false, MAX_ITER, errf, trace);
    }

    /** Overload for cases where progress reporting is not needed. */
    public EquilibriumResult minimize(CVMGibbsModel model, double[] moleFractions,
            double temperature, double tolerance) {
        return minimize(model, moleFractions, temperature, tolerance, null, null);
    }
}
