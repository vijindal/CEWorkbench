package org.ce.model.cvm;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.model.ProgressEvent;
import org.ce.model.cvm.CVMGibbsModel.ModelResult;

import java.util.concurrent.CancellationException;

/**
 * Numerical solver for CVM equilibrium calculations using Newton-Raphson.
 */
public class CVMSolver {

    private static final int MAX_ITER = 20;
    private static final double TOLX = 1.0e-12;

    /**
     * Result of CVM equilibrium calculation.
     */
    public static final class EquilibriumResult {

        /** Physics values at the equilibrium point. */
        public final ModelResult modelResult;

        /** Equilibrium non-point CVCF correlation functions (length = ncf). */
        public final double[] u;

        /** Convergence flag. Check before using modelResult. */
        public final boolean converged;

        /** Iteration count at convergence or failure. */
        public final int iterations;

        /** Final gradient norm ||∇G|| at exit. */
        public final double finalGradientNorm;

        public EquilibriumResult(ModelResult modelResult, double[] u, boolean converged,
                int iterations, double finalGradientNorm) {
            this.modelResult        = modelResult;
            this.u                  = u;
            this.converged          = converged;
            this.iterations         = iterations;
            this.finalGradientNorm  = finalGradientNorm;
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

        double errf = 0;
        ModelResult current = null;

        for (int its = 0; its < MAX_ITER; its++) {
            if (Thread.currentThread().isInterrupted())
                throw new CancellationException();

            // Evaluate physics from model
            current = model.evaluate(u, moleFractions, temperature);

            // Calculate gradient norm (L1 norm)
            errf = 0;
            for (double g : current.Gu)
                errf += Math.abs(g);

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
                return new EquilibriumResult(current, u.clone(), true, its, errf);
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
                    return new EquilibriumResult(current, u.clone(), true, its, errf);
                }

            } catch (Exception e) {
                return new EquilibriumResult(current, u.clone(), false, its, errf);
            }
        }

        return new EquilibriumResult(current, u.clone(), false, MAX_ITER, errf);
    }

    /** Overload for cases where progress reporting is not needed. */
    public EquilibriumResult minimize(CVMGibbsModel model, double[] moleFractions,
            double temperature, double tolerance) {
        return minimize(model, moleFractions, temperature, tolerance, null, null);
    }
}
