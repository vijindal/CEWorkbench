package org.ce.calculation.engine.cvm;

import org.ce.model.cluster.LinearAlgebra;
import org.ce.calculation.engine.ProgressEvent;
import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMGibbsModel.ModelResult;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CancellationException;

/**
 * Numerical solver for CVM equilibrium calculations using Newton-Raphson.
 */
public class CVMSolver {

    private static final int MAX_ITER = 400;
    private static final double TOLX = 1.0e-12;

    public static final class EquilibriumResult {
        public final double[] u;
        public final ModelResult modelValues;
        public final int iterations;
        public final double finalGradientNorm;
        public final boolean converged;
        public final List<IterationSnapshot> trace;

        public EquilibriumResult(double[] u, ModelResult modelValues, int iterations, 
                                 double finalGradientNorm, boolean converged, 
                                 List<IterationSnapshot> trace) {
            this.u = u;
            this.modelValues = modelValues;
            this.iterations = iterations;
            this.finalGradientNorm = finalGradientNorm;
            this.converged = converged;
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
            this.G = G; this.H = H; this.S = S;
            this.gradientNorm = gradientNorm;
            this.u = u;
        }
    }

    /**
     * Minimizes the Gibbs free energy for the given model and conditions.
     */
    public EquilibriumResult minimize(CVMGibbsModel model, double[] moleFractions, 
                                     double temperature, double[] eci, double tolerance,
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
            if (Thread.currentThread().isInterrupted()) throw new CancellationException();

            // Evaluate physics from model
            current = model.evaluate(u, moleFractions, temperature, eci);
            
            // Calculate gradient norm (L1 norm for simplicity and consistency with legacy)
            errf = 0;
            for (double g : current.Gu) errf += Math.abs(g);
            
            trace.add(new IterationSnapshot(its, current.G, current.H, current.S, errf, u.clone()));

            // Progress reporting
            if (progressSink != null) {
                progressSink.accept(String.format("  iter %3d  |∇G| = %.3e  G = %11.4f  H = %11.4f  S = %9.6f",
                        its, errf, current.G, current.H, current.S));
            }
            if (eventSink != null) {
                eventSink.accept(new ProgressEvent.CvmIteration(its, current.G, errf, current.H, current.S, u));
            }

            // Convergence check
            if (errf <= tolerance) {
                return new EquilibriumResult(u, current, its, errf, true, trace);
            }

            try {
                // Newton step: Guu * p = -Gu
                double[] negGu = new double[n];
                for (int i = 0; i < n; i++) negGu[i] = -current.Gu[i];
                
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
                    current = model.evaluate(u, moleFractions, temperature, eci);
                    return new EquilibriumResult(u, current, its, errf, true, trace);
                }

            } catch (Exception e) {
                return new EquilibriumResult(u, current, its, errf, false, trace);
            }
        }

        return new EquilibriumResult(u, current, MAX_ITER, errf, false, trace);
    }

    /** Overload for cases where progress reporting is not needed. */
    public EquilibriumResult minimize(CVMGibbsModel model, double[] moleFractions, 
                                     double temperature, double[] eci, double tolerance) {
        return minimize(model, moleFractions, temperature, eci, tolerance, null, null);
    }
}
