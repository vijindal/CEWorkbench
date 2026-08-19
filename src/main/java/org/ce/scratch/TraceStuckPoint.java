package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.List;
import java.util.function.Consumer;

/** Diagnostic: trace the Newton-Raphson iteration path at a known-stuck near-edge composition. */
public class TraceStuckPoint {
    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);

        Consumer<String> nullSink = s -> {};
        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow();

        SystemId id = new SystemId("Nb-Ti-V", "BCC_A2", "T");
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());
        ModelSession sess = builder.build(id, EngineConfig.CVM, nullSink);

        double T = 1273.0;
        // matches the Mathematica reference point: Nb (dilute) = 1/60, Ti = 1/3, V = 0.65
        double[] x = {1.0 / 60, 0.333333, 0.65};

        List<String> cfNames = CvCfBasis.getNonPointCfNames("BCC_A2", "T", 3);
        int iv3ab = cfNames.indexOf("v3AB");
        int iv21ab = cfNames.indexOf("v21AB");
        System.out.println("cfNames=" + cfNames);

        Consumer<ProgressEvent> eventSink = ev -> {
            if (ev instanceof ProgressEvent.CvmIteration it) {
                String v3ab = (it.cfs != null && iv3ab >= 0) ? String.format("%.6f", it.cfs[iv3ab]) : "?";
                String v21ab = (it.cfs != null && iv21ab >= 0) ? String.format("%.6f", it.cfs[iv21ab]) : "?";
                System.out.printf("it=%3d  G=%14.4f  gradNorm=%16.6f  v3AB=%s  v21AB=%s%n",
                        it.iteration, it.gibbsEnergy, it.gradientNorm, v3ab, v21ab);
            }
        };

        ThermodynamicResult r = thermoWorkflow.runCalculation(sess,
                new ThermodynamicWorkflow.Request(T, x, Property.GIBBS_ENERGY, nullSink, eventSink,
                        0, 0, 0, null));

        System.out.println();
        System.out.println("converged=" + r.converged + " iterations=" + r.iterations
                + " finalGradientNorm=" + r.finalGradientNorm + " G=" + r.gibbsEnergy);
    }
}
