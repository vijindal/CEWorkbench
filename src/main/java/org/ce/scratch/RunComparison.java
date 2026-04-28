package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.thermo.ThermodynamicWorkflow;
import org.ce.model.ThermodynamicResult;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace;

import java.util.function.Consumer;

public class RunComparison {
    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();
        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);

        Consumer<String> nullSink = s -> {};
        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow();

        // ── 4-component: Nb-Ti-V-Zr at (0.25, 0.25, 0.25, 0.25), T=1000K ──
        SystemId quatId = new SystemId("Nb-Ti-V-Zr", "BCC_A2", "T");
        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());
        ModelSession mcsSess4 = builder.build(quatId, EngineConfig.MCS, nullSink);
        ModelSession cvmSess4 = builder.build(quatId, EngineConfig.CVM, nullSink);

        double T = 1000.0;
        double[] x4 = {0.25, 0.25, 0.25, 0.25};

        ThermodynamicResult rCvm4 = thermoWorkflow.runCalculation(cvmSess4,
                new ThermodynamicWorkflow.Request(T, x4, Property.ENTHALPY, nullSink, null,
                        0, 0, 0, null));

        long t0 = System.currentTimeMillis();
        ThermodynamicResult rMcs4 = thermoWorkflow.runCalculation(mcsSess4,
                new ThermodynamicWorkflow.Request(T, x4, Property.ENTHALPY, nullSink, null,
                        6, 2000, 4000, null));
        long elapsed = System.currentTimeMillis() - t0;

        System.out.println("=== CVM vs MCS: Nb-Ti-V-Zr BCC_A2 T-model, T=1000K, x=(0.25,0.25,0.25,0.25) ===");
        System.out.println();
        System.out.printf("  CVM  H = %12.4f J/mol%n", rCvm4.enthalpy);
        System.out.printf("  MCS  H = %12.4f J/mol   (L=6, equil=2000, avg=4000, %d ms)%n", rMcs4.enthalpy, elapsed);
        System.out.printf("  |ΔH| = %12.4f J/mol%n", Math.abs(rCvm4.enthalpy - rMcs4.enthalpy));
        System.out.println();

        // CF comparison
        double[] cvmCFs = rCvm4.optimizedCFs;
        double[] mcsCFs = rMcs4.avgCFs;

        int ncf = (cvmCFs != null && mcsCFs != null) ? Math.min(cvmCFs.length, mcsCFs.length) : 0;
        if (ncf == 0) {
            System.out.println("  (CFs not available for comparison)");
        } else {
            double maxDelta = 0.0;
            System.out.printf("  %-5s  %-12s  %-12s  %-12s%n", "CF#", "CVM", "MCS", "|Δ|");
            System.out.println("  " + "-".repeat(48));
            for (int i = 0; i < ncf; i++) {
                double delta = Math.abs(cvmCFs[i] - mcsCFs[i]);
                maxDelta = Math.max(maxDelta, delta);
                System.out.printf("  %-5d  %+12.6f  %+12.6f  %12.6f%n", i, cvmCFs[i], mcsCFs[i], delta);
            }
            System.out.println("  " + "-".repeat(48));
            System.out.printf("  Max |ΔCF| = %.6f%n", maxDelta);
        }

        System.out.println();

        // ── Binary sanity check: Nb-Ti at (0.5, 0.5) ──
        SystemId binId = new SystemId("Nb-Ti", "BCC_A2", "T");
        ModelSession mcsSess2 = builder.build(binId, EngineConfig.MCS, nullSink);
        ModelSession cvmSess2 = builder.build(binId, EngineConfig.CVM, nullSink);
        double[] x2 = {0.5, 0.5};

        ThermodynamicResult rCvm2 = thermoWorkflow.runCalculation(cvmSess2,
                new ThermodynamicWorkflow.Request(T, x2, Property.ENTHALPY, nullSink, null,
                        0, 0, 0, null));
        ThermodynamicResult rMcs2 = thermoWorkflow.runCalculation(mcsSess2,
                new ThermodynamicWorkflow.Request(T, x2, Property.ENTHALPY, nullSink, null,
                        6, 2000, 4000, null));

        System.out.println("=== Sanity check: Nb-Ti, x=(0.5,0.5) ===");
        System.out.printf("  CVM  H = %12.4f J/mol%n", rCvm2.enthalpy);
        System.out.printf("  MCS  H = %12.4f J/mol%n", rMcs2.enthalpy);
        System.out.printf("  |ΔH| = %12.4f J/mol%n", Math.abs(rCvm2.enthalpy - rMcs2.enthalpy));
    }
}
