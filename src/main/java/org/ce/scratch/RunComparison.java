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
        
        SystemId systemId = new SystemId("Nb-Ti", "BCC_A2", "T");
        
        // Suppress all verbose logs for clean timing
        Consumer<String> nullSink = s -> {};

        ModelSession.Builder builder = new ModelSession.Builder(context.getHamiltonianStore());
        ModelSession mcsSession = builder.build(systemId, EngineConfig.MCS, nullSink);
        ModelSession cvmSession = builder.build(systemId, EngineConfig.CVM, nullSink);

        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow();

        double temp = 1000.0;
        double[][] compositions = {
            {0.1, 0.9},
            {0.3, 0.7},
            {0.5, 0.5},
            {0.7, 0.3},
            {0.9, 0.1}
        };

        System.out.println("=== MCS PERFORMANCE BENCHMARK (Nb-Ti, T=1000K, L=6, equil=100, avg=500) ===");
        System.out.printf(" %-12s | %-12s | %-12s | %-10s | %-12s%n",
                "Composition", "CVM H", "MCS H", "MCS Time", "ΔH");
        System.out.println("--------------------------------------------------------------------------");

        long totalMcsTime = 0;
        for (double[] x : compositions) {
            // Run CVM
            ThermodynamicWorkflow.Request reqCvm = new ThermodynamicWorkflow.Request(
                    temp, x, Property.ENTHALPY, nullSink, null,
                    0, 0, 0, null);
            ThermodynamicResult rCvm = thermoWorkflow.runCalculation(cvmSession, reqCvm);

            // Run MCS with timing
            ThermodynamicWorkflow.Request reqMcs = new ThermodynamicWorkflow.Request(
                    temp, x, Property.ENTHALPY, nullSink, null,
                    6, 100, 500, null);
            long t0 = System.currentTimeMillis();
            ThermodynamicResult rMcs = thermoWorkflow.runCalculation(mcsSession, reqMcs);
            long elapsed = System.currentTimeMillis() - t0;
            totalMcsTime += elapsed;

            double deltaH = Math.abs(rCvm.enthalpy - rMcs.enthalpy);
            System.out.printf(" [%.1f, %.1f]   | %10.2f   | %10.2f   | %6d ms  | %8.2f%n",
                    x[0], x[1], rCvm.enthalpy, rMcs.enthalpy, elapsed, deltaH);
        }
        System.out.println("--------------------------------------------------------------------------");
        System.out.printf(" TOTAL MCS TIME: %d ms (avg: %d ms/composition)%n",
                totalMcsTime, totalMcsTime / compositions.length);
        System.out.println("==========================================================================");
    }
}
