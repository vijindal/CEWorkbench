package org.ce.calculation.workflow.thermo;

import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Static scan utilities for thermodynamic parameter sweeps.
 *
 * <p>All ANALYSIS-mode calculations funnel through {@link #scan}, which executes
 * a 2D grid over any two varying parameters (T or Xi). A single-point calculation
 * is a 1×1 grid; a line scan is N×1 or 1×M — no separate code path is needed.</p>
 */
public class ScanWorkflows {

    private ScanWorkflows() {}

    /**
     * Defines one varying parameter dimension for a scan.
     * When {@code start == end} the loop runs exactly once (fixed dimension).
     */
    public record Varying(String label, boolean isTemp, int compIndex,
                          double start, double end, double step) {}

    /**
     * Derives the full n-component composition array from n-1 independent components.
     *
     * @param xIndep independent components [x2, x3, ..., xn]
     * @param session active session providing component count
     * @return full normalized array [x1, x2, ..., xn]
     */
    public static double[] deriveComposition(double[] xIndep, ModelSession session) {
        int nTotal = session.numComponents();
        double[] full = new double[nTotal];
        double sumIndep = 0;
        for (int i = 0; i < xIndep.length && i + 1 < nTotal; i++) {
            full[i + 1] = xIndep[i];
            sumIndep += xIndep[i];
        }
        full[0] = 1.0 - sumIndep;
        return full;
    }

    /**
     * Universal 2D grid scan over any two varying parameters (T or Xi).
     *
     * <p>Pass a sentinel {@code Varying} with {@code start == end} for a fixed
     * dimension — the loop runs once, producing a 1-row or 1-column grid.</p>
     *
     * @param workflow     calculation engine
     * @param session      pre-built model session
     * @param var1         outer (row) dimension
     * @param var2         inner (column) dimension
     * @param baseIndepComp base independent composition [x2, x3, ...]
     * @param baseT        base temperature (K)
     * @param property     thermodynamic property to calculate
     * @param mcsL         MCS lattice size
     * @param mcsNEquil    MCS equilibration sweeps
     * @param mcsNAvg      MCS averaging sweeps
     * @param strSink      optional text progress sink
     * @param eventSink    optional structured progress sink
     * @return row-major grid of results
     */
    public static List<List<ThermodynamicResult>> scan(
            ThermodynamicWorkflow workflow,
            ModelSession session,
            Varying var1,
            Varying var2,
            double[] baseIndepComp,
            double baseT,
            org.ce.calculation.CalculationDescriptor.Property property,
            int mcsL, int mcsNEquil, int mcsNAvg,
            Consumer<String> strSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        List<List<ThermodynamicResult>> grid = new ArrayList<>();

        for (double v1 = var1.start();
             var1.step() > 0 ? v1 <= var1.end() + 1e-9 : v1 >= var1.end() - 1e-9;
             v1 += Math.abs(var1.step()) > 1e-10 ? var1.step() : 1.0) {

            if (strSink != null && Math.abs(var1.step()) > 1e-10) {
                strSink.accept(String.format("Scanning %s = %.4g", var1.label(), v1));
            }

            List<ThermodynamicResult> row = new ArrayList<>();

            for (double v2 = var2.start();
                 var2.step() > 0 ? v2 <= var2.end() + 1e-9 : v2 >= var2.end() - 1e-9;
                 v2 += Math.abs(var2.step()) > 1e-10 ? var2.step() : 1.0) {

                double T = baseT;
                double[] xIndep = (baseIndepComp != null)
                        ? baseIndepComp.clone()
                        : new double[session.numComponents() - 1];

                if (var1.isTemp()) T = v1; else if (var1.compIndex() >= 0) xIndep[var1.compIndex()] = v1;
                if (var2.isTemp()) T = v2; else if (var2.compIndex() >= 0) xIndep[var2.compIndex()] = v2;

                row.add(workflow.runCalculation(session,
                        new ThermodynamicRequest(T, deriveComposition(xIndep, session),
                                property, strSink, eventSink, mcsL, mcsNEquil, mcsNAvg)));

                if (Math.abs(var2.step()) < 1e-10) break;
            }

            grid.add(row);
            if (Math.abs(var1.step()) < 1e-10) break;
        }
        return grid;
    }

    /**
     * Finite-size scaling scan: runs MCS at several supercell sizes L and
     * extrapolates to the thermodynamic limit.
     */
    public static ThermodynamicResult finiteSizeScan(
            ThermodynamicWorkflow workflow,
            ModelSession session,
            double temperature,
            double[] composition,
            org.ce.calculation.CalculationDescriptor.Property property,
            int nEquil,
            int nAvg,
            Consumer<String> strSink,
            Consumer<ProgressEvent> eventSink) throws Exception {

        int[] Ls = {12, 16, 24};
        ThermodynamicResult[] results = new ThermodynamicResult[Ls.length];

        for (int i = 0; i < Ls.length; i++) {
            if (strSink != null) strSink.accept("FSS: Running L=" + Ls[i]);
            results[i] = workflow.runCalculation(session,
                    new ThermodynamicRequest(temperature, composition, property,
                            strSink, eventSink, Ls[i], nEquil, nAvg));
        }

        ThermodynamicResult last = results[results.length - 1];
        return new ThermodynamicResult(
                results[0].temperature,
                results[0].composition,
                results[0].gibbsEnergy,
                last.enthalpy,
                Double.NaN,
                last.stdEnthalpy,
                last.heatCapacity,
                null,
                last.avgCFs,
                last.stdCFs
        );
    }
}
