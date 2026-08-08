package org.ce.model.mcs;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.ModelSpecifications;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.Conditions;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.CvCfBasis;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Compares CVM (fast, analytic, tetrahedron-truncated) against MCS (slower,
 * statistical, not limited by the cluster approximation's entropy expression)
 * for multicomponent systems, via the same public API real callers use
 * (CalculationService.calculate) — not a reimplementation of either engine.
 *
 * <p>Goal: check that MCS enthalpy tracks CVM's trend with temperature and
 * composition (both should show the same qualitative behavior — enthalpy
 * approaching 0 at pure elements, smooth variation with T and x), even though
 * MCS is expected to be quantitatively different (and more physically correct,
 * per USER_GUIDE) since it isn't limited by the tetrahedron entropy truncation.
 *
 * <p>Uses small L and modest sweep counts for MCS to keep runtime reasonable —
 * this is a trend check, not a converged-value benchmark.
 */
public class CvmVsMcsComparisonTest {

    public static void main(String[] args) throws Exception {
        CEWorkbenchContext appCtx = new CEWorkbenchContext();
        CalculationService service = appCtx.getCalculationService();

        CalculationService.McsParams mcsParams = new CalculationService.McsParams(4, 80, 150);

        System.out.println("=== CVM vs MCS: Temperature sweep + CF comparison, Nb-Ti-V equiatomic ===");
        temperatureSweep(service, "Nb-Ti-V", Map.of("Ti", 1.0 / 3.0, "V", 1.0 / 3.0),
                new double[] { 800, 1200 }, mcsParams);

        System.out.println();
        System.out.println("=== CVM vs MCS: Composition sweep + CF comparison, Nb-Ti-V at T=1000K (V fixed 0.2, varying Ti) ===");
        compositionSweep(service, "Nb-Ti-V", "Ti", "V", 0.2,
                new double[] { 0.1, 0.4, 0.7 }, 1000.0, mcsParams);
    }

    private static void temperatureSweep(CalculationService service, String elements,
            Map<String, Double> composition, double[] temps, CalculationService.McsParams mcsParams) throws Exception {

        ModelSpecifications cvmSpecs = new ModelSpecifications(elements, "BCC_A2", "T", EngineConfig.CVM);
        ModelSpecifications mcsSpecs = new ModelSpecifications(elements, "BCC_A2", "T", EngineConfig.MCS);
        ModelSession cvmSession = service.getOrBuildSession(cvmSpecs, null);
        ModelSession mcsSession = new ModelSession.Builder(appCtxHamStore(service))
                .build(mcsSession_systemId(elements), EngineConfig.MCS, cvmSession.cecEntry, null);

        List<String> cfNames = CvCfBasis.getNonPointCfNames("BCC_A2", "T", cvmSession.numComponents());

        System.out.printf("%-8s  %14s  %14s  %10s  %8s%n", "T (K)", "H_CVM", "H_MCS", "diff", "conv?");
        Double prevHCvm = null, prevHMcs = null;
        for (double T : temps) {
            Conditions cond = new Conditions(T, composition);
            ThermodynamicResult rCvm = service.calculate(cvmSession, cond, Property.GIBBS_ENERGY, null, null);
            ThermodynamicResult rMcs = service.calculate(mcsSession, cond, Property.ENTHALPY, mcsParams, null, null);

            String trendCvm = prevHCvm == null ? "" : (rCvm.enthalpy > prevHCvm ? "up" : "down");
            String trendMcs = prevHMcs == null ? "" : (rMcs.enthalpy > prevHMcs ? "up" : "down");
            String trendFlag = (!trendCvm.isEmpty() && !trendCvm.equals(trendMcs)) ? " <-- TREND MISMATCH" : "";

            System.out.printf("%-8.1f  %14.4f  %14.4f  %10.4f  %8s  %s/%s%s%n",
                    T, rCvm.enthalpy, rMcs.enthalpy, rMcs.enthalpy - rCvm.enthalpy,
                    rCvm.converged, trendCvm, trendMcs, trendFlag);
            printCfComparison(cfNames, rCvm.optimizedCFs, rMcs.avgCFs);

            prevHCvm = rCvm.enthalpy;
            prevHMcs = rMcs.enthalpy;
        }
    }

    /** Prints CVCF non-point correlation functions side by side, CVM (equilibrium) vs MCS (sampled). */
    private static void printCfComparison(List<String> cfNames, double[] cvmCfs, double[] mcsCfs) {
        if (cvmCfs == null || mcsCfs == null) {
            System.out.println("    [CF comparison unavailable: optimizedCFs or avgCFs is null]");
            return;
        }
        int n = Math.min(cfNames.size(), Math.min(cvmCfs.length, mcsCfs.length));
        System.out.printf("    %-10s  %10s  %10s  %10s%n", "CF", "CVM", "MCS", "diff");
        for (int i = 0; i < n; i++) {
            double diff = mcsCfs[i] - cvmCfs[i];
            System.out.printf("    %-10s  %10.5f  %10.5f  %10.5f%n", cfNames.get(i), cvmCfs[i], mcsCfs[i], diff);
        }
    }

    private static void compositionSweep(CalculationService service, String elements,
            String varyEl, String fixedEl, double fixedVal, double[] varyVals, double T,
            CalculationService.McsParams mcsParams) throws Exception {

        ModelSpecifications cvmSpecs = new ModelSpecifications(elements, "BCC_A2", "T", EngineConfig.CVM);
        ModelSession cvmSession = service.getOrBuildSession(cvmSpecs, null);
        ModelSession mcsSession = new ModelSession.Builder(appCtxHamStore(service))
                .build(mcsSession_systemId(elements), EngineConfig.MCS, cvmSession.cecEntry, null);

        List<String> cfNames = CvCfBasis.getNonPointCfNames("BCC_A2", "T", cvmSession.numComponents());

        System.out.printf("%-8s  %14s  %14s  %10s  %8s%n", varyEl + "=x", "H_CVM", "H_MCS", "diff", "conv?");
        for (double x : varyVals) {
            Map<String, Double> comp = new LinkedHashMap<>();
            comp.put(varyEl, x);
            comp.put(fixedEl, fixedVal);
            Conditions cond = new Conditions(T, comp);
            ThermodynamicResult rCvm = service.calculate(cvmSession, cond, Property.GIBBS_ENERGY, null, null);
            ThermodynamicResult rMcs = service.calculate(mcsSession, cond, Property.ENTHALPY, mcsParams, null, null);

            System.out.printf("%-8.2f  %14.4f  %14.4f  %10.4f  %8s%n",
                    x, rCvm.enthalpy, rMcs.enthalpy, rMcs.enthalpy - rCvm.enthalpy, rCvm.converged);
            printCfComparison(cfNames, rCvm.optimizedCFs, rMcs.avgCFs);
        }
    }

    // Small helpers to avoid duplicating session-building boilerplate — reach into
    // CEWorkbenchContext indirectly via the already-built CVM session's own store.
    private static org.ce.model.storage.DataStore.HamiltonianStore appCtxHamStore(CalculationService service) {
        // MCS sessions here reuse the CVM session's already-validated CECEntry directly
        // (build(systemId, engine, cecEntry, sink) bypasses the store entirely), so the
        // store instance itself is never touched — but the Builder constructor requires one.
        return new org.ce.model.storage.DataStore.HamiltonianStore(new org.ce.model.storage.Workspace());
    }

    private static org.ce.model.storage.Workspace.SystemId mcsSession_systemId(String elements) {
        return new org.ce.model.storage.Workspace.SystemId(elements, "BCC_A2", "T");
    }
}
