package org.ce.domain.engine;

import org.ce.domain.cluster.ClusCoordListResult;
import org.ce.domain.engine.mcs.MCResult;
import org.ce.domain.engine.mcs.MCSRunner;
import org.ce.domain.engine.mcs.MCSUpdate;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.domain.result.EquilibriumState;

import java.util.function.Consumer;

/**
 * Thermodynamic engine implementing Monte Carlo simulation (canonical ensemble).
 *
 * <p>Default parameters: L=4 supercell (128 BCC sites), 1000 equilibration sweeps,
 * 2000 averaging sweeps. ECI = a + b*T from CECEntry.</p>
 */
public class MCSEngine implements ThermodynamicEngine {

    private static final double GAS_CONSTANT      = 8.314;  // J/(mol·K)

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        ClusCoordListResult clusterData =
                input.clusterData.getDisorderedClusterResult().getDisClusterData();

        double[] eci = extractEci(input.cec.cecTerms, input.temperature);

        int L           = input.mcsL;
        int nEquil      = input.mcsNEquil;
        int nAvg        = input.mcsNAvg;
        int totalSweeps = nEquil + nAvg;
        int N           = 2 * L * L * L;  // BCC sites

        MCSRunner.Builder builder = MCSRunner.builder()
                .clusterData(clusterData)
                .eci(eci)
                .numComp(input.composition.length)
                .T(input.temperature)
                .composition(input.composition)
                .nEquil(nEquil)
                .nAvg(nAvg)
                .L(L)
                .seed(System.currentTimeMillis())
                .R(GAS_CONSTANT);

        Consumer<String> strSink = input.progressSink;
        Consumer<ProgressEvent> evtSink = input.eventSink;

        if (strSink != null || evtSink != null) {
            if (strSink != null) {
                strSink.accept(String.format("  MCS started: L=%d (%d sites), %d equil + %d avg sweeps",
                        L, N, nEquil, nAvg));
            }
            // EngineStart emitted here — before build().run() — so the chart clears before sweeps arrive
            if (evtSink != null) {
                evtSink.accept(new ProgressEvent.EngineStart("MCS", totalSweeps));
            }
            final int sitesCount = N;
            final int sweepCount = totalSweeps;
            builder.updateListener(mcUpdate -> {
                if (strSink != null && (mcUpdate.getStep() % 100 == 0 || mcUpdate.getStep() == sweepCount)) {
                    strSink.accept(String.format("  [%-5s] sweep %4d/%-4d  <E>/site=%9.5f  accept=%5.1f%%",
                            mcUpdate.getPhase(),
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate() * 100));
                }
                if (evtSink != null) {
                    evtSink.accept(new ProgressEvent.McSweep(
                            mcUpdate.getStep(), sweepCount,
                            mcUpdate.getE_total() / sitesCount,
                            mcUpdate.getAcceptanceRate(),
                            mcUpdate.getPhase() == MCSUpdate.Phase.EQUILIBRATION,
                            null));  // CFs no longer charted
                }
            });
        }

        MCResult result = builder.build().run();

        return new EquilibriumState(
                result.getTemperature(),
                result.getComposition(),
                result.getHmixPerSite(),   // enthalpy (energy of mixing per site)
                result.getHmixPerSite()    // free energy — G not directly from MCS; use Hmix as proxy
        );
    }

    /** Extracts temperature-dependent ECI: J(T) = a + b*T */
    private static double[] extractEci(CECTerm[] terms, double temperature) {
        double[] eci = new double[terms.length];
        for (int i = 0; i < terms.length; i++) {
            eci[i] = terms[i].a + terms[i].b * temperature;
        }
        return eci;
    }
}
