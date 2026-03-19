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
    private static final int    DEFAULT_L          = 4;
    private static final int    DEFAULT_EQUIL      = 1000;
    private static final int    DEFAULT_AVG        = 2000;

    @Override
    public EquilibriumState compute(ThermodynamicInput input) throws Exception {

        ClusCoordListResult clusterData =
                input.clusterData.getDisorderedClusterResult().getDisClusterData();

        double[] eci = extractEci(input.cec.cecTerms, input.temperature);

        int totalSweeps = DEFAULT_EQUIL + DEFAULT_AVG;
        int N           = 2 * DEFAULT_L * DEFAULT_L * DEFAULT_L;  // BCC sites

        MCSRunner.Builder builder = MCSRunner.builder()
                .clusterData(clusterData)
                .eci(eci)
                .numComp(input.composition.length)
                .T(input.temperature)
                .composition(input.composition)
                .nEquil(DEFAULT_EQUIL)
                .nAvg(DEFAULT_AVG)
                .L(DEFAULT_L)
                .seed(System.currentTimeMillis())
                .R(GAS_CONSTANT);

        if (input.progressSink != null) {
            Consumer<String> sink = input.progressSink;
            sink.accept(String.format("  MCS started: L=%d (%d sites), %d equil + %d avg sweeps",
                    DEFAULT_L, N, DEFAULT_EQUIL, DEFAULT_AVG));
            builder.updateListener(update -> {
                if (update.getStep() % 100 == 0 || update.getStep() == totalSweeps) {
                    sink.accept(String.format("  [%-5s] sweep %4d/%-4d  <E>/site=%9.5f  accept=%5.1f%%",
                            update.getPhase(),
                            update.getStep(), totalSweeps,
                            update.getE_total() / N,
                            update.getAcceptanceRate() * 100));
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
