package org.ce.model.mcs;

import java.util.Arrays;

/**
 * Immutable raw simulation output from MCEngine.
 *
 * <p>Contains only what the simulator produces: mean quantities from averaging sweeps,
 * acceptance rate, and raw time series arrays for calculation-layer post-processing
 * (statistical analysis, error bars, derived quantities).</p>
 *
 * <p>No statistics fields (τ_int, SEM, jackknife Cv). Those are computed by
 * {@link org.ce.calculation.workflow.thermo.MCSStatisticsProcessor}.</p>
 */
public class MCResult {

    // ===== Simulation Parameters =====
    private final double   temperature;
    private final double[] composition;
    private final long     nEquilSweeps;
    private final long     nAvgSweeps;
    private final int      supercellSize;
    private final int      nSites;

    // ===== Mean Quantities (from averaging sweeps) =====
    private final double   acceptRate;
    private final double   energyPerSite;      // ⟨E⟩/N
    private final double   hmixPerSite;         // ⟨H⟩/N from ECI
    private final double[] avgCFs;              // mean CVCF correlation functions

    // ===== Raw Time Series (for calculation-layer statistics) =====
    private final double[] seriesHmix;          // Hmix per sweep (J, not J/site)
    private final double[] seriesE;             // total energy per sweep (J)
    private final double[][] seriesCF;          // CVCF per sweep [numCF][nAvg]

    /**
     * Constructs MCResult from simulation parameters and raw time series.
     *
     * @param temperature in K
     * @param composition mole fractions
     * @param nEquilSweeps equilibration sweep count
     * @param nAvgSweeps averaging sweep count
     * @param supercellSize L (for BCC: N = 2L³)
     * @param nSites total sites in supercell (N)
     * @param acceptRate acceptance rate during averaging phase
     * @param energyPerSite ⟨E⟩/N (J/site)
     * @param hmixPerSite ⟨H⟩/N from ECI (J/site)
     * @param avgCFs mean CVCF correlation functions (length = numNonPointCfs)
     * @param seriesHmix raw Hmix time series from averaging sweeps (length = nAvgSweeps)
     * @param seriesE raw total energy time series from averaging sweeps (length = nAvgSweeps)
     * @param seriesCF raw CVCF per sweep (dimensions: [numNonPointCfs][nAvgSweeps])
     */
    public MCResult(
            double temperature,
            double[] composition,
            long nEquilSweeps,
            long nAvgSweeps,
            int supercellSize,
            int nSites,
            double acceptRate,
            double energyPerSite,
            double hmixPerSite,
            double[] avgCFs,
            double[] seriesHmix,
            double[] seriesE,
            double[][] seriesCF) {

        this.temperature    = temperature;
        this.composition    = composition.clone();
        this.nEquilSweeps   = nEquilSweeps;
        this.nAvgSweeps     = nAvgSweeps;
        this.supercellSize  = supercellSize;
        this.nSites         = nSites;
        this.acceptRate     = acceptRate;
        this.energyPerSite  = energyPerSite;
        this.hmixPerSite    = hmixPerSite;
        this.avgCFs         = avgCFs != null ? avgCFs.clone() : null;
        this.seriesHmix     = seriesHmix != null ? seriesHmix.clone() : null;
        this.seriesE        = seriesE != null ? seriesE.clone() : null;
        this.seriesCF       = seriesCF != null ? deepClone(seriesCF) : null;
    }

    /** Helper to deep-clone 2D array. */
    private static double[][] deepClone(double[][] arr) {
        if (arr == null) return null;
        double[][] copy = new double[arr.length][];
        for (int i = 0; i < arr.length; i++) {
            copy[i] = arr[i] != null ? arr[i].clone() : null;
        }
        return copy;
    }

    // ===== Getters for Simulation Parameters =====

    public double getTemperature() {
        return temperature;
    }

    public double[] getComposition() {
        return composition.clone();
    }

    public long getNEquilSweeps() {
        return nEquilSweeps;
    }

    public long getNAvgSweeps() {
        return nAvgSweeps;
    }

    public int getSupercellSize() {
        return supercellSize;
    }

    public int getNSites() {
        return nSites;
    }

    // ===== Getters for Mean Quantities =====

    public double getAcceptRate() {
        return acceptRate;
    }

    public double getEnergyPerSite() {
        return energyPerSite;
    }

    public double getHmixPerSite() {
        return hmixPerSite;
    }

    public double[] getAvgCFs() {
        return avgCFs != null ? avgCFs.clone() : null;
    }

    // ===== Getters for Raw Time Series (for post-processing) =====

    /**
     * Returns the raw Hmix time series from averaging sweeps.
     * Each element is total Hmix for one sweep (in J, not J/site).
     * Length = nAvgSweeps.
     */
    public double[] getSeriesHmix() {
        return seriesHmix != null ? seriesHmix.clone() : null;
    }

    /**
     * Returns the raw total energy time series from averaging sweeps.
     * Each element is total E for one sweep (in J).
     * Length = nAvgSweeps.
     */
    public double[] getSeriesE() {
        return seriesE != null ? seriesE.clone() : null;
    }

    /**
     * Returns the raw CVCF per sweep from averaging phase.
     * Dimensions: [numNonPointCfs][nAvgSweeps].
     * Each seriesCF[i] is the time series for CVCF variable i.
     */
    public double[][] getSeriesCF() {
        return seriesCF != null ? deepClone(seriesCF) : null;
    }

    @Override
    public String toString() {
        return "MCResult{T=" + temperature + ", x=" + Arrays.toString(composition)
             + ", <E>/site=" + String.format("%.4f", energyPerSite)
             + ", acceptRate=" + String.format("%.3f", acceptRate)
             + ", nAvg=" + nAvgSweeps + "}";
    }

    /** Diagnostic printing utilities for internal verification. */
    public static final class Debug {

        private Debug() {}

        public static void printMcsSummary(MCResult result) {
            System.out.println("============================================================");
            System.out.println("  MCS SIMULATION OUTPUT");
            System.out.println("============================================================");
            System.out.println(String.format("  Temperature   : %.2f K", result.getTemperature()));
            System.out.println(String.format("  Composition   : %s", Arrays.toString(result.getComposition())));
            System.out.println(String.format("  Acceptance    : %.2f%%", result.getAcceptRate() * 100));
            System.out.println(String.format("  Equilib       : %d sweeps", result.getNEquilSweeps()));
            System.out.println(String.format("  Averaging     : %d sweeps", result.getNAvgSweeps()));
            System.out.println("------------------------------------------------------------");
            System.out.println(String.format("  <E>/site      : %.6f J/mol (no error bar in raw output)",
                    result.getEnergyPerSite()));
            System.out.println(String.format("  <H>/site      : %.6f J/mol (no error bar in raw output)",
                    result.getHmixPerSite()));
            System.out.println("  (Correlation functions and error bars computed by");
            System.out.println("   MCSStatisticsProcessor in calculation layer)");
            System.out.println("------------------------------------------------------------");

            double[] cfs = result.getAvgCFs();
            if (cfs != null && cfs.length > 0) {
                System.out.println("  Mean Correlation Functions (CVCF basis):");
                for (int i = 0; i < cfs.length; i++) {
                    System.out.println(String.format("    CF[%2d] = %+.6f", i, cfs[i]));
                }
            } else {
                System.out.println("  No CVCF correlation functions measured.");
            }

            double[] seriesH = result.getSeriesHmix();
            if (seriesH != null) {
                System.out.println(String.format("  Raw time series available: Hmix (%d points), E (%d points)",
                        seriesH.length, result.getSeriesE() != null ? result.getSeriesE().length : 0));
            }

            System.out.println("============================================================");
        }

        public static void printEciInfo(double[] eciCvcf, double[] eciOrth,
                                        org.ce.model.cvm.CvCfBasis basis,
                                        int N, int[] orbitSizes) {
            System.out.println("============================================================");
            System.out.println("  MCS ECI & ORBIT DIAGNOSTIC");
            System.out.println("============================================================");
            if (basis == null) {
                System.out.println("  [ERROR] CvCfBasis is null; ECI names unknown.");
            } else {
                int ncf = basis.numNonPointCfs;
                System.out.println(String.format("  CVCF Basis: %s (K=%d), numNonPointCfs=%d",
                        basis.structurePhase, basis.numComponents, ncf));
                System.out.println("  Non-zero CVCF ECIs:");
                boolean found = false;
                for (int l = 0; l < ncf; l++) {
                    if (Math.abs(eciCvcf[l]) > 1e-12) {
                        System.out.println(String.format("    idx %2d: %-8s (= %-8s) : %+.6f J/mol",
                                l, basis.eciNames.get(l), basis.cfNames.get(l), eciCvcf[l]));
                        found = true;
                    }
                }
                if (!found) System.out.println("    [WARNING] All CVCF ECIs are ZERO!");
            }

            if (orbitSizes != null && N > 0) {
                System.out.println("  Orbit Multiplicities (from ClusterData):");
                for (int t = 0; t < orbitSizes.length; t++) {
                    double mult = (double) orbitSizes[t] / N;
                    System.out.println(String.format("    Orbit %2d: counts=%-6d mult=%6.2f", t, orbitSizes[t], mult));
                }
            }

            if (eciOrth != null) {
                System.out.println("  Transformed ECIs (Metropolis ECIs, orthogonal basis):");
                for (int t = 0; t < eciOrth.length; t++) {
                    if (Math.abs(eciOrth[t]) > 1e-12) {
                        System.out.println(String.format("    Orbit %2d : %+.6f", t, eciOrth[t]));
                    }
                }
            }
            System.out.println("============================================================");
        }
    }
}
