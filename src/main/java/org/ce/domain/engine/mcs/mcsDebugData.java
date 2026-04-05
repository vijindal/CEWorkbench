package org.ce.domain.engine.mcs;

import org.ce.domain.cluster.cvcf.CvCfBasis;
import java.util.Arrays;

/**
 * Temporary debug class to print MCS calculation results for internal verification.
 */
public final class mcsDebugData {

    private mcsDebugData() {}

    public static void printMcsSummary(MCResult result) {
        System.out.println("============================================================");
        System.out.println("  MCS CALCULATION DEBUG SUMMARY");
        System.out.println("============================================================");
        System.out.println(String.format("  Temperature   : %.2f K", result.getTemperature()));
        System.out.println(String.format("  Composition   : %s", Arrays.toString(result.getComposition())));
        System.out.println(String.format("  Acceptance    : %.2f%%", result.getAcceptRate() * 100));
        System.out.println(String.format("  Avg Sweeps    : %d", result.getNAvgSweeps()));
        System.out.println(String.format("  Tau_int       : %.2f sweeps", result.getTauInt()));
        System.out.println(String.format("  N_eff         : %d", result.getNEff()));
        System.out.println("------------------------------------------------------------");
        System.out.println(String.format("  <E>/site      : %.6f ± %.6f J/mol", 
                result.getEnergyPerSite(), result.getStdEnergyPerSite()));
        System.out.println(String.format("  <H>/site      : %.6f ± %.6f J/mol", 
                result.getHmixPerSite(), result.getStdHmixPerSite()));
        System.out.println(String.format("  Cv (Jackknife): %.6f ± %.6f J/(mol·K)", 
                result.getCvJackknife(), result.getCvStdErr()));
        System.out.println("------------------------------------------------------------");
        
        double[] cfs = result.getAvgCFs();
        double[] std = result.getStdCFs();
        if (cfs != null && cfs.length > 0) {
            System.out.println("  Mean Correlation Functions (CVCF basis):");
            for (int i = 0; i < cfs.length; i++) {
                double sigma = (std != null && i < std.length) ? std[i] : 0.0;
                System.out.println(String.format("    CF[%2d] = %+.6f ± %.6f", i, cfs[i], sigma));
            }
        } else {
            System.out.println("  No CVCF correlation functions measured.");
        }
        System.out.println("============================================================");
    }

    public static void printEciInfo(double[] eciCvcf, double[] eciOrth, CvCfBasis basis, 
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
