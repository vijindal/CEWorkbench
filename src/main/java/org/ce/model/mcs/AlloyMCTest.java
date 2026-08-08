package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.ThermodynamicResult;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;
import java.util.Arrays;
import java.util.List;
import org.ce.model.mcs.Embeddings;

/**
 * Quick test to verify AlloyMC constructor and state management.
 * Specifically targets the Nb-Ti BCC_A2 system.
 */
public class AlloyMCTest {
    public static void main(String[] args) {
        try {
            runTest();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void runTest() throws Exception {
        System.out.println("=== AlloyMC Multicore Parallel Test (L=8) ===");

        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);
        Workspace.SystemId id = new Workspace.SystemId("Nb-Ti-V", "BCC_A2", "T");

        System.out.println("\n[1] Building Shared ModelSession & Geometry...");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        // Build geometry once (expensive)
        MCSGeometry geo = MCSGeometry.build(session, 8, null);

        int nChains = 4;
        int nEquil = 100;
        int nAvg = 500;
        double temp = 1000.0;
        double[] x = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

        System.out.println(
                String.format("\n[2] Launching %d Parallel Chains (Equil=%d, Avg=%d)...", nChains, nEquil, nAvg));

        long startTime = System.currentTimeMillis();

        List<ThermodynamicResult> results = java.util.stream.IntStream.range(0, nChains)
                .parallel()
                .mapToObj(i -> {
                    System.out.println("  > Starting Chain " + i + "...");
                    AlloyMC engine = new AlloyMC(session, 8, null);
                    engine.setSeed(System.nanoTime() + i * 1000000L);
                    engine.setTemperature(temp);
                    engine.setComposition(x);
                    engine.run(nEquil, nAvg);

                    return new ThermodynamicResult(
                            temp, x, Double.NaN,
                            engine.getAverageEnergyPerSite(),
                            Double.NaN,
                            engine.getStdDevEnergyPerSite(),
                            Double.NaN,
                            null,
                            engine.getAverageCvcf(),
                            engine.getStdDevCvcf());
                })
                .collect(java.util.stream.Collectors.toList());

        long duration = System.currentTimeMillis() - startTime;
        System.out.println(String.format("\n[3] All Chains Complete in %.2f seconds.", duration / 1000.0));

        // 4. Aggregate Results
        aggregateAndPrint(results);
    }

    private static void aggregateAndPrint(List<ThermodynamicResult> results) {
        int m = results.size();
        double sumH = 0;
        double sumVarH = 0;
        int ncf = results.get(0).avgCFs.length;
        double[] sumCF = new double[ncf];
        double[] sumVarCF = new double[ncf];

        for (ThermodynamicResult r : results) {
            sumH += r.enthalpy;
            sumVarH += Math.pow(r.stdEnthalpy, 2);
            for (int i = 0; i < ncf; i++) {
                sumCF[i] += r.avgCFs[i];
                sumVarCF[i] += Math.pow(r.stdCFs[i], 2);
            }
        }

        double finalH = sumH / m;
        double finalStdH = Math.sqrt(sumVarH / m); // Pooled within-chain std

        System.out.println("\n=== Aggregated Results (Multicore) ===");
        System.out.format("H (J/mol) : %.4f \u00B1 %.4f (pooled \u03C3)\n", finalH, finalStdH);

        System.out.println("\nAggregated CFs (Top 5):");
        for (int i = 0; i < Math.min(5, ncf); i++) {
            double avg = sumCF[i] / m;
            double std = Math.sqrt(sumVarCF[i] / m);
            System.out.format("  [%02d] : %9.6f \u00B1 %7.6f\n", i, avg, std);
        }
    }

    private static void traceOrthoCFSteps(AlloyMC engine) {
        MCSGeometry geo = engine.getGeo();
        int ncf = engine.getNcf();
        int numComp = engine.getNumComp();
        double[] flatBasis = geo.getFlatBasisMatrix();
        int[] occ = engine.getConfig().getRawOcc();

        System.out.println("\n[TRACE] Calculation for all Orthogonal CFs (showing Emb 0 / Average):");

        int totalCols = ncf + numComp; // All CFs including points and constant
        int[][] allCfBasisIndices = geo.getCfBasisIndices();

        for (int l = 0; l < totalCols; l++) {
            int[] basisIndices = (l < allCfBasisIndices.length) ? allCfBasisIndices[l] : new int[0];
            StringBuilder def = new StringBuilder();
            for (int i = 0; i < basisIndices.length; i++) {
                if (i > 0)
                    def.append(" * ");
                def.append("phi" + basisIndices[i]);
            }
            if (basisIndices.length == 0)
                def.append("Constant");

            List<Embeddings.Embedding> embs = (l < geo.cfEmbeddings.size()) ? geo.cfEmbeddings.get(l) : null;

            if (embs != null && !embs.isEmpty()) {
                // Non-Point CF: Trace the FIRST embedding
                Embeddings.Embedding e = embs.get(0);
                int[] sites = e.getSiteIndices();
                int[] alphas = e.getAlphaIndices();
                double prod = 1.0;
                StringBuilder sb = new StringBuilder();
                for (int k = 0; k < sites.length; k++) {
                    double val = flatBasis[occ[sites[k]] * numComp + alphas[k]];
                    prod *= val;
                    sb.append(String.format("phi[%d](site %d)=%.1f ", alphas[k], sites[k], val));
                }
                System.out.println(
                        String.format("  Col %2d (nEmb=%4d): Def=%-18s | Emb 0: sites=%-12s, values=[%s] -> prod=%.1f",
                                l, embs.size(), def.toString(), Arrays.toString(sites), sb.toString().trim(), prod));
            } else {
                // Point or Constant CF: Trace as Site Average
                int k = (basisIndices.length > 0) ? basisIndices[0] : 0;
                double sum = 0.0;
                for (int s = 0; s < geo.nSites(); s++) {
                    sum += flatBasis[occ[s] * numComp + k];
                }
                double avg = sum / geo.nSites();

                if (basisIndices.length > 0) {
                    System.out.println(
                            String.format("  Col %2d (nEmb=%4d): Def=%-18s | Average of all %d sites -> result=%.1f",
                                    l, geo.nSites(), def.toString(), geo.nSites(), avg));
                } else {
                    System.out.println(String.format("  Col %2d (nEmb=%4d): Def=%-18s | Always constant -> result=%.1f",
                            l, 1, def.toString(), 1.0));
                }
            }
        }
    }

    private static void printSiteDecorations(AlloyMC engine) {
        int N = engine.getGeo().nSites();
        int[] config = engine.getConfig().getRawOcc();
        double[] basis = org.ce.model.cluster.ClusterMath.buildBasis(engine.getNumComp());

        System.out.println("Site Occupations & Basis Values (phi1=sigma, phi2=sigma^2) [Clipped to 10]:");
        for (int i = 0; i < Math.min(10, N); i++) {
            int occ = config[i];
            double s = basis[occ];
            System.out.println(String.format("  Site %3d: occ=%d (element=%s), sigma=%+.1f, sigma^2=%.1f",
                    i, occ, elementSymbol(occ), s, s * s));
        }
        if (N > 10)
            System.out.println("  ...");
    }

    private static String elementSymbol(int occ) {
        switch (occ) {
            case 0:
                return "Nb";
            case 1:
                return "Ti";
            case 2:
                return "V";
            default:
                return "?";
        }
    }

    private static void auditOrthoCFs(AlloyMC engine) {
        double[] uOrth = engine.getCorrelationFunctions();
        int ncf = engine.getNcf();
        int numComp = engine.getNumComp();

        System.out.println("\nMeasured Orthogonal CFs:");
        System.out.println("  Non-Point (0.." + (ncf - 1) + "): " + Arrays.toString(Arrays.copyOfRange(uOrth, 0, ncf)));
        System.out.println(
                "  Point CFs (phi1, phi2): " + Arrays.toString(Arrays.copyOfRange(uOrth, ncf, ncf + numComp - 1)));

        System.out.println("  Constant: " + uOrth[uOrth.length - 1]);

        System.out.println("\nFull CVCF Basis Vector (Recovered via Tinv):");
        double[] cvcf = engine.getCvcfCorrelationFunctions();
        List<String> cfNames = engine.getGeo().getBasis().cfNames;
        List<String> definitions = engine.getGeo().getBasis().cfDefinitions;

        for (int i = 0; i < cvcf.length; i++) {
            String name = (i < cfNames.size()) ? cfNames.get(i) : "unknown";
            String def = (i < definitions.size()) ? definitions.get(i) : "";
            System.out.println(String.format("  Col %2d [%-7s]: %-40s = %.4f", i, name, def, cvcf[i]));
        }
    }
}
