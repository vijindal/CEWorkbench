package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

/**
 * Diagnostic: reports per-CF-column embedding counts, maxEmbPerCol, DeltaScratch
 * array sizes, and average per-site affected-embedding count for a given system/L,
 * to find why deltaEExchangeCvcf costs scale disproportionately with system size
 * for larger component counts.
 */
public class EmbeddingScaleProbe {

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        String elements = args.length > 0 ? args[0] : "Nb-Ti-V-Zr";
        int L = args.length > 1 ? Integer.parseInt(args[1]) : 6;

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);

        System.out.println("=== Embedding Scale Probe (" + elements + " BCC_A2, L=" + L + ") ===");

        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis != null ? geo.basis.numNonPointCfs : geo.ncf;

        System.out.println("N sites = " + N + ", ncf = " + ncf);
        System.out.println();
        System.out.println("Per-CF-column embedding counts:");
        int maxEmbPerCol = 0;
        long totalEmb = 0;
        for (int l = 0; l < geo.cfEmbeddings.size(); l++) {
            int size = geo.cfEmbeddings.get(l) != null ? geo.cfEmbeddings.get(l).size() : 0;
            maxEmbPerCol = Math.max(maxEmbPerCol, size);
            totalEmb += size;
            System.out.printf("  cf[%2d]: %,d embeddings%n", l, size);
        }
        System.out.println();
        System.out.println("maxEmbPerCol = " + maxEmbPerCol);
        System.out.println("totalEmbeddings (sum over columns) = " + totalEmb);
        System.out.printf("DeltaScratch.seen size = ncf * maxEmbPerCol = %d * %d = %,d booleans (%.2f MB)%n",
                ncf, maxEmbPerCol, (long) ncf * maxEmbPerCol, (ncf * (long) maxEmbPerCol) / 1e6);

        // Per-site affected-embedding count via siteToCfIndex CSR structure.
        Embeddings.CsrSiteToCfIndex idx = geo.siteToCfIndex;
        long totalSiteEntries = 0;
        int maxPerSite = 0;
        for (int s = 0; s < N; s++) {
            int count = idx.offsets[s + 1] - idx.offsets[s];
            totalSiteEntries += count;
            maxPerSite = Math.max(maxPerSite, count);
        }
        System.out.println();
        System.out.printf("Per-site (cfCol,embIdx) entries: total=%,d, avg=%.1f, max=%d%n",
                totalSiteEntries, (double) totalSiteEntries / N, maxPerSite);
        System.out.println("(A single exchange touches up to 2x this many entries before dedup —"
                + " this bounds the per-call affected-embedding count.)");
    }
}
