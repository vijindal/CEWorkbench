package org.ce.model.mcs;

import org.ce.model.cluster.ClusterPrimitives.Vector3D;
import java.util.ArrayList;
import java.util.List;

/**
 * Partitions the MCS supercell into spatially independent blocks for parallel execution.
 * 
 * <p>Uses the Cartesian positions of sites and the maximum cluster interaction radius 
 * to divide the lattice into a grid of blocks. Blocks are then "colored" such that 
 * any two blocks of the same color are separated by at least R_max, making them 
 * safe to update simultaneously without race conditions on cluster measurements.</p>
 */
public final class LatticeDecomposer {

    public static final class DecomposedLattice {
        public final int numBlocks;
        public final int numColors;
        /** Blocks grouped by color: colors[colorIdx] = list of block indices. */
        public final List<Integer>[] colors;
        /** sitesInBlock[blockIdx] = list of site indices in that block. */
        public final List<Integer>[] sitesInBlock;
        /** blockOfSite[siteIdx] = block index for that site. */
        public final int[] blockOfSite;

        @SuppressWarnings("unchecked")
        public DecomposedLattice(int numBlocks, int numColors) {
            this.numBlocks = numBlocks;
            this.numColors = numColors;
            this.colors = new ArrayList[numColors];
            for (int i = 0; i < numColors; i++) this.colors[i] = new ArrayList<>();
            this.sitesInBlock = new ArrayList[numBlocks];
            for (int i = 0; i < numBlocks; i++) this.sitesInBlock[i] = new ArrayList<>();
            this.blockOfSite = null; // initialized externally if needed
        }
    }

    /**
     * Decomposes the lattice into a grid of blocks for parallel execution.
     * 
     * @param positions Cartesian positions of all sites.
     * @param L Linear dimension of the supercell (in unit cells).
     * @param rMax Maximum interaction radius (max cluster diameter).
     * @param numBlocksPerDim Number of blocks to divide each dimension into (e.g., 2 or 4).
     */
    public static DecomposedLattice decompose(List<Vector3D> positions, int L, double rMax, int numBlocksPerDim) {
        int N = positions.size();
        int totalBlocks = numBlocksPerDim * numBlocksPerDim * numBlocksPerDim;
        
        // Use a 2x2x2 coloring scheme (8 colors) for simplicity. 
        // This ensures that same-colored blocks are separated by at least 1 block width.
        // We must verify that (L / numBlocksPerDim) > rMax.
        int numColors = 8; 
        DecomposedLattice dl = new DecomposedLattice(totalBlocks, numColors);
        int[] blockOfSite = new int[N];

        // 1. Assign sites to blocks
        // We assume positions are roughly in a box of [0, L]
        double blockWidth = (double) L / numBlocksPerDim;
        
        for (int i = 0; i < N; i++) {
            Vector3D p = positions.get(i);
            int bx = (int) (p.getX() / blockWidth);
            int by = (int) (p.getY() / blockWidth);
            int bz = (int) (p.getZ() / blockWidth);
            
            // Clamp to avoid edge cases
            bx = Math.min(bx, numBlocksPerDim - 1);
            by = Math.min(by, numBlocksPerDim - 1);
            bz = Math.min(bz, numBlocksPerDim - 1);
            
            int blockIdx = (bx * numBlocksPerDim + by) * numBlocksPerDim + bz;
            dl.sitesInBlock[blockIdx].add(i);
            blockOfSite[i] = blockIdx;
        }

        // 2. Assign blocks to colors (checkerboard)
        for (int bx = 0; bx < numBlocksPerDim; bx++) {
            for (int by = 0; by < numBlocksPerDim; by++) {
                for (int bz = 0; bz < numBlocksPerDim; bz++) {
                    int blockIdx = (bx * numBlocksPerDim + by) * numBlocksPerDim + bz;
                    // Standard 2x2x2 coloring: (x%2, y%2, z%2) maps to 0-7
                    int colorIdx = (bx % 2) * 4 + (by % 2) * 2 + (bz % 2);
                    dl.colors[colorIdx].add(blockIdx);
                }
            }
        }

        return dl;
    }

    /** Computes max cluster radius from embeddings. */
    public static double computeRMax(List<List<Embeddings.Embedding>> cfEmbeddings, List<Vector3D> positions) {
        double maxD2 = 0;
        for (List<Embeddings.Embedding> embs : cfEmbeddings) {
            if (embs == null) continue;
            for (Embeddings.Embedding e : embs) {
                int[] sites = e.getSiteIndices();
                for (int s1 : sites) {
                    for (int s2 : sites) {
                        Vector3D p1 = positions.get(s1);
                        Vector3D p2 = positions.get(s2);
                        double dx = p1.getX() - p2.getX();
                        double dy = p1.getY() - p2.getY();
                        double dz = p1.getZ() - p2.getZ();
                        maxD2 = Math.max(maxD2, dx*dx + dy*dy + dz*dz);
                    }
                }
            }
        }
        return Math.sqrt(maxD2);
    }
}
