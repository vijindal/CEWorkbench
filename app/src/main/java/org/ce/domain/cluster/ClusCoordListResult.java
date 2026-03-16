package org.ce.domain.cluster;

import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.Site;
import java.util.List;

/**
 * Immutable result object returned by {@link org.ce.domain.identification.engine.ClusCoordListGenerator}.
 *
 * <p>Encapsulates the complete output of one run of the cluster-coordinate
 * list generation algorithm, which is the Java equivalent of the Mathematica
 * {@code genClusCoordList} function.  The fields mirror the Mathematica
 * return tuple:</p>
 *
 * <table border="1" cellpadding="4">
 *   <tr><th>Java field</th><th>Mathematica equivalent</th><th>Description</th></tr>
 *   <tr><td>{@code clusCoordList}</td><td>{@code clusCoordList}</td>
 *       <td>Canonical representative cluster for each type, sorted by
 *           descending size.</td></tr>
 *   <tr><td>{@code multiplicities}</td><td>{@code clusMList}</td>
 *       <td>Normalised orbit size for each cluster type
 *           (raw orbit size / number of point-cluster orbit members).</td></tr>
 *   <tr><td>{@code orbitList}</td><td>{@code subClusOrbitList}</td>
 *       <td>Full symmetry orbit (list of equivalent clusters) for each type.</td></tr>
 *   <tr><td>{@code rcList}</td><td>{@code rc}</td>
 *       <td>Per-sublattice site counts for each cluster type.</td></tr>
 *   <tr><td>{@code tc}</td><td>{@code tc}</td>
 *       <td>Total number of distinct cluster types found.</td></tr>
 *   <tr><td>{@code numPointSubClusFound}</td><td>(internal counter)</td>
 *       <td>Number of point sub-clusters found (used for normalisation).</td></tr>
 * </table>
 *
 * @author  CVM Project
 * @version 1.0
 * @see     org.ce.domain.identification.engine.ClusCoordListGenerator
 */
public class ClusCoordListResult {

    private final List<Cluster>          clusCoordList;
    private final List<Double>           multiplicities;
    private final List<List<Cluster>>    orbitList;
    private final List<List<Integer>>    rcList;
    private final int                    tc;
    private final int                    numPointSubClusFound;

    /**
     * Constructs a fully populated result.
     *
     * @param clusCoordList       canonical representatives; must not be {@code null}
     * @param multiplicities      normalised multiplicities (same length as clusCoordList)
     * @param orbitList           orbits (same length as clusCoordList)
     * @param rcList              per-sublattice site counts (same length as clusCoordList)
     * @param tc                  total distinct cluster-type count
     * @param numPointSubClusFound number of point sub-clusters used for normalisation
     */
    public ClusCoordListResult(List<Cluster>       clusCoordList,
                               List<Double>        multiplicities,
                               List<List<Cluster>> orbitList,
                               List<List<Integer>> rcList,
                               int                 tc,
                               int                 numPointSubClusFound) {
        this.clusCoordList        = clusCoordList;
        this.multiplicities       = multiplicities;
        this.orbitList            = orbitList;
        this.rcList               = rcList;
        this.tc                   = tc;
        this.numPointSubClusFound = numPointSubClusFound;
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    /**
     * Returns the canonical representative cluster for each cluster type,
     * sorted by descending number of sites.
     *
     * @return list of representative {@link Cluster} objects; never {@code null}
     */
    public List<Cluster> getClusCoordList() { return clusCoordList; }

    /**
     * Returns the normalised multiplicity for each cluster type.
     *
     * <p>Normalisation is performed by dividing the raw orbit size by the
     * total number of point-cluster orbit members, matching the Mathematica
     * convention.</p>
     *
     * @return list of normalised multiplicities; same length as {@link #getClusCoordList()}
     */
    public List<Double> getMultiplicities() { return multiplicities; }

    /**
     * Returns the full symmetry orbit for each cluster type.
     *
     * @return list of orbits; each orbit is a list of symmetry-equivalent
     *         {@link Cluster} objects; never {@code null}
     */
    public List<List<Cluster>> getOrbitList() { return orbitList; }

    /**
     * Returns the per-sublattice site counts for each cluster type.
     *
     * <p>For a disordered (single-sublattice) cluster, each inner list
     * has exactly one element equal to the total number of sites in the cluster.</p>
     *
     * @return list of site-count lists; never {@code null}
     */
    public List<List<Integer>> getRcList() { return rcList; }

    /**
     * Returns the total number of distinct cluster types that were found.
     *
     * @return {@code tc} â‰¥ 0
     */
    public int getTc() { return tc; }

    /**
     * Returns the number of point sub-clusters found during the generation,
     * used internally for multiplicity normalisation.
     *
     * @return count â‰¥ 0
     */
    public int getNumPointSubClusFound() { return numPointSubClusFound; }

    // -------------------------------------------------------------------------
    // Debug
    // -------------------------------------------------------------------------

    /**
     * Prints a structured debug summary of this result to standard output.
     *
     * <p>Output format:</p>
     * <pre>
     * [ClusCoordListResult]
     *   total cluster types (tc)     : 7
     *   point sub-clusters found     : 2
     *   type[0] : sites=4, mult=3.0000, rc=[4], orbit size=6
     *     sites: (0.0,0.0,0.0) (0.5,0.5,0.5) ...
     *   type[1] : ...
     * </pre>
     */
    public void printDebug() {
        System.out.println("[ClusCoordListResult]");
        System.out.println("  total cluster types (tc)  : " + tc);
        System.out.println("  point sub-clusters found  : " + numPointSubClusFound);
        System.out.println("  number of entries         : " + clusCoordList.size());
        for (int i = 0; i < clusCoordList.size(); i++) {
            List<Site> sites = clusCoordList.get(i).getAllSites();
            System.out.printf("  type[%d] : sites=%d, mult=%.4f, rc=%s, orbitSize=%d%n",
                    i,
                    sites.size(),
                    multiplicities.get(i),
                    rcList.get(i),
                    orbitList.get(i).size());
            StringBuilder sb = new StringBuilder("    positions:");
            for (Site s : sites) {
                sb.append(" ").append(s.getPosition());
                if (s.getSymbol() != null)
                    sb.append("(").append(s.getSymbol()).append(")");
            }
            System.out.println(sb);
        }
    }
}



