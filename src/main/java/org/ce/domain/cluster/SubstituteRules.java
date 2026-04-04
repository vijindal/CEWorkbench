package org.ce.domain.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Mapping from site-operator products to CF indices.
 * 
 * <p>Includes geometry-augmented rules so that site-operator products from 
 * any equivalent set of sites (same inter-site distances, same basis indices) 
 * map to the correct CF.</p>
 */
public final class SubstituteRules {

    private final Map<SiteOpProductKey, CFIndex> rules;

    private SubstituteRules(Map<SiteOpProductKey, CFIndex> rules) {
        this.rules = Collections.unmodifiableMap(new LinkedHashMap<>(rules));
    }

    /**
     * Builds substitute rules with geometry augmentation.
     *
     * @param cfSiteOpList per-CF site-operator products from CF identification
     * @param siteList     site positions (from {@link SiteListBuilder}); must
     *                     not be {@code null}
     * @return substitute rules mapping site-op products to CF indices
     */
    public static SubstituteRules build(
            List<List<List<List<SiteOp>>>> cfSiteOpList,
            List<Position> siteList) {

        if (cfSiteOpList == null) {
            throw new IllegalArgumentException("cfSiteOpList must not be null");
        }
        if (siteList == null || siteList.isEmpty()) {
            throw new IllegalArgumentException("siteList must not be null or empty");
        }

        Map<SiteOpProductKey, CFIndex> rules = new LinkedHashMap<>();

        // --------------------------------------------------------------------
        // Step 1 – representative rules from CF identification
        // --------------------------------------------------------------------
        int maxBody = 0;
        int maxBasis = 0;

        for (int t = 0; t < cfSiteOpList.size(); t++) {
            List<List<List<SiteOp>>> typeGroups = cfSiteOpList.get(t);
            for (int j = 0; j < typeGroups.size(); j++) {
                List<List<SiteOp>> group = typeGroups.get(j);
                for (int k = 0; k < group.size(); k++) {
                    List<SiteOp> ops = group.get(k);
                    SiteOpProductKey key = new SiteOpProductKey(ops);
                    CFIndex index = new CFIndex(t, j, k);
                    CFIndex existing = rules.putIfAbsent(key, index);
                    if (existing != null && !existing.equals(index)) {
                        throw new IllegalStateException(
                                "Conflicting CF mapping for key " + key + ": "
                                        + existing + " vs " + index);
                    }
                    maxBody = Math.max(maxBody, ops.size());
                    for (SiteOp op : ops) {
                        maxBasis = Math.max(maxBasis, op.getBasisIndex());
                    }
                }
            }
        }

        // --------------------------------------------------------------------
        // Step 2 – augment with all geometry-equivalent products
        // --------------------------------------------------------------------
        if (maxBody > 0) {
            augmentWithGeometry(rules, siteList, maxBody, maxBasis);
        }

        return new SubstituteRules(rules);
    }

    public CFIndex lookup(List<SiteOp> ops) {
        return rules.get(new SiteOpProductKey(ops));
    }

    public Map<SiteOpProductKey, CFIndex> getRules() {
        return rules;
    }

    // =====================================================================
    // Geometry augmentation
    // =====================================================================

    /**
     * For every possible combination of sites and basis assignments,
     * computes the geometry signature and adds a mapping if it matches
     * a known CF representative.
     */
    private static void augmentWithGeometry(
            Map<SiteOpProductKey, CFIndex> rules,
            List<Position> siteList,
            int maxBody,
            int maxBasis) {

        // Build geometry map from existing representative rules
        Map<String, CFIndex> geoMap = new LinkedHashMap<>();
        for (Map.Entry<SiteOpProductKey, CFIndex> entry : rules.entrySet()) {
            String sig = geometrySignature(entry.getKey().getOps(), siteList);
            geoMap.putIfAbsent(sig, entry.getValue());
        }

        int nSites = siteList.size();

        // Generate all n-site combinations (n = 1 .. maxBody) with every
        // basis-index assignment, and map those whose geometry matches a CF.
        for (int n = 1; n <= Math.min(maxBody, nSites); n++) {
            final int nn = n;
            forEachCombination(nSites, n, combo -> {
                forEachBasisTuple(nn, maxBasis, bases -> {
                    List<SiteOp> ops = new ArrayList<>(nn);
                    for (int i = 0; i < nn; i++) {
                        ops.add(new SiteOp(combo[i], bases[i]));
                    }
                    SiteOpProductKey key = new SiteOpProductKey(ops);
                    if (!rules.containsKey(key)) {
                        String sig = geometrySignature(ops, siteList);
                        CFIndex cfIdx = geoMap.get(sig);
                        if (cfIdx != null) {
                            rules.put(key, cfIdx);
                        }
                    }
                });
            });
        }
    }

    /**
     * Returns a canonical string encoding the geometry and basis pattern
     * of a site-operator product.  Two products are equivalent iff they
     * have the same sorted basis indices <em>and</em> the same sorted
     * pairwise inter-site distances.
     */
    static String geometrySignature(List<SiteOp> ops, List<Position> siteList) {
        int n = ops.size();

        // Sorted basis indices
        List<Integer> bases = new ArrayList<>(n);
        for (SiteOp op : ops) {
            bases.add(op.getBasisIndex());
        }
        bases.sort(Integer::compareTo);

        // Sorted pairwise distances (rounded to microsecond-level to avoid
        // floating-point noise)
        List<Long> dists = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            Position pi = siteList.get(ops.get(i).getSiteIndex());
            for (int j = i + 1; j < n; j++) {
                Position pj = siteList.get(ops.get(j).getSiteIndex());
                dists.add(Math.round(pi.distance(pj) * 1_000_000));
            }
        }
        dists.sort(Long::compareTo);

        return bases + "|" + dists;
    }

    // =====================================================================
    // Combination generators
    // =====================================================================

    /** Calls {@code action} for each {@code k}-combination of indices 0..n-1. */
    private static void forEachCombination(int n, int k, Consumer<int[]> action) {
        int[] combo = new int[k];
        generateCombos(combo, 0, 0, n, k, action);
    }

    private static void generateCombos(
            int[] combo, int start, int depth, int n, int k,
            Consumer<int[]> action) {
        if (depth == k) {
            action.accept(combo);
            return;
        }
        for (int i = start; i < n; i++) {
            combo[depth] = i;
            generateCombos(combo, i + 1, depth + 1, n, k, action);
        }
    }

    /** Calls {@code action} for each {@code k}-tuple of basis indices 1..maxBasis. */
    private static void forEachBasisTuple(int k, int maxBasis, Consumer<int[]> action) {
        int[] tuple = new int[k];
        generateBasisTuples(tuple, 0, k, maxBasis, action);
    }

    private static void generateBasisTuples(
            int[] tuple, int depth, int k, int maxBasis,
            Consumer<int[]> action) {
        if (depth == k) {
            action.accept(tuple);
            return;
        }
        for (int b = 1; b <= maxBasis; b++) {
            tuple[depth] = b;
            generateBasisTuples(tuple, depth + 1, k, maxBasis, action);
        }
    }
}
