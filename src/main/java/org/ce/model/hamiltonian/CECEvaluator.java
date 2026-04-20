package org.ce.model.hamiltonian;

import org.ce.model.cvm.CvCfBasis;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Unified evaluator for Cluster Expansion Coefficients (CEC).
 * Maps named Hamiltonian terms to basis correlation functions.
 */
public class CECEvaluator {

    private static final Logger LOG = Logger.getLogger(CECEvaluator.class.getName());

    public static double[] evaluate(CECEntry cec, double temperature, CvCfBasis basis, String label) {
        return evaluate(cec, temperature, basis, label, null);
    }

    /**
     * Evaluates ECI at a given temperature and maps them to the basis.
     * 
     * @param cec         Hamiltonian data
     * @param temperature system temperature in Kelvin
     * @param basis       target basis (defining name mapping and ordering)
     * @param label       label for logging (e.g., "CVM", "MCS")
     * @param sink        optional sink for user-visible output
     * @return array of evaluated ECIs indexed by basis.cfNames (0..numNonPointCfs-1)
     */
    public static double[] evaluate(CECEntry cec, double temperature, CvCfBasis basis, String label, java.util.function.Consumer<String> sink) {
        int ncf = basis.numNonPointCfs;
        double[] eci = new double[ncf];
        boolean[] mapped = new boolean[ncf];

        if (cec == null || cec.cecTerms == null || cec.cecTerms.length == 0) {
            LOG.warning(String.format("[%s] No CEC terms provided for %s—all interactions will be zero!", 
                    label, basis.structurePhase));
            return eci;
        }

        LOG.info(String.format("[%s-MAPPING] Loading CEC: elements=%s, notes=%s, ref=%s", 
                label, cec.elements, cec.notes, cec.reference));

        // 1. Build an optimized lookup map for this specific basis session
        Map<String, Integer> lookup = buildLookupMap(basis);

        // 2. Evaluate and map CEC terms into the ECI vector
        int loadCount = 0;
        java.util.List<String> unmatchedCecs = new java.util.ArrayList<>();
        
        for (CECEntry.CECTerm term : cec.cecTerms) {
            String termName = term.name;
            if (termName == null) continue;

            Integer idx = resolveIndex(termName, lookup, basis.totalCfs());

            if (idx != null && idx >= 0 && idx < ncf) {
                double val = term.a + term.b * temperature;
                eci[idx] = val;
                mapped[idx] = true;
                loadCount++;
                
                LOG.info(String.format("[%s-EXTRACT] %-8s (idx %d): file_a = %.2f -> eci[%d] = %.6f", 
                        label, termName, idx, term.a, idx, val));
            } else {
                unmatchedCecs.add(termName);
            }
        }

        // 3. Quality Report
        reportMappingQuality(label, basis, mapped, unmatchedCecs, loadCount, ncf);

        // 4. User-visible ECI report
        if (sink != null) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("  [%s Model] Evaluated ECIs at T=%.1f K:\n", label, temperature));
            for (int i = 0; i < ncf; i++) {
                if (mapped[i]) {
                    sb.append(String.format("    - %-12s: %10.6f J/mol\n", basis.cfNames.get(i), eci[i]));
                }
            }
            sink.accept(sb.toString().trim());
        }

        return eci;
    }

    /** Builds a case-insensitive map containing canonical names and e/v aliases. */
    private static Map<String, Integer> buildLookupMap(CvCfBasis basis) {
        Map<String, Integer> map = new HashMap<>(basis.numNonPointCfs * 3);
        for (int i = 0; i < basis.numNonPointCfs; i++) {
            String name = basis.cfNames.get(i).toLowerCase();
            map.put(name, i);
            
            // Normalize "v" basis functions to also accept "e" Hamiltonian terms
            if (name.startsWith("v")) {
                map.put("e" + name.substring(1), i);
            } else if (name.startsWith("e")) {
                map.put("v" + name.substring(1), i);
            }
        }
        return map;
    }

    /** Resolves a CEC term name to a basis index using the map or legacy fallbacks. */
    private static Integer resolveIndex(String name, Map<String, Integer> lookup, int maxTcf) {
        String key = name.toLowerCase();
        
        // 1. Map-based lookup (Handles Case-insensitive, Direct, and E/V aliases)
        Integer idx = lookup.get(key);
        if (idx != null) return idx;

        // 2. CF_i fallback (Legacy orthogonal nomenclature)
        if (key.startsWith("cf_")) {
            try {
                int col = Integer.parseInt(key.substring(3));
                if (col < maxTcf) return col;
            } catch (NumberFormatException ignored) {}
        }

        return null;
    }

    private static void reportMappingQuality(String label, CvCfBasis basis, boolean[] mapped, 
                                           java.util.List<String> unmatched, int loadCount, int ncf) {
        java.util.List<String> unmappedCfs = new java.util.ArrayList<>();
        for (int i = 0; i < ncf; i++) {
            if (!mapped[i]) unmappedCfs.add(basis.cfNames.get(i));
        }

        if (!unmappedCfs.isEmpty()) {
            LOG.warning(String.format("[%s-MAPPING] Unmapped CFs (missing from %s Hamiltonian): %s", 
                    label, basis.structurePhase, unmappedCfs));
        }
        if (!unmatched.isEmpty()) {
            LOG.warning(String.format("[%s-MAPPING] Unmatched CEC terms (not found in basis): %s", 
                    label, unmatched));
        }
        
        LOG.info(String.format("[%s-SUMMARY] ✓ Extracted %d ECIs (expected %d non-point CFs)", 
                label, loadCount, ncf));
    }
}
