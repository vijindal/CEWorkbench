package org.ce.model.hamiltonian;

import org.ce.model.cluster.cvcf.CvCfBasis;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Unified evaluator for Cluster Expansion Coefficients (CEC).
 * Maps named Hamiltonian terms to basis correlation functions.
 */
public class CECEvaluator {

    private static final Logger LOG = Logger.getLogger(CECEvaluator.class.getName());

    /**
     * Evaluates ECI at a given temperature and maps them to the basis.
     * 
     * @param cec         Hamiltonian data
     * @param temperature system temperature in Kelvin
     * @param basis       target basis (defining name mapping and ordering)
     * @param label       label for logging (e.g., "CVM", "MCS")
     * @return array of evaluated ECIs indexed by basis.cfNames (0..numNonPointCfs-1)
     */
    public static double[] evaluate(CECEntry cec, double temperature, CvCfBasis basis, String label) {
        int ncf = basis.numNonPointCfs;
        double[] eci = new double[ncf];
        boolean[] mapped = new boolean[ncf];
        List<String> unmatchedCecs = new ArrayList<>();

        if (cec == null || cec.cecTerms == null || cec.cecTerms.length == 0) {
            LOG.warning(String.format("[%s] No CEC terms provided for %s—all interactions will be zero!", 
                    label, basis.structurePhase));
            return eci;
        }

        LOG.info(String.format("[%s-MAPPING] Loading CEC: elements=%s, notes=%s, ref=%s", 
                label, cec.elements, cec.notes, cec.reference));

        int loadCount = 0;
        for (CECEntry.CECTerm term : cec.cecTerms) {
            String termName = term.name;
            if (termName == null) continue;

            int idx = findCfIndex(termName, basis);

            if (idx >= 0 && idx < ncf) {
                double val = term.a + term.b * temperature;
                eci[idx] = val;
                mapped[idx] = true;
                loadCount++;
                
                // Detailed extraction log (matches former MCSEngine debug style)
                LOG.info(String.format("[%s-EXTRACT] %-8s (idx %d): file_a = %.2f -> eci[%d] = %.6f", 
                        label, basis.cfNames.get(idx), idx, term.a, idx, val));
            } else {
                unmatchedCecs.add(termName);
            }
        }

        // --- Quality Report ---
        List<String> unmappedCfs = new ArrayList<>();
        for (int i = 0; i < ncf; i++) {
            if (!mapped[i]) unmappedCfs.add(basis.cfNames.get(i));
        }

        if (!unmappedCfs.isEmpty()) {
            LOG.warning(String.format("[%s-MAPPING] Unmapped CFs (missing from %s Hamiltonian): %s", 
                    label, cec.elements, unmappedCfs));
        }
        if (!unmatchedCecs.isEmpty()) {
            LOG.warning(String.format("[%s-MAPPING] Unmatched CEC terms (not found in %s basis): %s", 
                    label, basis.structurePhase, unmatchedCecs));
        }
        
        LOG.info(String.format("[%s-SUMMARY] ✓ Extracted %d ECIs (expected %d non-point CFs)", 
                label, loadCount, ncf));

        return eci;
    }

    /**
     * Flexible CF name matcher based on CVMEngine and MCSEngine legacy logic.
     */
    private static int findCfIndex(String cecName, CvCfBasis basis) {
        if (cecName == null) return -1;

        // 1. Direct match (e.g., e4AB -> e4AB)
        int idx = basis.indexOfCf(cecName);
        if (idx >= 0) return idx;

        // 2. e -> v translation (common in CVCF basis: v4AB vs e4AB)
        if (cecName.startsWith("e")) {
            String vName = "v" + cecName.substring(1);
            idx = basis.indexOfCf(vName);
            if (idx >= 0) return idx;
        }

        // 3. Fallback: Case-insensitive search 
        // Handles legacy nomenclature and cross-platform naming differences
        for (int i = 0; i < basis.numNonPointCfs; i++) {
            String bName = basis.cfNames.get(i);
            if (bName.equalsIgnoreCase(cecName)) return i;
            if (cecName.startsWith("e") && bName.equalsIgnoreCase("v" + cecName.substring(1))) return i;
            if (cecName.startsWith("v") && bName.equalsIgnoreCase("e" + cecName.substring(1))) return i;
        }

        // 4. CF_i fallback (Legacy orthogonal nomenclature)
        if (cecName.startsWith("CF_")) {
            try {
                int col = Integer.parseInt(cecName.substring(3));
                if (col < basis.totalCfs()) return col;
            } catch (Exception ignored) {}
        }

        return -1;
    }
}
