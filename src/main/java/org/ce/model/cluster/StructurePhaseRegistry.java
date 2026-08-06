package org.ce.model.cluster;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Single source of truth for the ordered-phase → disordered-parent (HSP)
 * structure mapping used by cluster/CF identification.
 *
 * <p>
 * A cluster file is named {@code <structure>-<model>.txt} and a symmetry
 * file {@code <structure>-SG.txt} (see {@code Workspace} javadoc). The
 * structure part of both filenames MUST match — the pipeline result is
 * physically meaningless otherwise. This registry lets the GUI/CLI derive
 * the disordered structure, and both symmetry files, from a single
 * ordered-phase selection instead of asking the user to pick four
 * independent files that must agree.
 * </p>
 *
 * <h2>How to add a new ordered/disordered structure pair</h2>
 * <p>
 * Add one line to {@link #PARENT_OF}: {@code "ORDERED_STRUCTURE", "DISORDERED_PARENT"}.
 * A disordered (parent) structure does not need an entry — any structure
 * absent from this map is treated as its own parent (self-referential),
 * which is correct for HSP structures like {@code BCC_A2}, {@code FCC_A1},
 * {@code HCP_A3}. Do not add self-mappings; they are the default.
 * </p>
 */
public final class StructurePhaseRegistry {

    private StructurePhaseRegistry() {
    }

    /**
     * Ordered structure → disordered parent (HSP) structure.
     * Add new ordered structures here. See class javadoc for instructions.
     */
    private static final Map<String, String> PARENT_OF = new LinkedHashMap<>();
    static {
        PARENT_OF.put("BCC_B2", "BCC_A2");
        PARENT_OF.put("FCC_L12", "FCC_A1");
        // BCC_D03 and FCC_L10 were removed: their cluster/symmetry files
        // produced degenerate identification results (tcf<=2, ncf<=1 —
        // effectively no usable non-point correlation functions) when run
        // through Stage 0-2, and were pulled from inputs/clus + inputs/sym.
        //
        // DIA_A4, CUB_B32, HEX_D019, ORC_B19, BCT_D022, CUB_C15-8p were also
        // removed: their -SG.txt files declare a non-identity transformation
        // matrix (rotation != identity and/or translation != 0), meaning
        // they are NOT standalone disordered (HSP) parents — the file
        // encodes a superlattice relationship to some other, larger
        // reference structure that was never identified. Self-mapping them
        // (treating them as their own parent) fails Stage 1b classification.
        // CUB_C15-8p additionally hangs (combinatorial blowup: 8-site
        // maximal cluster x 192 symmetry ops).
        //
        // Re-add only after re-deriving/confirming the correct parent
        // structure and re-verifying Stage 0-2 produces a sane, non-hanging
        // result.
    }

    /**
     * Returns the disordered (HSP) parent structure for the given (possibly
     * ordered) structure. Structures not present in {@link #PARENT_OF} are
     * their own parent (already disordered/HSP). Strips a trailing
     * {@code _CVCF} suffix before lookup.
     */
    public static String parentOf(String structure) {
        if (structure == null)
            return null;
        String base = structure.replace("_CVCF", "");
        return PARENT_OF.getOrDefault(base, base);
    }

    /** True if {@code structure} has a distinct disordered parent (i.e. is an ordered phase). */
    public static boolean isOrdered(String structure) {
        if (structure == null)
            return false;
        return PARENT_OF.containsKey(structure.replace("_CVCF", ""));
    }
}
