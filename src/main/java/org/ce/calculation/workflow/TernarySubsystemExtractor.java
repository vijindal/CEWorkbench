package org.ce.calculation.workflow;

import org.ce.model.hamiltonian.CECEntry;

import java.util.List;

/**
 * Extracts a genuine ternary {@link CECEntry} for three elements from a
 * higher-order (quaternary+) system's CVCF Hamiltonian.
 *
 * <p>Thin wrapper over {@link SubsystemExtractor} (3-of-N case); see that
 * class for the extraction rationale and the three CVCF term-name spellings
 * it recognises. Sibling of {@link BinarySubsystemExtractor} (2-of-N).</p>
 */
public final class TernarySubsystemExtractor {

    private TernarySubsystemExtractor() {}

    /**
     * Extracts the ternary sub-Hamiltonian for {@code elementA}/{@code elementB}/
     * {@code elementC} from {@code parentEntry} (a quaternary+ CVCF {@link CECEntry}).
     *
     * @param parentEntry     the higher-order system's loaded Hamiltonian
     * @param parentElements  the parent system's canonical element order
     *                        (index defines the A,B,C,... letter assignment
     *                        CVCF term names use)
     * @param elementA        one of the triple's elements (any order)
     * @param elementB        another
     * @param elementC        the third
     * @return a new {@link CECEntry} for the 3-component system, containing
     *         only the point + pair + triangle + tetrahedron terms that refer
     *         exclusively to this triple, renamed to the ternary basis's own
     *         letter assignment
     * @throws IllegalArgumentException if any element is not in
     *         {@code parentElements}, or no matching terms are found
     */
    public static CECEntry extractTernary(CECEntry parentEntry, List<String> parentElements,
            String elementA, String elementB, String elementC) {
        return SubsystemExtractor.extract(parentEntry, parentElements,
                List.of(elementA, elementB, elementC));
    }
}
