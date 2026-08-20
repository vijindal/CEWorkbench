package org.ce.calculation.workflow;

import org.ce.model.hamiltonian.CECEntry;

import java.util.List;

/**
 * Extracts a genuine binary {@link CECEntry} for one unlike-element pair from
 * a higher-order (ternary+) system's CVCF Hamiltonian.
 *
 * <p>Thin wrapper over {@link SubsystemExtractor} (2-of-N case); see that
 * class for the extraction rationale and the three CVCF term-name spellings
 * it recognises.</p>
 *
 * <p>This exists so a ternary/quaternary composition-grid scan can fall back
 * to solving the actual binary CVM problem at a binary edge, rather than
 * evaluating the parent Hamiltonian at a zero-composition point.</p>
 */
public final class BinarySubsystemExtractor {

    private BinarySubsystemExtractor() {}

    /**
     * Extracts the binary sub-Hamiltonian for {@code elementA}/{@code elementB}
     * from {@code parentEntry} (a ternary+ CVCF {@link CECEntry}).
     *
     * @param parentEntry     the higher-order system's loaded Hamiltonian
     * @param parentElements  the parent system's canonical element order
     *                        (index defines the A,B,C,... letter assignment
     *                        CVCF term names use)
     * @param elementA        one of the pair's elements (any order)
     * @param elementB        the other
     * @return a new {@link CECEntry} for the 2-component system, containing
     *         only the point + pair + triangle + tetrahedron terms that refer
     *         exclusively to this pair, renamed to the binary basis's own
     *         letter assignment
     * @throws IllegalArgumentException if either element is not in
     *         {@code parentElements}, or no matching terms are found
     */
    public static CECEntry extractBinary(CECEntry parentEntry, List<String> parentElements,
            String elementA, String elementB) {
        return SubsystemExtractor.extract(parentEntry, parentElements, List.of(elementA, elementB));
    }
}
