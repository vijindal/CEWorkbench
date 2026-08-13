package org.ce.calculation.workflow;

import org.ce.model.hamiltonian.CECEntry;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Extracts a genuine binary {@link CECEntry} for one unlike-element pair from
 * a higher-order (ternary+) system's CVCF Hamiltonian.
 *
 * <p>Per the CV-basis formalism this codebase implements (Jindal &amp; Lele
 * 2025, CALPHAD 89, 102825, Eq. 30), binary-cluster CECs are inherited
 * <em>unchanged</em> into any higher-order system containing that pair — no
 * transformation is needed, only extraction of the terms whose CVCF name
 * refers exclusively to the two requested elements (e.g. {@code e21AB},
 * {@code e22AB}, {@code e3AB}, {@code e4AB} for the pair mapped to letters
 * A,B in the parent system), renamed to the binary system's own letter
 * assignment (always A,B for a 2-component system).</p>
 *
 * <p>This exists so a ternary composition-grid scan can fall back to solving
 * the actual binary CVM problem at a binary edge, rather than evaluating the
 * ternary Hamiltonian at a zero-composition point — the latter still runs
 * the full ternary cluster machinery and can be numerically fragile there
 * (see CLAUDE.md's near-edge non-convergence discussion). A true binary
 * solve has a much smaller configuration space and is far more robust.</p>
 */
public final class BinarySubsystemExtractor {

    private BinarySubsystemExtractor() {}

    /** CVCF term-name pattern: a cluster-type prefix followed by 1-4 letters (species). */
    private static final Pattern TERM_NAME = Pattern.compile("^([a-zA-Z0-9]*?)([A-Z]+)(\\d*)$");

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

        int idxA = indexOfIgnoreCase(parentElements, elementA);
        int idxB = indexOfIgnoreCase(parentElements, elementB);
        if (idxA < 0 || idxB < 0) {
            throw new IllegalArgumentException("Pair [" + elementA + "," + elementB
                    + "] not found in parent system elements " + parentElements);
        }

        // Binary system's own alphabetical letter assignment: the lower parent
        // index gets 'A', matching CvCfBasis's convention (elements ordered
        // A,B,C,... by their position in the system's element string).
        int loIdx = Math.min(idxA, idxB), hiIdx = Math.max(idxA, idxB);
        char parentLetterLo = (char) ('A' + loIdx);
        char parentLetterHi = (char) ('A' + hiIdx);
        String parentLower = parentElements.get(loIdx);
        String parentHigher = parentElements.get(hiIdx);

        List<CECEntry.CECTerm> binaryTerms = new ArrayList<>();
        for (CECEntry.CECTerm term : parentEntry.cecTerms) {
            String renamed = renameIfPairOnly(term.name, parentLetterLo, parentLetterHi);
            if (renamed == null) continue; // term references a 3rd/4th element, or isn't a pair-only term

            CECEntry.CECTerm copy = new CECEntry.CECTerm();
            copy.name = renamed;
            copy.description = term.description;
            copy.numSites = term.numSites;
            copy.multiplicity = term.multiplicity;
            copy.a = term.a;
            copy.b = term.b;
            binaryTerms.add(copy);
        }

        if (binaryTerms.isEmpty()) {
            throw new IllegalArgumentException("No binary-pair CVCF terms found for ["
                    + elementA + "," + elementB + "] in the parent Hamiltonian's cecTerms.");
        }

        CECEntry binary = new CECEntry();
        binary.elements = parentLower + "-" + parentHigher;
        binary.structurePhase = parentEntry.structurePhase;
        binary.model = parentEntry.model;
        binary.cecUnits = parentEntry.cecUnits;
        binary.cecTerms = binaryTerms.toArray(new CECEntry.CECTerm[0]);
        binary.ncf = binaryTerms.size();
        binary.notes = "Binary subsystem [" + parentLower + "-" + parentHigher
                + "] extracted from parent Hamiltonian " + parentEntry.elements
                + " (binary CECs are inherited unchanged per the CV-basis formalism).";
        return binary;
    }

    /**
     * If {@code termName}'s species-letter suffix refers exclusively to
     * {@code loLetter}/{@code hiLetter} (in either order, any count of each -
     * covering pair/triangle/tetrahedron like-pair or unlike-pair terms),
     * returns the term renamed to the binary basis (lo -> 'A', hi -> 'B').
     * Returns null if the term references any other letter (i.e. it's a
     * ternary/quaternary-only interaction) or doesn't match the expected
     * "prefix + LETTERS + optional trailing digits" shape.
     */
    private static String renameIfPairOnly(String termName, char loLetter, char hiLetter) {
        Matcher m = TERM_NAME.matcher(termName);
        if (!m.matches()) return null;

        String prefix = m.group(1);
        String letters = m.group(2);
        String suffix = m.group(3);

        for (int i = 0; i < letters.length(); i++) {
            char c = letters.charAt(i);
            if (c != loLetter && c != hiLetter) return null; // references a 3rd/4th element
        }
        // Point CFs (e.g. "xA") are per-element, not per-pair; skip them here -
        // the binary basis derives its own point CFs from composition.
        if (letters.length() == 1) return null;

        StringBuilder renamed = new StringBuilder(prefix);
        for (int i = 0; i < letters.length(); i++) {
            renamed.append(letters.charAt(i) == loLetter ? 'A' : 'B');
        }
        renamed.append(suffix);
        return renamed.toString();
    }

    private static int indexOfIgnoreCase(List<String> list, String s) {
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).equalsIgnoreCase(s)) return i;
        }
        return -1;
    }
}
