package org.ce.calculation.workflow;

import org.ce.model.cvm.CvCfBasis;
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
     * Positional name, e.g. {@code CF_45} — the spelling the Nb-Ti-V-Zr
     * Hamiltonian actually stores. {@code CECEvaluator} resolves these as an
     * index into the basis, and so does this class.
     */
    private static final Pattern CF_INDEX = Pattern.compile("^(?i:cf)_(\\d+)$");

    /**
     * The paper's own concrete spelling: {@code e2<Pair><shell>} with element
     * symbols rather than letters, e.g. {@code e2NbTi1}, {@code e3NbZr},
     * {@code e4VZr} (Jindal &amp; Lele 2025, Table 17, "CECs for the
     * Nb-Ti-V-Zr system in the CVCF basis"). Note the shell trails the pair
     * here, where our internal names put it first ({@code v21AB}).
     */
    private static final Pattern ELEMENT_NAMED =
            Pattern.compile("^([a-zA-Z]+?)(\\d?)((?:[A-Z][a-z]?)+)(\\d*)$");

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

        // A Hamiltonian may spell its terms three ways (see canonicalName);
        // resolve each to the internal letter form first, so the pair filter
        // below has exactly one shape to reason about.
        List<String> basisNames = basisNamesOrNull(parentEntry, parentElements.size());

        List<CECEntry.CECTerm> binaryTerms = new ArrayList<>();
        for (CECEntry.CECTerm term : parentEntry.cecTerms) {
            String canonical = canonicalName(term.name, parentElements, basisNames);
            String renamed = renameIfPairOnly(canonical, parentLetterLo, parentLetterHi);
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
            String sample = parentEntry.cecTerms.length == 0 ? "(none)"
                    : parentEntry.cecTerms[0].name;
            throw new IllegalArgumentException("No binary-pair CVCF terms found for ["
                    + elementA + "," + elementB + "] in the parent Hamiltonian's cecTerms"
                    + " (first term is named '" + sample + "'; recognised spellings are"
                    + " the internal letter form v21AB, the paper's e2NbTi1, and"
                    + " positional CF_<index>).");
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


    /**
     * The parent basis's canonical CF names, or null if the combination is not
     * registered. Only needed to resolve positional {@code CF_<index>} names.
     */
    private static List<String> basisNamesOrNull(CECEntry entry, int numComponents) {
        try {
            String model = entry.model == null ? "T" : entry.model.replace("_CVCF", "");
            return CvCfBasis.getNonPointCfNames(entry.structurePhase, model, numComponents);
        } catch (RuntimeException e) {
            return null;   // unregistered combination; letter/element names still work
        }
    }

    /**
     * Resolves a stored term name to the internal letter form ({@code v21AB}).
     *
     * <p>Three spellings occur in practice and all mean the same cluster:</p>
     * <ul>
     *   <li><b>Internal letter form</b> — {@code v21AB}, {@code e21AB},
     *       {@code v3AD}: species as A,B,C,... letters, shell first. Returned
     *       unchanged.</li>
     *   <li><b>Positional</b> — {@code CF_45}: an index into the basis, which
     *       is what the Nb-Ti-V-Zr Hamiltonian stores. Resolved through
     *       {@link CvCfBasis} so this class and {@code CECEvaluator} cannot
     *       disagree about what index 45 means.</li>
     *   <li><b>The paper's element-named form</b> — {@code e2NbTi1},
     *       {@code e3NbZr} (Table 17). Element symbols instead of letters, and
     *       the shell <em>trails</em> the pair. Translated back to letters via
     *       the parent's element order.</li>
     * </ul>
     *
     * <p>Returns the input unchanged when no rule applies; the caller's pair
     * filter then rejects it, as it would any unrecognised name.</p>
     */
    private static String canonicalName(String name, List<String> parentElements,
            List<String> basisNames) {
        if (name == null) return "";

        Matcher cf = CF_INDEX.matcher(name);
        if (cf.matches() && basisNames != null) {
            int idx = Integer.parseInt(cf.group(1));
            if (idx >= 0 && idx < basisNames.size()) return basisNames.get(idx);
            return name;
        }

        // Already in letter form? Then every letter is a valid species letter.
        Matcher letters = TERM_NAME.matcher(name);
        if (letters.matches() && allSpeciesLetters(letters.group(2), parentElements.size())) {
            return name;
        }

        Matcher en = ELEMENT_NAMED.matcher(name);
        if (en.matches()) {
            String prefix = en.group(1);          // "e" / "v"
            String order = en.group(2);           // cluster order digit, e.g. "2"
            String symbols = en.group(3);         // "NbTi"
            String shell = en.group(4);           // trailing shell, e.g. "1"

            StringBuilder out = new StringBuilder(prefix).append(order).append(shell);
            StringBuilder mapped = new StringBuilder();
            for (Matcher sm = ELEMENT_SYMBOL.matcher(symbols); sm.find(); ) {
                int i = indexOfIgnoreCase(parentElements, sm.group());
                if (i < 0) return name;           // not an element of this system
                mapped.append((char) ('A' + i));
            }
            if (mapped.length() == 0) return name;
            return out.append(mapped).toString();
        }
        return name;
    }

    /** True if every character is a species letter within this system's range. */
    private static boolean allSpeciesLetters(String letters, int numComponents) {
        for (int i = 0; i < letters.length(); i++) {
            int off = letters.charAt(i) - 'A';
            if (off < 0 || off >= numComponents) return false;
        }
        return true;
    }

    /** One element symbol: an uppercase letter optionally followed by a lowercase one. */
    private static final Pattern ELEMENT_SYMBOL = Pattern.compile("[A-Z][a-z]?");

    private static int indexOfIgnoreCase(List<String> list, String s) {
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).equalsIgnoreCase(s)) return i;
        }
        return -1;
    }
}
