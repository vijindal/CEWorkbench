package org.ce.calculation.workflow;

import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEntry;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Extracts a genuine lower-order {@link CECEntry} for a subset of elements
 * from a higher-order system's CVCF Hamiltonian.
 *
 * <p>Per the CV-basis formalism this codebase implements (Jindal &amp; Lele
 * 2025, CALPHAD 89, 102825, Eq. 30), a cluster's CEC is inherited
 * <em>unchanged</em> into any higher-order system containing every element
 * that cluster references — no transformation is needed, only extraction of
 * the terms whose CVCF name refers exclusively to the requested element
 * subset, renamed to that subsystem's own letter assignment (A,B,... in the
 * subsystem's own element order).</p>
 *
 * <p>This is the engine behind {@link BinarySubsystemExtractor} (2-of-N),
 * used so a composition-grid scan can fall back to solving the true
 * lower-order CVM problem at a boundary of the parent composition space,
 * rather than evaluating the parent Hamiltonian at a zero-composition point
 * — CVM's cluster probabilities are continuous functions of composition,
 * and the parent solve can be numerically fragile near such a boundary
 * even though the physical answer only involves fewer elements (see
 * CLAUDE.md's near-edge non-convergence discussion). The true lower-order
 * problem has a much smaller configuration space and solves far more
 * reliably. Implemented for an arbitrary subset size (not just 2-of-N)
 * since the renaming/canonicalisation logic doesn't care how many elements
 * are kept.</p>
 */
final class SubsystemExtractor {

    private SubsystemExtractor() {}

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

    /** One element symbol: an uppercase letter optionally followed by a lowercase one. */
    private static final Pattern ELEMENT_SYMBOL = Pattern.compile("[A-Z][a-z]?");

    /**
     * Extracts the sub-Hamiltonian for {@code subsetElements} from
     * {@code parentEntry} (a higher-order CVCF {@link CECEntry}).
     *
     * @param parentEntry     the higher-order system's loaded Hamiltonian
     * @param parentElements  the parent system's canonical element order
     *                        (index defines the A,B,C,... letter assignment
     *                        CVCF term names use)
     * @param subsetElements  2 or more of the parent's elements, in the
     *                        order they should be assigned A,B,C,... in the
     *                        extracted subsystem (any subset, any order)
     * @return a new {@link CECEntry} for the subsystem, containing only the
     *         point + pair + triangle + tetrahedron terms that refer
     *         exclusively to elements in {@code subsetElements}, renamed to
     *         the subsystem's own letter assignment
     * @throws IllegalArgumentException if any element is not in
     *         {@code parentElements}, {@code subsetElements} has fewer than
     *         2 or more entries than {@code parentElements}, or no matching
     *         terms are found
     */
    static CECEntry extract(CECEntry parentEntry, List<String> parentElements, List<String> subsetElements) {
        if (subsetElements.size() < 2 || subsetElements.size() > parentElements.size()) {
            throw new IllegalArgumentException("Subset must have between 2 and " + parentElements.size()
                    + " elements, got " + subsetElements);
        }

        // Subsystem's own letter assignment: elements ordered by their position
        // in the PARENT system get consecutive letters A,B,C,... in that same
        // relative order, matching CvCfBasis's convention (elements ordered
        // A,B,C,... by position in the system's element string). The order of
        // subsetElements as passed in only decides output naming (elements
        // field), not which parent letter maps to which subsystem letter.
        List<Integer> parentIdx = new ArrayList<>(subsetElements.size());
        for (String el : subsetElements) {
            int idx = indexOfIgnoreCase(parentElements, el);
            if (idx < 0) {
                throw new IllegalArgumentException("Element '" + el
                        + "' not found in parent system elements " + parentElements);
            }
            parentIdx.add(idx);
        }
        List<Integer> sortedParentIdx = new ArrayList<>(parentIdx);
        sortedParentIdx.sort(null);

        // parentLetterToSubLetter[parent index] = subsystem letter, or '\0' if excluded.
        char[] parentLetterToSubLetter = new char[parentElements.size()];
        for (int i = 0; i < parentLetterToSubLetter.length; i++) parentLetterToSubLetter[i] = '\0';
        for (int subLetter = 0; subLetter < sortedParentIdx.size(); subLetter++) {
            parentLetterToSubLetter[sortedParentIdx.get(subLetter)] = (char) ('A' + subLetter);
        }

        StringBuilder subsystemName = new StringBuilder();
        for (int i = 0; i < sortedParentIdx.size(); i++) {
            if (i > 0) subsystemName.append('-');
            subsystemName.append(parentElements.get(sortedParentIdx.get(i)));
        }

        // A Hamiltonian may spell its terms three ways (see canonicalName);
        // resolve each to the internal letter form first, so the subset filter
        // below has exactly one shape to reason about.
        List<String> basisNames = basisNamesOrNull(parentEntry, parentElements.size());

        List<CECEntry.CECTerm> subTerms = new ArrayList<>();
        for (CECEntry.CECTerm term : parentEntry.cecTerms) {
            String canonical = canonicalName(term.name, parentElements, basisNames);
            String renamed = renameIfSubsetOnly(canonical, parentLetterToSubLetter);
            if (renamed == null) continue; // term references an excluded element, or isn't a multi-site term

            CECEntry.CECTerm copy = new CECEntry.CECTerm();
            copy.name = renamed;
            copy.description = term.description;
            copy.numSites = term.numSites;
            copy.multiplicity = term.multiplicity;
            copy.a = term.a;
            copy.b = term.b;
            subTerms.add(copy);
        }

        if (subTerms.isEmpty()) {
            String sample = parentEntry.cecTerms.length == 0 ? "(none)"
                    : parentEntry.cecTerms[0].name;
            throw new IllegalArgumentException("No matching CVCF terms found for subset "
                    + subsystemName + " in the parent Hamiltonian's cecTerms"
                    + " (first term is named '" + sample + "'; recognised spellings are"
                    + " the internal letter form v21AB, the paper's e2NbTi1, and"
                    + " positional CF_<index>).");
        }

        CECEntry sub = new CECEntry();
        sub.elements = subsystemName.toString();
        sub.structurePhase = parentEntry.structurePhase;
        sub.model = parentEntry.model;
        sub.cecUnits = parentEntry.cecUnits;
        sub.cecTerms = subTerms.toArray(new CECEntry.CECTerm[0]);
        sub.ncf = subTerms.size();
        sub.notes = "Subsystem [" + subsystemName + "] extracted from parent Hamiltonian "
                + parentEntry.elements + " (CECs are inherited unchanged per the CV-basis formalism).";
        return sub;
    }

    /**
     * If {@code termName}'s species-letter suffix refers exclusively to
     * letters present in {@code parentLetterToSubLetter} (any count of each -
     * covering pair/triangle/tetrahedron like- or unlike-species terms),
     * returns the term renamed to the subsystem's own letters. Returns null
     * if the term references any excluded letter (i.e. it's an interaction
     * outside the requested subset), doesn't match the expected
     * "prefix + LETTERS + optional trailing digits" shape, or is a
     * single-letter point CF (per-element, not per-cluster; the subsystem
     * derives its own point CFs from composition).
     */
    private static String renameIfSubsetOnly(String termName, char[] parentLetterToSubLetter) {
        Matcher m = TERM_NAME.matcher(termName);
        if (!m.matches()) return null;

        String prefix = m.group(1);
        String letters = m.group(2);
        String suffix = m.group(3);

        if (letters.length() == 1) return null; // point CF; skip

        StringBuilder renamed = new StringBuilder(prefix);
        for (int i = 0; i < letters.length(); i++) {
            int off = letters.charAt(i) - 'A';
            if (off < 0 || off >= parentLetterToSubLetter.length) return null;
            char subLetter = parentLetterToSubLetter[off];
            if (subLetter == '\0') return null; // references an excluded element
            renamed.append(subLetter);
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
     * <p>Returns the input unchanged when no rule applies; the caller's
     * subset filter then rejects it, as it would any unrecognised name.</p>
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

    private static int indexOfIgnoreCase(List<String> list, String s) {
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).equalsIgnoreCase(s)) return i;
        }
        return -1;
    }
}
