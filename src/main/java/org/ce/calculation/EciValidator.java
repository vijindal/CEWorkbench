package org.ce.calculation;

import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Strict validation of a Hamiltonian's ECI terms against a CVCF basis.
 *
 * <p><b>Why this exists.</b> {@link CECEvaluator} matches ECI terms to basis
 * correlation functions <em>by name</em>. A name that resolves to nothing is
 * silently skipped, leaving that interaction at {@code 0.0} while the calculation
 * still reports success — the only trace is a log warning. A human editing
 * {@code hamiltonian.json} might catch that; a machine API caller receiving
 * plausible-looking numbers from a partially-zeroed Hamiltonian will not.</p>
 *
 * <p>This validator runs <em>before</em> any calculation and reports both failure
 * modes: names that matched nothing ({@code unmatched}) and basis CFs that nothing
 * supplied ({@code unmapped}). It deliberately reuses
 * {@link CECEvaluator#buildLookupMap} and {@link CECEvaluator#resolveIndex} rather
 * than reimplementing the alias rules — a divergent second copy could accept a name
 * the evaluator later drops, recreating the exact hazard being guarded against.</p>
 */
public final class EciValidator {

    private EciValidator() {}

    /** Outcome of validating a Hamiltonian against a basis. */
    public static final class Result {
        /** Supplied term names that resolved to no basis CF. */
        public final List<String> unmatched;
        /** Basis (non-point) CF names that no supplied term mapped onto. */
        public final List<String> unmapped;
        /** All non-point CF names the basis expects, in basis order. */
        public final List<String> expected;

        Result(List<String> unmatched, List<String> unmapped, List<String> expected) {
            this.unmatched = List.copyOf(unmatched);
            this.unmapped  = List.copyOf(unmapped);
            this.expected  = List.copyOf(expected);
        }

        /** True when every supplied name matched and every basis CF was covered. */
        public boolean isValid() {
            return unmatched.isEmpty() && unmapped.isEmpty();
        }

        public String message() {
            StringBuilder sb = new StringBuilder("ECI validation failed.");
            if (!unmatched.isEmpty())
                sb.append(" Unmatched term names (not in basis): ").append(unmatched).append('.');
            if (!unmapped.isEmpty())
                sb.append(" Missing ECIs (basis CFs left unset, would default to 0.0): ")
                  .append(unmapped).append('.');
            sb.append(" Expected: ").append(expected).append('.');
            return sb.toString();
        }
    }

    /**
     * Validates that {@code cec}'s terms map exactly onto {@code basis}'s non-point CFs.
     *
     * @param cec   Hamiltonian to check; null or empty terms yields all-unmapped
     * @param basis target basis defining expected names and ordering
     */
    public static Result validate(CECEntry cec, CvCfBasis basis) {
        return validate(cec,
                basis.cfNames.subList(0, basis.numNonPointCfs),
                CECEvaluator.buildLookupMap(basis),
                basis.totalCfs());
    }

    /**
     * Validates against a bare list of expected non-point CF names, without needing a
     * fully generated {@link CvCfBasis}.
     *
     * <p>Used for pre-flight checks (e.g. the JSON API) where rejecting a bad
     * Hamiltonian before paying for basis generation is the point. The alias map is
     * built with {@link CECEvaluator#buildAliasMap} — the same routine the evaluator
     * uses — so acceptance here implies mapping there.</p>
     *
     * @param expectedNonPointCfNames from {@code CvCfBasis.getNonPointCfNames(...)}
     */
    public static Result validate(CECEntry cec, List<String> expectedNonPointCfNames) {
        return validate(cec,
                expectedNonPointCfNames,
                CECEvaluator.buildAliasMap(expectedNonPointCfNames),
                expectedNonPointCfNames.size());
    }

    private static Result validate(CECEntry cec, List<String> expected,
                                   Map<String, Integer> lookup, int maxTcf) {
        int ncf = expected.size();
        boolean[] mapped = new boolean[ncf];
        List<String> unmatched = new ArrayList<>();

        if (cec != null && cec.cecTerms != null) {
            for (CECEntry.CECTerm term : cec.cecTerms) {
                if (term.name == null) {
                    unmatched.add("<null name>");
                    continue;
                }
                Integer idx = CECEvaluator.resolveIndex(term.name, lookup, maxTcf);
                if (idx != null && idx >= 0 && idx < ncf) {
                    mapped[idx] = true;
                } else {
                    unmatched.add(term.name);
                }
            }
        }

        List<String> unmapped = new ArrayList<>();
        for (int i = 0; i < ncf; i++) {
            if (!mapped[i]) unmapped.add(expected.get(i));
        }

        return new Result(unmatched, unmapped, List.copyOf(expected));
    }
}
