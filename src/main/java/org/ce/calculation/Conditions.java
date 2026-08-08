package org.ce.calculation;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A single-point condition set: temperature + composition by element symbol.
 *
 * <p>Mirrors pycalphad's {@code conditions} dict ({@code {v.T: 1000, v.X('Ti'): 0.5}}) —
 * composition is named, never positional. At most K-1 of the K elements need be
 * given; the omitted one is derived so the composition sums to 1.0.</p>
 */
public record Conditions(double temperature, Map<String, Double> composition) {

    public Conditions {
        composition = Map.copyOf(composition);
    }

    public static Conditions of(double temperature, Map<String, Double> composition) {
        return new Conditions(temperature, composition);
    }

    /**
     * Resolves composition against the canonical element order to a full-K array.
     * The element not present in {@link #composition} is derived as {@code 1 - sum(given)}.
     * Throws on unknown element names, out-of-range fractions, duplicate entries
     * (impossible via Map, kept for symmetry with error messages), or an
     * under/fully-determined system that doesn't resolve to exactly one derived value.
     */
    public double[] resolveComposition(List<String> canonicalOrder) {
        int K = canonicalOrder.size();
        if (K == 0)
            throw new IllegalStateException("Empty canonical element order");
        if (composition.isEmpty())
            throw new IllegalArgumentException(
                    "No composition specified. Provide at least " + (K - 1) + " of " + K
                    + " element=fraction pairs for system " + canonicalOrder + ".");

        double[] full = new double[K];
        boolean[] given = new boolean[K];

        for (var e : composition.entrySet()) {
            int idx = indexOfIgnoreCase(canonicalOrder, e.getKey());
            if (idx < 0)
                throw new IllegalArgumentException(
                        "Unknown element '" + e.getKey() + "' for system " + canonicalOrder
                        + ". Valid: " + canonicalOrder);
            double v = e.getValue();
            if (Double.isNaN(v) || v < 0.0 || v > 1.0)
                throw new IllegalArgumentException(
                        "Mole fraction for '" + e.getKey() + "' must be in [0,1], got " + v);
            full[idx] = v;
            given[idx] = true;
        }

        int missing = 0, missingIdx = -1;
        for (int i = 0; i < K; i++) if (!given[i]) { missing++; missingIdx = i; }
        double sum = 0;
        for (int i = 0; i < K; i++) if (given[i]) sum += full[i];

        if (missing == 0) {
            if (Math.abs(sum - 1.0) > 1e-6)
                throw new IllegalArgumentException(
                        "Fully-specified composition must sum to 1.0, got " + sum + ": " + composition);
        } else if (missing == 1) {
            if (sum > 1.0 + 1e-6)
                throw new IllegalArgumentException(
                        "Specified fractions sum to " + sum + " > 1; no room for derived element '"
                        + canonicalOrder.get(missingIdx) + "'");
            full[missingIdx] = 1.0 - sum;
        } else {
            List<String> missingNames = new ArrayList<>();
            for (int i = 0; i < K; i++) if (!given[i]) missingNames.add(canonicalOrder.get(i));
            throw new IllegalArgumentException(
                    "Underdetermined: " + missing + " elements unspecified " + missingNames
                    + ", need at most 1 unspecified. Specify " + (K - 1) + " of " + K + " fractions.");
        }
        return full;
    }

    private static int indexOfIgnoreCase(List<String> list, String sym) {
        for (int i = 0; i < list.size(); i++)
            if (list.get(i).equalsIgnoreCase(sym)) return i;
        return -1;
    }
}
