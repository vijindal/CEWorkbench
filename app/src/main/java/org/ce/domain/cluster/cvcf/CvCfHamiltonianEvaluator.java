package org.ce.domain.cluster.cvcf;

import java.util.List;
import java.util.Map;

/**
 * Computes the mixing enthalpy Hmc directly in the CVCF basis (eq. 30 of
 * Jindal &amp; Lele 2025, CALPHAD).
 *
 * <p>For BCC_A2, Hmc is the sum of binary, ternary, and quaternary interaction
 * terms, where each term is an effective ECI (multiplicity absorbed) times a
 * cluster variable value:</p>
 *
 * <pre>
 * Hmc = Σ_{P&lt;Q}       (e21PQ·v21PQ + e22PQ·v22PQ + e3PQ·v3PQ + e4PQ·v4PQ)
 *     + Σ_{P&lt;Q&lt;R}     (e3PQR1·v3PQR1 + e3PQR2·v3PQR2 + e3PQR3·v3PQR3
 *                     + e4PQR1·v4PQR1 + e4PQR2·v4PQR2 + e4PQR3·v4PQR3)
 *     + Σ_{P&lt;Q&lt;R&lt;S}  (e4PQRS1·v4PQRS1 + e4PQRS2·v4PQRS2 + e4PQRS3·v4PQRS3)
 * </pre>
 *
 * <p>Effective ECIs already absorb multiplicity (user provides them directly).
 * Missing ECIs default to 0 — this is the CEC inheritance mechanism: binary
 * CECs are used as-is in ternary/quaternary systems with zero for higher-order
 * interaction terms unless explicitly provided.</p>
 */
public final class CvCfHamiltonianEvaluator {

    private CvCfHamiltonianEvaluator() {}

    /**
     * Computes Hmc from CVCF cluster variable values and effective ECIs.
     *
     * @param cvCfValues  map from CF name (e.g., "v21NbTi") to its current CV value
     * @param ecis        map from CF name to effective ECI (multiplicity absorbed)
     * @param componentPairs   ordered pairs [P,Q] with P &lt; Q (alphabetical)
     * @param componentTriples ordered triples [P,Q,R] with P &lt; Q &lt; R
     * @param componentQuads   ordered quads [P,Q,R,S] with P &lt; Q &lt; R &lt; S
     * @return Hmc in the same units as the ECIs
     */
    public static double computeH(
            Map<String, Double> cvCfValues,
            Map<String, Double> ecis,
            List<String[]> componentPairs,
            List<String[]> componentTriples,
            List<String[]> componentQuads) {

        double H = 0.0;

        // Binary terms: Σ_{P<Q} (e21PQ·v21PQ + e22PQ·v22PQ + e3PQ·v3PQ + e4PQ·v4PQ)
        for (String[] pq : componentPairs) {
            String suffix = pq[0] + pq[1];
            H += eci(ecis, "v21" + suffix) * cv(cvCfValues, "v21" + suffix);
            H += eci(ecis, "v22" + suffix) * cv(cvCfValues, "v22" + suffix);
            H += eci(ecis, "v3"  + suffix) * cv(cvCfValues, "v3"  + suffix);
            H += eci(ecis, "v4"  + suffix) * cv(cvCfValues, "v4"  + suffix);
        }

        // Ternary terms: Σ_{P<Q<R} (triangle1,2,3 + tetra1,2,3 terms)
        for (String[] pqr : componentTriples) {
            String suffix = pqr[0] + pqr[1] + pqr[2];
            H += eci(ecis, "v3"  + suffix + "1") * cv(cvCfValues, "v3"  + suffix + "1");
            H += eci(ecis, "v3"  + suffix + "2") * cv(cvCfValues, "v3"  + suffix + "2");
            H += eci(ecis, "v3"  + suffix + "3") * cv(cvCfValues, "v3"  + suffix + "3");
            H += eci(ecis, "v4"  + suffix + "1") * cv(cvCfValues, "v4"  + suffix + "1");
            H += eci(ecis, "v4"  + suffix + "2") * cv(cvCfValues, "v4"  + suffix + "2");
            H += eci(ecis, "v4"  + suffix + "3") * cv(cvCfValues, "v4"  + suffix + "3");
        }

        // Quaternary terms: Σ_{P<Q<R<S} (tetra1,2,3 terms)
        for (String[] pqrs : componentQuads) {
            String suffix = pqrs[0] + pqrs[1] + pqrs[2] + pqrs[3];
            H += eci(ecis, "v4" + suffix + "1") * cv(cvCfValues, "v4" + suffix + "1");
            H += eci(ecis, "v4" + suffix + "2") * cv(cvCfValues, "v4" + suffix + "2");
            H += eci(ecis, "v4" + suffix + "3") * cv(cvCfValues, "v4" + suffix + "3");
        }

        return H;
    }

    /** Returns the ECI for the given name, or 0 if not present (CEC inheritance). */
    private static double eci(Map<String, Double> ecis, String name) {
        return ecis.getOrDefault(name, 0.0);
    }

    /** Returns the CV value for the given name, or 0 if not present. */
    private static double cv(Map<String, Double> cvCfValues, String name) {
        return cvCfValues.getOrDefault(name, 0.0);
    }
}
