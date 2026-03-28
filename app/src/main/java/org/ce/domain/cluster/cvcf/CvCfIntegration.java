package org.ce.domain.cluster.cvcf;

import org.ce.domain.cluster.CMatrixResult;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Encapsulates CVCF basis integration for CVMPhaseModel.
 *
 * <p>Handles:
 * <ul>
 *   <li>Transformation of C-matrix from old orthogonal basis to CVCF basis</li>
 *   <li>Extraction of component pairs, triples, and quads (alphabetically ordered)</li>
 *   <li>Mapping between flat ECI array and CVCF CF name-based ECI map</li>
 * </ul>
 */
public final class CvCfIntegration {

    private final CvCfBasis basis;
    private final CMatrixResult transformedCmat;
    private final List<String[]> componentPairs;
    private final List<String[]> componentTriples;
    private final List<String[]> componentQuads;

    /**
     * Constructor (package-private for testing).
     */
    CvCfIntegration(
            CvCfBasis basis,
            CMatrixResult transformedCmat,
            List<String[]> componentPairs,
            List<String[]> componentTriples,
            List<String[]> componentQuads) {

        this.basis = basis;
        this.transformedCmat = transformedCmat;
        this.componentPairs = componentPairs;
        this.componentTriples = componentTriples;
        this.componentQuads = componentQuads;
    }

    /**
     * Creates a CVCF integration for binary BCC_A2.
     *
     * @param oldCmatResult the C-matrix in old orthogonal basis
     * @param componentA    name of component A (e.g., "Nb")
     * @param componentB    name of component B (e.g., "Ti")
     * @return CvCfIntegration with transformed data
     */
    public static CvCfIntegration forBccA2Binary(
            CMatrixResult oldCmatResult,
            String componentA,
            String componentB) {

        // Components in alphabetical order
        String[] sorted = sortComponents(componentA, componentB);
        CvCfBasis basis = BccA2TModelCvCfTransformations.binaryBasis();
        CMatrixResult transformedCmat = CvCfBasisTransformer.transform(oldCmatResult, basis);

        // Component pairs: only [A, B]
        List<String[]> pairs = new ArrayList<>();
        pairs.add(new String[]{sorted[0], sorted[1]});

        return new CvCfIntegration(basis, transformedCmat, pairs, new ArrayList<>(), new ArrayList<>());
    }

    /**
     * Creates a CVCF integration for ternary BCC_A2.
     *
     * @param oldCmatResult the C-matrix in old orthogonal basis
     * @param componentA    name of component A
     * @param componentB    name of component B
     * @param componentC    name of component C
     * @return CvCfIntegration with transformed data
     */
    public static CvCfIntegration forBccA2Ternary(
            CMatrixResult oldCmatResult,
            String componentA,
            String componentB,
            String componentC) {

        // Components in alphabetical order
        String[] sorted = sortComponents(componentA, componentB, componentC);
        CvCfBasis basis = BccA2TModelCvCfTransformations.ternaryBasis();
        CMatrixResult transformedCmat = CvCfBasisTransformer.transform(oldCmatResult, basis);

        // Generate all component pairs (P < Q)
        List<String[]> pairs = new ArrayList<>();
        pairs.add(new String[]{sorted[0], sorted[1]});
        pairs.add(new String[]{sorted[0], sorted[2]});
        pairs.add(new String[]{sorted[1], sorted[2]});

        // Generate all component triples (P < Q < R)
        List<String[]> triples = new ArrayList<>();
        triples.add(new String[]{sorted[0], sorted[1], sorted[2]});

        return new CvCfIntegration(basis, transformedCmat, pairs, triples, new ArrayList<>());
    }

    /**
     * Creates a CVCF integration for quaternary BCC_A2.
     *
     * @param oldCmatResult the C-matrix in old orthogonal basis
     * @param componentA    name of component A
     * @param componentB    name of component B
     * @param componentC    name of component C
     * @param componentD    name of component D
     * @return CvCfIntegration with transformed data
     */
    public static CvCfIntegration forBccA2Quaternary(
            CMatrixResult oldCmatResult,
            String componentA,
            String componentB,
            String componentC,
            String componentD) {

        // Components in alphabetical order
        String[] sorted = sortComponents(componentA, componentB, componentC, componentD);
        CvCfBasis basis = BccA2TModelCvCfTransformations.quaternaryBasis();
        CMatrixResult transformedCmat = CvCfBasisTransformer.transform(oldCmatResult, basis);

        // Generate all component pairs (P < Q)
        List<String[]> pairs = new ArrayList<>();
        pairs.add(new String[]{sorted[0], sorted[1]});
        pairs.add(new String[]{sorted[0], sorted[2]});
        pairs.add(new String[]{sorted[0], sorted[3]});
        pairs.add(new String[]{sorted[1], sorted[2]});
        pairs.add(new String[]{sorted[1], sorted[3]});
        pairs.add(new String[]{sorted[2], sorted[3]});

        // Generate all component triples (P < Q < R)
        List<String[]> triples = new ArrayList<>();
        triples.add(new String[]{sorted[0], sorted[1], sorted[2]});
        triples.add(new String[]{sorted[0], sorted[1], sorted[3]});
        triples.add(new String[]{sorted[0], sorted[2], sorted[3]});
        triples.add(new String[]{sorted[1], sorted[2], sorted[3]});

        // Generate all component quads (P < Q < R < S)
        List<String[]> quads = new ArrayList<>();
        quads.add(new String[]{sorted[0], sorted[1], sorted[2], sorted[3]});

        return new CvCfIntegration(basis, transformedCmat, pairs, triples, quads);
    }

    // =========================================================================
    // Public accessors
    // =========================================================================

    public CvCfBasis getBasis() {
        return basis;
    }

    public CMatrixResult getTransformedCmat() {
        return transformedCmat;
    }

    public List<String[]> getComponentPairs() {
        return componentPairs;
    }

    public List<String[]> getComponentTriples() {
        return componentTriples;
    }

    public List<String[]> getComponentQuads() {
        return componentQuads;
    }

    /**
     * Converts a flat ECI array (CVCF basis order) to a Map with CVCF CF names.
     *
     * <p>The input array must have length ≥ numNonPointCfs. The ECI values
     * are mapped to CF names in order: eci[i] corresponds to the i-th
     * non-point CF name from basis.cfNames.</p>
     *
     * @param flatEciArray the ECI array in CVCF basis order (length ≥ numNonPointCfs)
     * @return Map from CVCF CF name to ECI value (only non-point CFs)
     */
    public Map<String, Double> flatEciToMap(double[] flatEciArray) {
        Map<String, Double> result = new HashMap<>();

        int numNonPoint = basis.numNonPointCfs;
        if (flatEciArray == null || flatEciArray.length < numNonPoint) {
            throw new IllegalArgumentException(
                    "ECI array too short: got " + (flatEciArray == null ? 0 : flatEciArray.length)
                    + ", expected >= " + numNonPoint);
        }

        // Map each non-point ECI value to its CF name
        List<String> cfNames = basis.cfNames;
        for (int i = 0; i < numNonPoint; i++) {
            String cfName = cfNames.get(i);
            result.put(cfName, flatEciArray[i]);
        }

        return result;
    }

    /**
     * Converts an ECI Map (by CF name) to a flat array in CVCF basis order.
     *
     * <p>Inverse of {@link #flatEciToMap}. Missing ECIs default to 0 (CEC inheritance).</p>
     *
     * @param eciMap map from CF name to ECI value
     * @return flat ECI array with length = numNonPointCfs
     */
    public double[] mapEciToFlat(Map<String, Double> eciMap) {
        double[] result = new double[basis.numNonPointCfs];
        List<String> cfNames = basis.cfNames;

        for (int i = 0; i < basis.numNonPointCfs; i++) {
            String cfName = cfNames.get(i);
            result[i] = eciMap.getOrDefault(cfName, 0.0);
        }

        return result;
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    /** Sorts component names alphabetically. */
    private static String[] sortComponents(String... names) {
        String[] sorted = names.clone();
        java.util.Arrays.sort(sorted);
        return sorted;
    }
}
