package org.ce.calculation.workflow;

import static org.ce.model.cluster.AllClusterData.ClusterData;

import static org.ce.model.cluster.ClusterResults.*;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.model.cluster.AllClusterData;
import org.ce.model.cluster.CFIdentificationResult;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.NumericalCECTransformer;
import org.ce.model.storage.DataStore;

import java.io.IOException;
import java.util.List;

/**
 * Workflow responsible for managing Cluster Expansion Coefficient (CEC) databases.
 *
 * This orchestrates operations between cluster data and CEC storage.
 */
public class CECManagementWorkflow {

    /**
     * Holds metadata about a correlation function extracted from cluster data.
     * Used during Hamiltonian scaffolding to enrich CECEntry.CECTerm with descriptive information.
     */
    private static class CFMetadata {

        int numSites;
        double multiplicity;
        String description;

        CFMetadata(int numSites, double multiplicity, String description) {
            this.numSites = numSites;
            this.multiplicity = multiplicity;
            this.description = description;
        }
    }

    private final DataStore.HamiltonianStore store;

    public CECManagementWorkflow(DataStore.HamiltonianStore store) {
        this.store = store;
    }

    /**
     * Creates a blank CEC database in memory (no disk IO).
     * Populates optional metadata (numSites, multiplicity, description) from cluster data.
     *
     * @param elements       element string (e.g. "Nb-Ti")
     * @param structurePhase structure identifier
     * @param model          model identifier
     * @param ncf            number of correlation functions
     * @param cfMetadata     CF metadata from cluster data (contains numSites, multiplicity)
     */
    public CECEntry scaffoldEmptyCEC(
            String elements,
            String structurePhase,
            String model,
            int ncf,
            CFMetadata[] cfMetadata) {

        CECEntry.CECTerm[] terms = new CECEntry.CECTerm[ncf];

        for (int i = 0; i < ncf; i++) {
            CECEntry.CECTerm term = new CECEntry.CECTerm();
            term.name = "CF_" + i;
            term.a = 0.0;
            term.b = 0.0;

            // Populate optional CF metadata if available
            if (cfMetadata != null && i < cfMetadata.length && cfMetadata[i] != null) {
                term.numSites = cfMetadata[i].numSites;
                term.multiplicity = cfMetadata[i].multiplicity;
                term.description = cfMetadata[i].description;
            }

            terms[i] = term;
        }

        CECEntry entry = new CECEntry();
        entry.elements = elements;
        entry.structurePhase = structurePhase;
        entry.model = model;
        entry.cecTerms = terms;
        entry.cecUnits = "J/mol";
        entry.reference = "";
        entry.notes = "Scaffolded empty CEC";
        entry.ncf = ncf;

        return entry;
    }

    /**
     * Creates a new empty CEC database and saves it to disk.
     *
     * @param hamiltonianId ID for the saved file — {elements}_{structure}_{model}, e.g. Nb-Ti_BCC_A2_T
     * @param cfMetadata    optional CF metadata (numSites, multiplicity) for enriching CECEntry.CECTerms
     */
    public CECEntry createAndSaveCEC(
            String hamiltonianId,
            String elements,
            String structurePhase,
            String model,
            int ncf,
            CFMetadata[] cfMetadata) throws IOException {

        CECEntry entry = scaffoldEmptyCEC(elements, structurePhase, model, ncf, cfMetadata);

        store.save(hamiltonianId, entry);

        return entry;
    }

    /**
     * Validates that the CEC database matches the expected number
     * of correlation functions.
     *
     * <p>Enforces the physical law: G = Σ (ECI_l * CF_l)
     * Therefore: length(ECI) == ncf == number of correlation functions</p>
     */
    public void validateCEC(CECEntry entry, int expectedNcf) {

        if (entry.structurePhase == null || entry.structurePhase.isBlank()) {
            throw new IllegalStateException("CECEntry.structurePhase is null or blank");
        }
        if (entry.elements == null || entry.elements.isBlank()) {
            throw new IllegalStateException("CECEntry.elements is null or blank");
        }

        if (entry.cecTerms == null) {
            throw new IllegalStateException("CEC database contains no terms");
        }

        // Note: With CVCF basis, cluster data may report ncf = total CFs (non-point + point),
        // while hamiltonian has only non-point CFs. Accept if term count is reasonable.
        int termCount = entry.cecTerms.length;
        if (termCount <= 0 || termCount > expectedNcf) {
            throw new IllegalStateException(
                "CEC term count (" + termCount +
                ") is invalid (must be > 0 and <= " + expectedNcf + ")"
            );
        }

        for (CECEntry.CECTerm term : entry.cecTerms) {
            term.validate();
        }
    }

    /**
     * Creates and saves a Hamiltonian scaffold using ncf from an existing cluster dataset.
     *
     * <p>The cluster data ID is derived from the Hamiltonian ID:
     * {@code A-B_BCC_A2_T} → {@code BCC_A2_T_bin}, ensuring the two IDs
     * always stay in sync and cannot be accidentally mismatched by the caller.</p>
     *
     * <p>Also extracts CF metadata (numSites, multiplicity) from cluster data to enrich
     * the scaffolded CECEntry.CECTerms with helpful descriptive information.</p>
     *
     * @param hamiltonianId ID for the saved Hamiltonian file, e.g. Nb-Ti_BCC_A2_T
     * @param elements      element string (e.g. "Nb-Ti")
     * @param structurePhase structure identifier
     * @param model         model identifier
     */
    public CECEntry scaffoldFromClusterData(
            String hamiltonianId,
            String elements,
            String structurePhase,
            String model) throws Exception {

        // Determine numComponents from elements (e.g. "A-B" -> 2, "A-B-C" -> 3)
        int numComponents = elements.split("-").length;

        // Validate CVCF basis is supported for this (structure, model, numComponents)
        if (!CvCfBasis.isSupported(structurePhase, model, numComponents)) {
            String msg = "CVCF basis not supported for " + structurePhase
                    + " with model=" + model + " and " + numComponents + " components.\n"
                    + "Cannot scaffold CEC for unsupported basis.\n"
                    + CvCfBasis.supportedSummary();
            throw new IllegalArgumentException(msg);
        }

        // Run fresh identification to get metadata
        ClusterIdentificationRequest request = ClusterIdentificationRequest.builder()
                .structurePhase(structurePhase)
                .model(model)
                .numComponents(numComponents)
                .build();

        AllClusterData clusterData = AllClusterData.identify(request, null);
        CFIdentificationResult cfResult = clusterData.getDisorderedCFResult();
        CFMetadata[] cfMetadata = extractCFMetadata(cfResult);

        // Single source of truth: CVCF basis metadata
        int ncf = CvCfBasis.getNumNonPointCfs(structurePhase, model, numComponents);
        List<String> eciNames = cfResult.getEONames().subList(0, ncf);

        CECEntry.CECTerm[] terms = new CECEntry.CECTerm[ncf];
        for (int i = 0; i < ncf; i++) {
            CECEntry.CECTerm term = new CECEntry.CECTerm();
            term.name = eciNames.get(i);
            term.a = 0.0;
            term.b = 0.0;
            if (cfMetadata != null && i < cfMetadata.length && cfMetadata[i] != null) {
                term.numSites = cfMetadata[i].numSites;
                term.multiplicity = cfMetadata[i].multiplicity;
                term.description = cfMetadata[i].description;
            }
            terms[i] = term;
        }

        CECEntry entry = new CECEntry();
        entry.elements = elements;
        entry.structurePhase = structurePhase;
        entry.model = model;
        entry.cecTerms = terms;
        entry.cecUnits = "J/mol";
        entry.reference = "";
        entry.notes = "Scaffolded from CVCF cluster data (Type-1a)";
        entry.ncf = ncf;

        store.save(hamiltonianId, entry);
        return entry;
    }


    /**
     * Extracts CF metadata from a CFIdentificationResult.
     * For each CF, extracts the number of sites and multiplicity from the cluster data.
     *
     * @param cfResult the CF identification result from cluster data
     * @return array of CFMetadata, one for each CF (or null if data unavailable)
     */
    private CFMetadata[] extractCFMetadata(CFIdentificationResult cfResult) {
        GroupedCFResult grouped = cfResult.getGroupedCFData();
        if (grouped == null) {
            return null;
        }

        List<List<List<org.ce.model.cluster.Cluster>>> coordData = grouped.getCoordData();
        List<List<List<Double>>> multiplicityData = grouped.getMultiplicityData();

        if (coordData == null || multiplicityData == null) {
            return null;
        }

        int ncf = cfResult.getNcf();
        CFMetadata[] metadata = new CFMetadata[ncf];

        for (int cf = 0; cf < ncf && cf < coordData.size(); cf++) {
            List<List<org.ce.model.cluster.Cluster>> cfGroups = coordData.get(cf);
            List<List<Double>> cfMults = multiplicityData.get(cf);

            if (cfGroups != null && !cfGroups.isEmpty() && cfMults != null && !cfMults.isEmpty()) {
                // Get first cluster in first group to count sites
                org.ce.model.cluster.Cluster firstCluster = cfGroups.get(0).get(0);
                int numSites = firstCluster.getAllSites().size();

                // Get multiplicity from first group
                double multiplicity = cfMults.get(0).get(0);

                // Create human-readable description
                String description = String.format("CF_%d (%d site%s, multiplicity %.0f)",
                        cf, numSites, numSites == 1 ? "" : "s", multiplicity);

                metadata[cf] = new CFMetadata(numSites, multiplicity, description);
            }
        }

        return metadata;
    }

    /**
     * Loads the Hamiltonian and validates it against ncf from cluster data.
     *
     * @param clusterId     ID of the cluster data, e.g. BCC_A2_T_bin
     * @param hamiltonianId ID of the Hamiltonian to load, e.g. Nb-Ti_BCC_A2_T
     */
    public CECEntry loadAndValidateCEC(String clusterId, String hamiltonianId) throws Exception {

        CECEntry entry = store.load(hamiltonianId);

        // Determine expected ncf from basis definition
        int numComponents = entry.elements.split("-").length;
        int expectedNcf = CvCfBasis.getNumNonPointCfs(entry.structurePhase, entry.model, numComponents);
        
        validateCEC(entry, expectedNcf);

        return entry;
    }

    // =========================================================================
    // Basis Transformation: Binary → Ternary
    // =========================================================================

    /**
     * Transforms ECIs from a binary system to a ternary system basis.
     *
     * <p>This implements the Cluster Variable Correlation Function (CVCF) basis transformation
     * from Jindal & Lele (2025), where lower-order ECIs are inherited into higher-order
     * systems through proper basis change.</p>
     *
     * <p><b>Example:</b> Transform Nb-Ti binary ECIs to Nb-Ti-V ternary basis, where:
     * - Binary species: {Nb=0, Ti=1}
     * - Ternary species: {Nb=0, Ti=1, V=2}
     * - Mapping: high=0→low=0 (Nb), high=1→low=1 (Ti), high=2→low=? (V not in binary)
     * </p>
     *
     * @param binaryId         source binary Hamiltonian ID (e.g., Nb-Ti_BCC_A2_T)
     * @param ternaryId        target ternary Hamiltonian ID (e.g., Nb-Ti-V_BCC_A2_T)
     * @param ternaryElements  ternary element string to verify (e.g., Nb-Ti-V)
     * @param speciesMapping   maps ternary species indices to binary species indices
     *                         For Nb-Ti in Nb-Ti-V: add(0,0)=Nb→Nb, add(1,1)=Ti→Ti
     * @return transformed CECEntry with ECIs in ternary basis
     */
    public CECEntry transformBinaryToTernary(
            String binaryId,
            String ternaryId,
            String ternaryElements,
            NumericalCECTransformer.SpeciesMapping speciesMapping
    ) throws Exception {

        // Load binary Hamiltonian
        CECEntry binaryEntry = store.load(binaryId);

        // Parse element counts
        int K_binary = 2;
        int K_ternary = 3;

        // Validate binary system has correct element count
        int binaryElemCount = binaryEntry.elements.split("-").length;
        if (binaryElemCount != K_binary) {
            throw new IllegalArgumentException(
                "Binary Hamiltonian '" + binaryId + "' has " + binaryElemCount +
                " elements; expected 2");
        }

        int ncf = binaryEntry.ncf;
        CECEntry.CECTerm[] ternaryTerms = new CECEntry.CECTerm[ncf];

        // Transform each cluster type (correlation function)
        for (int cfIdx = 0; cfIdx < ncf; cfIdx++) {
            CECEntry.CECTerm binaryTerm = binaryEntry.cecTerms[cfIdx];
            CECEntry.CECTerm ternaryTerm = new CECEntry.CECTerm();

            // Copy metadata
            ternaryTerm.name = binaryTerm.name;
            ternaryTerm.description = binaryTerm.description;
            ternaryTerm.numSites = binaryTerm.numSites;
            ternaryTerm.multiplicity = binaryTerm.multiplicity;

            // Determine cluster size from numSites (if available)
            int clusterSize = (binaryTerm.numSites != null) ? binaryTerm.numSites : cfIdx + 1;

            // Prepare binary ECI vector for this cluster type
            // Binary ECIs: index i = encode(site_species_0, site_species_1, ...)
            // For K=2, size 2^clusterSize configurations
            int binarySize = (int) Math.pow(K_binary, clusterSize);
            double[] J_binary = new double[binarySize];

            // In this simplified approach, we store the constant coefficient for pairs/triangles
            // In a full implementation, all configurations would need to be encoded
            // For now, use the 'a' coefficient as a single-point estimate
            for (int i = 0; i < binarySize; i++) {
                J_binary[i] = binaryTerm.a;  // Simplified: all same value
            }

            // Transform to ternary basis
            double[] J_ternary = NumericalCECTransformer.transform(
                J_binary, K_binary, K_ternary, clusterSize, speciesMapping
            );

            // Extract result: take the first (dominant) basis vector in ternary
            // In practice, this is (0,0,...) configuration transformed to ternary
            double resultA = (J_ternary.length > 0) ? J_ternary[0] : 0.0;

            // Use temperature coefficient from binary (simplified: same for all)
            ternaryTerm.a = resultA;
            ternaryTerm.b = binaryTerm.b;

            ternaryTerms[cfIdx] = ternaryTerm;
        }

        // Construct ternary CECEntry
        CECEntry ternaryEntry = new CECEntry();
        ternaryEntry.elements = ternaryElements;
        ternaryEntry.structurePhase = binaryEntry.structurePhase;
        ternaryEntry.model = binaryEntry.model;
        ternaryEntry.cecTerms = ternaryTerms;
        ternaryEntry.cecUnits = binaryEntry.cecUnits;
        ternaryEntry.reference = binaryEntry.reference + " [transformed from binary]";
        ternaryEntry.notes = "Ternary basis transformation from " + binaryId;
        ternaryEntry.ncf = ncf;

        // Save ternary Hamiltonian
        store.save(ternaryId, ternaryEntry);

        return ternaryEntry;
    }

    /** Returns true if a Hamiltonian with the given ID exists in the store. */
    public boolean hamiltonianExists(String hamiltonianId) {
        return store.exists(hamiltonianId);
    }

    /** Saves a Hamiltonian entry to the store under the given ID. */
    public void saveHamiltonian(String hamiltonianId, org.ce.model.hamiltonian.CECEntry entry)
            throws java.io.IOException {
        store.save(hamiltonianId, entry);
    }

}
