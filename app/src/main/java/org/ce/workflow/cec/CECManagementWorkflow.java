package org.ce.workflow.cec;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.GroupedCFResult;
import org.ce.storage.ClusterDataStore;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.storage.HamiltonianStore;

import java.io.IOException;
import java.util.List;

/**
 * Workflow responsible for managing Cluster Expansion Coefficient (CEC) databases.
 *
 * This orchestrates operations between cluster data and CEC storage.
 */
public class CECManagementWorkflow {

    public final HamiltonianStore store;

    private final ClusterDataStore clusterStore;

    public CECManagementWorkflow(
            HamiltonianStore store,
            ClusterDataStore clusterStore) {

        this.store = store;
        this.clusterStore = clusterStore;
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

        CECTerm[] terms = new CECTerm[ncf];

        for (int i = 0; i < ncf; i++) {
            CECTerm term = new CECTerm();
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
     * @param cfMetadata    optional CF metadata (numSites, multiplicity) for enriching CECTerms
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

        if (entry.cecTerms == null) {
            throw new IllegalStateException("CEC database contains no terms");
        }

        if (entry.cecTerms.length != expectedNcf) {
            throw new IllegalStateException(
                "CEC term count (" + entry.cecTerms.length +
                ") does not match expected number of correlation functions (" +
                expectedNcf + ")"
            );
        }

        if (entry.ncf != expectedNcf) {
            throw new IllegalStateException(
                "CEC ncf field (" + entry.ncf +
                ") does not match expected ncf (" + expectedNcf + ")"
            );
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
     * the scaffolded CECTerms with helpful descriptive information.</p>
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

        String clusterId = deriveClusterId(hamiltonianId);

        if (!clusterStore.exists(clusterId)) {
            throw new IllegalStateException(
                "Cluster data '" + clusterId + "' not found. Run type1a first to generate cluster_data.json.");
        }

        AllClusterData clusterData = clusterStore.load(clusterId);
        CFIdentificationResult cfResult = clusterData.getDisorderedCFResult();
        int ncf = cfResult.getNcf();

        // Extract CF metadata to enrich CECTerms
        CFMetadata[] cfMetadata = extractCFMetadata(cfResult);

        return createAndSaveCEC(hamiltonianId, elements, structurePhase, model, ncf, cfMetadata);
    }

    /**
     * Derives the cluster data ID from a Hamiltonian ID.
     *
     * <p>Format: {@code {elements}_{structure}_{model}} → {@code {structure}_{model}_{ncomp}}<br>
     * Example: {@code "A-B_BCC_A2_T"} → {@code "BCC_A2_T_bin"}</p>
     */
    private static String deriveClusterId(String hamiltonianId) {
        int sep = hamiltonianId.indexOf('_');
        if (sep < 0) {
            throw new IllegalArgumentException("Invalid hamiltonianId format: " + hamiltonianId);
        }
        String elements      = hamiltonianId.substring(0, sep);
        String structureModel = hamiltonianId.substring(sep + 1);
        int ncomp = elements.split("-").length;
        String suffix = switch (ncomp) {
            case 2 -> "bin";
            case 3 -> "tern";
            case 4 -> "quat";
            default -> throw new IllegalArgumentException(
                    "No ncomp suffix for " + ncomp + " components (hamiltonianId=" + hamiltonianId + ")");
        };
        return structureModel + "_" + suffix;
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

        List<List<List<org.ce.domain.cluster.Cluster>>> coordData = grouped.getCoordData();
        List<List<List<Double>>> multiplicityData = grouped.getMultiplicityData();

        if (coordData == null || multiplicityData == null) {
            return null;
        }

        int ncf = cfResult.getNcf();
        CFMetadata[] metadata = new CFMetadata[ncf];

        for (int cf = 0; cf < ncf && cf < coordData.size(); cf++) {
            List<List<org.ce.domain.cluster.Cluster>> cfGroups = coordData.get(cf);
            List<List<Double>> cfMults = multiplicityData.get(cf);

            if (cfGroups != null && !cfGroups.isEmpty() && cfMults != null && !cfMults.isEmpty()) {
                // Get first cluster in first group to count sites
                org.ce.domain.cluster.Cluster firstCluster = cfGroups.get(0).get(0);
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

        AllClusterData clusterData = clusterStore.load(clusterId);

        int ncf = clusterData.getDisorderedCFResult().getNcf();

        CECEntry entry = store.load(hamiltonianId);

        validateCEC(entry, ncf);

        return entry;
    }

}
