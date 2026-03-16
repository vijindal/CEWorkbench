package org.ce.workflow.cec;

import org.ce.domain.cluster.AllClusterData;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.cec.CECEntry;
import org.ce.storage.cec.CECTerm;
import org.ce.storage.cec.HamiltonianStore;

import java.io.IOException;

/**
 * Workflow responsible for managing Cluster Expansion Coefficient (CEC) databases.
 *
 * This orchestrates operations between cluster data and CEC storage.
 */
public class CECManagementWorkflow {

    private final HamiltonianStore store;

    private final ClusterDataStore clusterStore;

    public CECManagementWorkflow(
            HamiltonianStore store,
            ClusterDataStore clusterStore) {

        this.store = store;
        this.clusterStore = clusterStore;
    }

    /**
     * Creates a blank CEC database for a system.
     *
     * @param systemId system identifier
     * @param elements element string (e.g. "Nb-Ti")
     * @param structurePhase structure identifier
     * @param model model identifier
     * @param ncf number of correlation functions
     */
    public CECEntry scaffoldEmptyCEC(
            String systemId,
            String elements,
            String structurePhase,
            String model,
            int ncf) {

        CECTerm[] terms = new CECTerm[ncf];

        for (int i = 0; i < ncf; i++) {
            CECTerm term = new CECTerm();
            term.name = "CF_" + i;
            term.a = 0.0;
            term.b = 0.0;
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
     */
    public CECEntry createAndSaveCEC(
            String systemId,
            String elements,
            String structurePhase,
            String model,
            int ncf) throws IOException {

        CECEntry entry = scaffoldEmptyCEC(
                systemId,
                elements,
                structurePhase,
                model,
                ncf
        );

        store.save(systemId, entry);

        return entry;
    }

    /**
     * Validates that the CEC database matches the expected number
     * of correlation functions.
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
    }

    /**
     * Creates and saves a CEC database automatically using
     * cluster identification results.
     */
    public CECEntry scaffoldFromClusterData(
            String systemId,
            String elements,
            String structurePhase,
            String model) throws Exception {

        AllClusterData clusterData = clusterStore.load(systemId);

        int ncf = clusterData.getDisorderedCFResult().getNcf();

        return createAndSaveCEC(
                systemId,
                elements,
                structurePhase,
                model,
                ncf
        );
    }

    /**
     * Loads the CEC database and validates it against
     * the number of correlation functions in cluster data.
     */
    public CECEntry loadAndValidateCEC(String systemId) throws Exception {

        AllClusterData clusterData = clusterStore.load(systemId);

        int ncf = clusterData.getDisorderedCFResult().getNcf();

        CECEntry entry = store.load(systemId);

        validateCEC(entry, ncf);

        return entry;
    }

}
