package org.ce.model.storage;

import org.ce.model.hamiltonian.CECEntry;

/**
 * Storage class responsible for loading and saving Cluster Expansion Coefficient (CEC) databases.
 */
public class HamiltonianStore extends DataStore<CECEntry> {

    public HamiltonianStore(Workspace workspace) {
        super(workspace::hamiltonianFile, CECEntry.class);
    }

}
