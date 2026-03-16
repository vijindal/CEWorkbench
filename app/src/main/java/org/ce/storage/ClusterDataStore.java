package org.ce.storage;

import org.ce.domain.cluster.ClusterData;

/**
 * Storage service for cluster data.
 */
public class ClusterDataStore {

    private final Workspace workspace;

    public ClusterDataStore(Workspace workspace) {
        this.workspace = workspace;
    }

    public ClusterData load(String systemId) {

        // Placeholder implementation.
        // Will later load JSON cluster data.

        return null;
    }

}
