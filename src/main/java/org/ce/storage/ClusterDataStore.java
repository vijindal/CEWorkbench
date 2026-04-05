package org.ce.storage;

import org.ce.domain.cluster.AllClusterData;

/**
 * Storage class responsible for loading and saving cluster identification results.
 */
public class ClusterDataStore extends DataStore<AllClusterData> {

    public ClusterDataStore(Workspace workspace) {
        super(workspace::clusterDataFile, AllClusterData.class);
    }

}
