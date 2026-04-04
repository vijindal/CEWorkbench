package org.ce.storage;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.ce.domain.cluster.AllClusterData;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Storage class responsible for loading and saving
 * cluster identification results.
 */
public class ClusterDataStore {

    private final Workspace workspace;

    private final ObjectMapper mapper = new ObjectMapper();

    public ClusterDataStore(Workspace workspace) {
        this.workspace = workspace;
    }

    private Path getClusterFile(String clusterId) {
        return workspace.clusterDataFile(clusterId);
    }

    public boolean exists(String clusterId) {
        return Files.exists(getClusterFile(clusterId));
    }

    public AllClusterData load(String clusterId) throws IOException {

        Path file = getClusterFile(clusterId);

        if (!Files.exists(file)) {
            throw new IOException("Cluster data not found: " + file);
        }

        return mapper.readValue(file.toFile(), AllClusterData.class);
    }

    public void save(String clusterId, AllClusterData data) throws IOException {

        Path file = getClusterFile(clusterId);

        Files.createDirectories(file.getParent());

        mapper.writerWithDefaultPrettyPrinter()
                .writeValue(file.toFile(), data);
    }
}
