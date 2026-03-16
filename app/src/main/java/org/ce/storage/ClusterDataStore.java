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

    private Path getClusterFile(String systemId) {
        return workspace.systemDir(systemId)
                .resolve("all_cluster_data.json");
    }

    public boolean exists(String systemId) {
        return Files.exists(getClusterFile(systemId));
    }

    public AllClusterData load(String systemId) throws IOException {

        Path file = getClusterFile(systemId);

        if (!Files.exists(file)) {
            throw new IOException("Cluster data not found: " + file);
        }

        return mapper.readValue(file.toFile(), AllClusterData.class);
    }

    public void save(String systemId, AllClusterData data) throws IOException {

        Path file = getClusterFile(systemId);

        Files.createDirectories(file.getParent());

        mapper.writerWithDefaultPrettyPrinter()
                .writeValue(file.toFile(), data);
    }
}
