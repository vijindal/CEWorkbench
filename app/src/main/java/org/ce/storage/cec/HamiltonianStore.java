package org.ce.storage.cec;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.ce.storage.Workspace;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Storage class responsible for loading and saving
 * Cluster Expansion Coefficient (CEC) databases.
 */
public class HamiltonianStore {

    private final Workspace workspace;

    private final ObjectMapper mapper = new ObjectMapper();

    public HamiltonianStore(Workspace workspace) {
        this.workspace = workspace;
    }

    private Path getCECFile(String systemId) {
        return workspace.cecFile(systemId);
    }

    public boolean exists(String systemId) {
        return Files.exists(getCECFile(systemId));
    }

    /**
     * Loads the CEC database for the given system.
     */
    public CECEntry load(String systemId) throws IOException {

        Path file = getCECFile(systemId);

        if (!Files.exists(file)) {
            throw new IOException("CEC database not found: " + file);
        }

        return mapper.readValue(file.toFile(), CECEntry.class);
    }

    /**
     * Saves the CEC database for the given system.
     */
    public void save(String systemId, CECEntry entry) throws IOException {

        Path file = getCECFile(systemId);

        Files.createDirectories(file.getParent());

        mapper.writerWithDefaultPrettyPrinter()
              .writeValue(file.toFile(), entry);
    }

}
