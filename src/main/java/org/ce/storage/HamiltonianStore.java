package org.ce.storage;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.ce.domain.hamiltonian.CECEntry;

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

    private Path getHamiltonianFile(String hamiltonianId) {
        return workspace.hamiltonianFile(hamiltonianId);
    }

    public boolean exists(String hamiltonianId) {
        return Files.exists(getHamiltonianFile(hamiltonianId));
    }

    /**
     * Loads the Hamiltonian (ECI parameters) for the given ID.
     * ID convention: {elements}_{structure}_{model}, e.g. Nb-Ti_BCC_A2_T
     */
    public CECEntry load(String hamiltonianId) throws IOException {

        Path file = getHamiltonianFile(hamiltonianId);
        System.out.println("  [STORAGE] Loading Hamiltonian from: " + file.toAbsolutePath());

        if (!Files.exists(file)) {
            throw new IOException("Hamiltonian not found: " + file);
        }

        return mapper.readValue(file.toFile(), CECEntry.class);
    }

    /**
     * Saves the Hamiltonian (ECI parameters) for the given ID.
     * ID convention: {elements}_{structure}_{model}, e.g. Nb-Ti_BCC_A2_T
     */
    public void save(String hamiltonianId, CECEntry entry) throws IOException {

        Path file = getHamiltonianFile(hamiltonianId);

        Files.createDirectories(file.getParent());

        mapper.writerWithDefaultPrettyPrinter()
              .writeValue(file.toFile(), entry);
    }

}
