package org.ce.storage;

import java.nio.file.Path;

/**
 * Defines filesystem locations used by the application.
 *
 * The workspace is the root directory where user data
 * (CEC database, results, etc.) are stored.
 */
public class Workspace {

    private final Path root;

    public Workspace(Path root) {
        this.root = root;
    }

    /**
     * Root workspace directory.
     */
    public Path getRoot() {
        return root;
    }

    /**
     * Directory containing all system databases.
     */
    public Path systemsDir() {
        return root.resolve("data").resolve("systems");
    }

    /**
     * Directory for a specific system.
     */
    public Path systemDir(String systemId) {
        return systemsDir().resolve(systemId);
    }

    /**
     * Path to the CEC database file.
     */
    public Path cecFile(String systemId) {
        return systemDir(systemId).resolve("cec.json");
    }
}
