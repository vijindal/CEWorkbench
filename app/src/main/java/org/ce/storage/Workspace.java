package org.ce.storage;

import java.nio.file.Path;

/**
 * Represents the CEWorkbench workspace.
 * Responsible for resolving storage locations.
 */
public class Workspace {

    private final Path root;

    public Workspace(Path root) {
        this.root = root;
    }

    public Path getRoot() {
        return root;
    }

    public Path getClusterDataDirectory() {
        return root.resolve("data").resolve("cluster");
    }

    public Path getHamiltonianDirectory() {
        return root.resolve("data").resolve("hamiltonian");
    }
}
