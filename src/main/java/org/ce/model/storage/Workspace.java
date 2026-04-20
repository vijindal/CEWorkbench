package org.ce.model.storage;

import java.nio.file.Path;

/**
 * Defines all filesystem locations used by CE Workbench.
 *
 * <p>Default root: {@code ~/CEWorkbench/}. Pass a different root to the
 * constructor to override (e.g., for tests or a user-selected project folder).</p>
 *
 * <h2>Layout and ID conventions</h2>
 * <pre>
 *   ~/CEWorkbench/
 *   │
 *   ├─ inputs/                          ← Type-1 INPUT files (cluster geometry, symmetry)
 *   │   ├─ clus/
 *   │   │   ├─ BCC_A2-T.txt
 *   │   │   └─ BCC_B2-T.txt
 *   │   └─ sym/
 *   │       ├─ BCC_A2-SG.txt
 *   │       └─ BCC_B2-SG.txt
 *   │
 *   └─ hamiltonians/                    ← ECI parameters  (element-specific)
 *         ID = {elements}_{structure}_{model}
 *         e.g. Nb-Ti_BCC_A2_T, Al-Nb_BCC_A2_T
 *       └─ Nb-Ti_BCC_A2_T/
 *           └─ hamiltonian.json
 * </pre>
 */
public class Workspace {

    /** Default workspace directory: ~/CEWorkbench */
    public static final Path DEFAULT_ROOT =
            Path.of(System.getProperty("user.home"), "CEWorkbench");

    /** Project-local workspace (development): ./data/CEWorkbench */
    public static final Path PROJECT_LOCAL_ROOT = Path.of("data", "CEWorkbench");

    private final Path root;

    /** Creates a workspace at the given root directory. */
    public Workspace(Path root) {
        this.root = root;
    }

    /** Creates a workspace at the project-local location if it exists, else uses default. */
    public Workspace() {
        this(java.nio.file.Files.exists(PROJECT_LOCAL_ROOT) ? PROJECT_LOCAL_ROOT : DEFAULT_ROOT);
    }

    /** Root workspace directory. */
    public Path getRoot() {
        return root;
    }

    // -------------------------------------------------------------------------
    // Type-1 inputs  (cluster geometry + symmetry files)
    // -------------------------------------------------------------------------

    /**
     * Base directory for Type-1 input files: {@code <root>/inputs/}.
     * Mirrors the classpath layout — {@code inputs/clus/} and {@code inputs/sym/}.
     */
    public Path inputsDir() {
        return root.resolve("inputs");
    }

    // -------------------------------------------------------------------------
    // Hamiltonians  (element-specific ECI parameters)
    // ID convention: {elements}_{structure}_{model}  e.g. Nb-Ti_BCC_A2_T
    // -------------------------------------------------------------------------

    /** Directory for a Hamiltonian: {@code <root>/hamiltonians/<hamiltonianId>/}. */
    public Path hamiltonianDir(String hamiltonianId) {
        return root.resolve("hamiltonians").resolve(hamiltonianId);
    }

    /** Path to the Hamiltonian JSON: {@code <root>/hamiltonians/<hamiltonianId>/hamiltonian.json}. */
    public Path hamiltonianFile(String hamiltonianId) {
        return hamiltonianDir(hamiltonianId).resolve("hamiltonian.json");
    }

    // =========================================================================
    // System Identity & ID Derivation
    // =========================================================================

    /**
     * Encapsulates the three-part system identity (elements, structure, model)
     * and derives the storage IDs used by this Workspace.
     */
    public record SystemId(String elements, String structure, String model) {

        /**
         * Cluster data ID (element-agnostic): {@code {structure}_{model}_{ncomp}}
         * e.g. {@code BCC_A2_T_bin}
         */
        public String clusterId() {
            String modelPart = model.trim();
            if (modelPart.equalsIgnoreCase("CVCF") || modelPart.toUpperCase().endsWith("_CVCF")) {
                return structure + "_CVCF_" + ncompSuffix();
            }
            return structure + "_" + model + "_" + ncompSuffix();
        }

        /**
         * Hamiltonian ID (element-specific): {@code {elements}_{structure}_{model}}
         * e.g. {@code Nb-Ti_BCC_A2_T}
         */
        public String hamiltonianId() {
            return elements + "_" + structure + "_" + model;
        }

        /** Returns true only if all three fields are non-blank. */
        public boolean isComplete() {
            return elements != null && !elements.isBlank() &&
                   structure != null && !structure.isBlank() &&
                   model != null && !model.isBlank();
        }

        private String ncompSuffix() {
            return ncompSuffix(elements.split("-").length);
        }

        /**
         * Returns the canonical cluster-ID suffix for the given component count.
         */
        public static String ncompSuffix(int ncomp) {
            return switch (ncomp) {
                case 2 -> "bin";
                case 3 -> "tern";
                case 4 -> "quat";
                default -> throw new IllegalArgumentException(
                        "No ncomp suffix for " + ncomp + " components");
            };
        }
    }
}
