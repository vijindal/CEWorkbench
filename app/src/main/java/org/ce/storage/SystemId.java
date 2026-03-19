package org.ce.storage;

/**
 * Encapsulates the three-part system identity (elements, structure, model)
 * and derives the two storage IDs from them.
 *
 * <p>This is the single canonical source of ID-derivation logic, used by
 * all UI panels (GUI and CLI) and the workflow layer.</p>
 *
 * <ul>
 *   <li>{@code clusterId}     = {@code {structure}_{model}_{ncomp}}  e.g. {@code BCC_A2_T_bin}</li>
 *   <li>{@code hamiltonianId} = {@code {elements}_{structure}_{model}} e.g. {@code Nb-Ti_BCC_A2_T}</li>
 * </ul>
 */
public class SystemId {

    public final String elements;
    public final String structure;
    public final String model;

    public SystemId(String elements, String structure, String model) {
        this.elements  = elements;
        this.structure = structure;
        this.model     = model;
    }

    /**
     * Cluster data ID (element-agnostic): {@code {structure}_{model}_{ncomp}}
     * e.g. {@code BCC_A2_T_bin}
     */
    public String clusterId() {
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
        return !elements.isBlank() && !structure.isBlank() && !model.isBlank();
    }

    private String ncompSuffix() {
        int ncomp = elements.split("-").length;
        return switch (ncomp) {
            case 2 -> "bin";
            case 3 -> "tern";
            case 4 -> "quat";
            default -> throw new IllegalArgumentException(
                    "No ncomp suffix for " + ncomp + " components");
        };
    }
}
