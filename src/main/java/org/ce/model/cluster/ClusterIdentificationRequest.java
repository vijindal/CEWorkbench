package org.ce.model.cluster;

import static org.ce.model.cluster.ClusterPrimitives.*;
import org.ce.model.storage.InputLoader;
import org.ce.model.storage.Workspace;

import java.io.File;
import java.nio.file.Path;

/**
 * Configuration request for cluster and correlation function identification.
 *
 * <p>This class encapsulates all parameters needed to run the identification
 * workflow, including cluster files, symmetry groups, transformation matrices,
 * and system parameters.</p>
 */
public class ClusterIdentificationRequest {

    private final String disorderedClusterFile;
    private final String orderedClusterFile;
    private final String disorderedSymmetryGroup;
    private final String orderedSymmetryGroup;
    private final double[][] transformationMatrix;
    private final Vector3D translationVector;
    private final int numComponents;
    private final String structurePhase;
    private final String model;

    private ClusterIdentificationRequest(Builder builder) {
        this.disorderedClusterFile = builder.disorderedClusterFile;
        this.orderedClusterFile = builder.orderedClusterFile;
        this.disorderedSymmetryGroup = builder.disorderedSymmetryGroup;
        this.orderedSymmetryGroup = builder.orderedSymmetryGroup;
        this.transformationMatrix = builder.transformationMatrix;
        this.translationVector = builder.translationVector;
        this.numComponents = builder.numComponents;
        this.structurePhase = builder.structurePhase;
        this.model = builder.model;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    public String getDisorderedClusterFile() {
        return disorderedClusterFile;
    }

    public String getOrderedClusterFile() {
        return orderedClusterFile;
    }

    public String getDisorderedSymmetryGroup() {
        return disorderedSymmetryGroup;
    }

    public String getOrderedSymmetryGroup() {
        return orderedSymmetryGroup;
    }

    public double[][] getTransformationMatrix() {
        return transformationMatrix;
    }

    public Vector3D getTranslationVector() {
        return translationVector;
    }

    public int getNumComponents() {
        return numComponents;
    }

    public String getStructurePhase() {
        return structurePhase;
    }

    public String getModel() {
        return model;
    }

    // =========================================================================
    // Builder
    // =========================================================================

    public static Builder builder() {
        return new Builder();
    }

    /**
     * Creates a standard identification request from a high-level system identity.
     * Automatically derives component count and standard file paths.
     */
    public static ClusterIdentificationRequest fromSystem(String elements, String structure, String model) {
        int ncomp = elements.split("-").length;
        return builder()
                .numComponents(ncomp)
                .structurePhase(structure)
                .model(model)
                .build();
    }

    /**
     * Builder for constructing ClusterIdentificationRequest with fluent API.
     */
    public static class Builder {
        private String disorderedClusterFile;
        private String orderedClusterFile;
        private String disorderedSymmetryGroup;
        private String orderedSymmetryGroup;
        private double[][] transformationMatrix;
        private Vector3D translationVector;
        private int numComponents;
        private String structurePhase = "BCC_A2";
        private String model = "T";

        public Builder disorderedClusterFile(String file) {
            this.disorderedClusterFile = file;
            return this;
        }

        public Builder orderedClusterFile(String file) {
            this.orderedClusterFile = file;
            return this;
        }

        public Builder disorderedSymmetryGroup(String group) {
            this.disorderedSymmetryGroup = group;
            return this;
        }

        public Builder orderedSymmetryGroup(String group) {
            this.orderedSymmetryGroup = group;
            return this;
        }

        public Builder transformationMatrix(double[][] matrix) {
            this.transformationMatrix = matrix;
            return this;
        }

        public Builder translationVector(Vector3D vector) {
            this.translationVector = vector;
            return this;
        }

        public Builder numComponents(int numComp) {
            this.numComponents = numComp;
            return this;
        }

        public Builder structurePhase(String sp) {
            this.structurePhase = sp;
            return this;
        }

        public Builder model(String m) {
            this.model = m;
            return this;
        }

        public ClusterIdentificationRequest build() {
            // If the caller didn't set orderedClusterFile directly, derive it
            // from structurePhase+model (the CLI/CECManagementWorkflow path).
            // Either way, orderedClusterFile is now the single source of truth
            // for the ordered structure/model — deriveDisorderedParent() below
            // re-parses it fresh, so orderedSymmetryGroup is never computed
            // from a stale/default structurePhase value.
            if (orderedClusterFile == null || orderedClusterFile.isBlank()) {
                if (structurePhase == null || model == null) {
                    throw new IllegalArgumentException(
                            "Must set either orderedClusterFile, or both structurePhase and model.");
                }
                String base = structurePhase.replace("_CVCF", "");
                String mod = model.replace("_CVCF", "");
                this.orderedClusterFile = "clus/" + base + "-" + mod + ".txt";
            }

            // Derive the ordered symmetry group, disordered parent (structure,
            // model, cluster file, symmetry group) from orderedClusterFile via
            // StructurePhaseRegistry — the single source of truth for which
            // disordered structure an ordered structure belongs to. This is
            // the ONLY place this derivation should happen; GUI and CLI both
            // go through this Builder so they cannot diverge.
            //
            // Callers may still set disorderedClusterFile/disorderedSymmetryGroup
            // explicitly (e.g. tests, advanced use) to bypass derivation.
            if (disorderedClusterFile == null || disorderedClusterFile.isBlank()
                    || disorderedSymmetryGroup == null || disorderedSymmetryGroup.isBlank()) {
                deriveDisorderedParent();
            }

            validate();
            // Auto-extract transformation matrix and translation vector from symmetry group files if not set
            if (transformationMatrix == null || translationVector == null) {
                extractTransformationFromSymmetryGroup();
            }
            return new ClusterIdentificationRequest(this);
        }

        /**
         * Derives the disordered-parent cluster file and symmetry group from
         * {@link #orderedClusterFile}, using {@link StructurePhaseRegistry} for
         * the structure mapping. The disordered parent MUST use the same
         * maximal-cluster model letter as the ordered phase (e.g. both "T" or
         * both "TO") — the model selects which maximal clusters are used for
         * identification, so substituting a different model would build
         * cluster/CF data from the wrong maximal clusters entirely. This
         * method never falls back to a different model; it fails clearly
         * instead.
         */
        private void deriveDisorderedParent() {
            if (orderedClusterFile == null || orderedClusterFile.isBlank()) {
                throw new IllegalArgumentException(
                        "Cannot derive disordered parent: orderedClusterFile is not set. "
                                + "Provide structurePhase+model, or set orderedClusterFile explicitly.");
            }

            String name = orderedClusterFile;
            int slash = name.lastIndexOf('/');
            if (slash >= 0)
                name = name.substring(slash + 1);
            if (name.toLowerCase().endsWith(".txt"))
                name = name.substring(0, name.length() - 4);

            int dash = name.lastIndexOf('-');
            if (dash < 0) {
                throw new IllegalArgumentException(
                        "Cannot parse ordered cluster file '" + orderedClusterFile
                                + "' as <structure>-<model>; expected e.g. BCC_B2-T.txt");
            }
            String orderedStructure = name.substring(0, dash);
            String orderedModel = name.substring(dash + 1);

            // Keep structurePhase/model consistent with the parsed ordered
            // file regardless of whether the caller set them explicitly or
            // only supplied orderedClusterFile.
            this.structurePhase = orderedStructure;
            this.model = orderedModel;

            String disorderedStructure = StructurePhaseRegistry.parentOf(orderedStructure);

            Path inputsDir = new Workspace().inputsDir();
            File disorderedClusterHandle = inputsDir.resolve("clus")
                    .resolve(disorderedStructure + "-" + orderedModel + ".txt").toFile();
            if (!disorderedClusterHandle.exists()) {
                throw new ClusterIdentificationException("Stage 0",
                        "Disordered parent '" + disorderedStructure + "' has no '-" + orderedModel
                                + "' maximal-cluster model file (looked for " + disorderedClusterHandle
                                + "). The disordered parent must use the same cluster model as ordered structure '"
                                + orderedStructure + "-" + orderedModel
                                + "'. Add that file, or pick an ordered cluster file whose model the parent "
                                + "structure already has.");
            }

            this.disorderedClusterFile = "clus/" + disorderedStructure + "-" + orderedModel + ".txt";
            this.disorderedSymmetryGroup = disorderedStructure + "-SG";

            if (this.orderedSymmetryGroup == null || this.orderedSymmetryGroup.isBlank()) {
                this.orderedSymmetryGroup = orderedStructure + "-SG";
            }

            File orderedSymHandle = inputsDir.resolve("sym").resolve(orderedSymmetryGroup + ".txt").toFile();
            if (!orderedSymHandle.exists()) {
                throw new ClusterIdentificationException("Stage 0",
                        "Missing symmetry file for ordered structure: " + orderedSymHandle);
            }
            File disorderedSymHandle = inputsDir.resolve("sym").resolve(disorderedSymmetryGroup + ".txt").toFile();
            if (!disorderedSymHandle.exists()) {
                throw new ClusterIdentificationException("Stage 0",
                        "Missing symmetry file for disordered parent: " + disorderedSymHandle);
            }
        }

        private void validate() {
            if (disorderedClusterFile == null || disorderedClusterFile.isBlank()) {
                throw new IllegalArgumentException("disorderedClusterFile must not be blank");
            }
            if (orderedClusterFile == null || orderedClusterFile.isBlank()) {
                throw new IllegalArgumentException("orderedClusterFile must not be blank");
            }
            if (disorderedSymmetryGroup == null || disorderedSymmetryGroup.isBlank()) {
                throw new IllegalArgumentException("disorderedSymmetryGroup must not be blank");
            }
            if (orderedSymmetryGroup == null || orderedSymmetryGroup.isBlank()) {
                throw new IllegalArgumentException("orderedSymmetryGroup must not be blank");
            }
            if (numComponents < 2) {
                throw new IllegalArgumentException("numComponents must be >= 2");
            }
        }

        private void extractTransformationFromSymmetryGroup() {
            try {
                // Load the disordered symmetry group to extract transformation matrix and vector
                SpaceGroup disorderedSG = InputLoader.parseSpaceGroup(disorderedSymmetryGroup);
                this.transformationMatrix = disorderedSG.getRotateMat();
                double[] translateMat = disorderedSG.getTranslateMat();
                this.translationVector = new Vector3D(translateMat[0], translateMat[1], translateMat[2]);
            } catch (Exception e) {
                throw new RuntimeException(
                    "Failed to extract transformation matrix and vector from symmetry group '" +
                    disorderedSymmetryGroup + "': " + e.getMessage(), e);
            }
        }
    }
}
