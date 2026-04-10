package org.ce.calculation.workflow;

import static org.ce.model.cluster.ClusterPrimitives.*;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.storage.InputLoader;

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
    private final String engineType;

    private static String sanitize(String id) {
        if (id == null) return null;
        return id.replace("_CVCF", "");
    }

    /**
     * Constructor for CVM calculation.
     */
    public ClusterIdentificationRequest(org.ce.calculation.engine.ThermodynamicEngine.Input input) {
        this(input, "CVM");
    }

    /**
     * Constructor for specific engine types (e.g., MCS).
     */
    public ClusterIdentificationRequest(org.ce.calculation.engine.ThermodynamicEngine.Input input, String engineType) {
        org.ce.model.hamiltonian.CECEntry cec = input.cec;
        this.structurePhase = cec.structurePhase;
        this.model = sanitize(cec.model);
        this.numComponents = input.composition.length;
        this.engineType = engineType;

        String base = sanitize(this.structurePhase);
        String mod = sanitize(this.model);
        
        this.disorderedClusterFile = "clus/" + base + "-" + mod + ".txt";
        this.orderedClusterFile = this.disorderedClusterFile;
        this.disorderedSymmetryGroup = base + "-SG";
        this.orderedSymmetryGroup = this.disorderedSymmetryGroup;

        // Extract transformation from symmetry group immediately
        try {
            org.ce.model.cluster.SpaceGroup sg = org.ce.model.storage.InputLoader.parseSpaceGroup(this.disorderedSymmetryGroup);
            this.transformationMatrix = sg.getRotateMat();
            double[] translateMat = sg.getTranslateMat();
            this.translationVector = new org.ce.model.cluster.ClusterPrimitives.Vector3D(translateMat[0], translateMat[1], translateMat[2]);
        } catch (Exception e) {
            throw new RuntimeException("Failed to extract transformation from symmetry group: " + this.disorderedSymmetryGroup, e);
        }
    }

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
        this.engineType = builder.engineType;
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

    public String getEngineType() {
        return engineType;
    }

    // =========================================================================
    // Builder
    // =========================================================================

    public static Builder builder() {
        return new Builder();
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
        private String engineType = "CVM";

        public Builder engineType(String type) {
            this.engineType = type;
            return this;
        }

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
            validate();
            // Auto-extract transformation matrix and translation vector from symmetry group files if not set
            if (transformationMatrix == null || translationVector == null) {
                extractTransformationFromSymmetryGroup();
            }
            return new ClusterIdentificationRequest(this);
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
