package org.ce.workflow;

import org.ce.domain.cluster.Vector3D;

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

    private ClusterIdentificationRequest(Builder builder) {
        this.disorderedClusterFile = builder.disorderedClusterFile;
        this.orderedClusterFile = builder.orderedClusterFile;
        this.disorderedSymmetryGroup = builder.disorderedSymmetryGroup;
        this.orderedSymmetryGroup = builder.orderedSymmetryGroup;
        this.transformationMatrix = builder.transformationMatrix;
        this.translationVector = builder.translationVector;
        this.numComponents = builder.numComponents;
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

        public ClusterIdentificationRequest build() {
            validate();
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
            if (transformationMatrix == null || transformationMatrix.length != 3 || transformationMatrix[0].length != 3) {
                throw new IllegalArgumentException("transformationMatrix must be a 3x3 matrix");
            }
            if (translationVector == null) {
                throw new IllegalArgumentException("translationVector must not be null");
            }
            if (numComponents < 2) {
                throw new IllegalArgumentException("numComponents must be >= 2");
            }
        }
    }
}
