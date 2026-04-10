package org.ce.model.cluster;

/**
 * Domain-level immutable input contract for CVM engine/model execution.
 *
 * <p>This bundles the precomputed Stage 1-3 topology data plus system metadata,
 * allowing CVM classes to remain decoupled from legacy workbench containers.</p>
 */
public final class CVMModelInput {

    private final String systemId;
    private final String systemName;
    private final int numComponents;

    private final ClusterIdentificationResult stage1;
    private final CFIdentificationResult stage2;
    private final CMatrix.Result stage3;

    public CVMModelInput(
            String systemId,
            String systemName,
            int numComponents,
            ClusterIdentificationResult stage1,
            CFIdentificationResult stage2,
            CMatrix.Result stage3) {

        if (systemId == null || systemId.isBlank()) {
            throw new IllegalArgumentException("systemId must not be blank");
        }
        if (systemName == null || systemName.isBlank()) {
            throw new IllegalArgumentException("systemName must not be blank");
        }
        if (numComponents < 2) {
            throw new IllegalArgumentException("numComponents must be >= 2");
        }
        if (stage1 == null || stage2 == null || stage3 == null) {
            throw new IllegalArgumentException("Stage 1/2/3 data must be non-null");
        }
        if (stage3.getCmat() == null || stage3.getLcv() == null || stage3.getWcv() == null
                || stage3.getCfBasisIndices() == null) {
            throw new IllegalArgumentException("Stage 3 is missing required fields (cmat/lcv/wcv/cfBasisIndices)");
        }
        if (stage3.getCfBasisIndices().length != stage2.getTcf()) {
            throw new IllegalArgumentException("Stage 3 cfBasisIndices length mismatch: got "
                    + stage3.getCfBasisIndices().length + ", expected " + stage2.getTcf());
        }

        this.systemId = systemId;
        this.systemName = systemName;
        this.numComponents = numComponents;
        this.stage1 = stage1;
        this.stage2 = stage2;
        this.stage3 = stage3;
    }

    public String getSystemId() {
        return systemId;
    }

    public String getSystemName() {
        return systemName;
    }

    public int getNumComponents() {
        return numComponents;
    }

    public ClusterIdentificationResult getStage1() {
        return stage1;
    }

    public CFIdentificationResult getStage2() {
        return stage2;
    }

    public CMatrix.Result getStage3() {
        return stage3;
    }
}


