package org.ce.domain.engine.cvm;

import org.ce.domain.cluster.CFIdentificationResult;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.ClusterIdentificationResult;

/**
 * Domain-level immutable input contract for CVM engine/model execution.
 *
 * <p>This bundles the precomputed Stage 1-3 topology data plus system metadata,
 * allowing CVM classes to remain decoupled from legacy workbench containers.</p>
 */
public class CVMInput {

    private final ClusterIdentificationResult stage1;
    private final CFIdentificationResult stage2;
    private final CMatrixResult stage3;

    private final String systemId;
    private final String systemName;
    private final int numComponents;

    public CVMInput(
            ClusterIdentificationResult stage1,
            CFIdentificationResult stage2,
            CMatrixResult stage3,
            String systemId,
            String systemName,
            int numComponents) {

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

        this.stage1 = stage1;
        this.stage2 = stage2;
        this.stage3 = stage3;
        this.systemId = systemId;
        this.systemName = systemName;
        this.numComponents = numComponents;
    }

    public ClusterIdentificationResult getStage1() {
        return stage1;
    }

    public CFIdentificationResult getStage2() {
        return stage2;
    }

    public CMatrixResult getStage3() {
        return stage3;
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
}
