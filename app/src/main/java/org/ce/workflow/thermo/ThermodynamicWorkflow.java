package org.ce.workflow.thermo;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.ClusCoordListResult;
import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.ClusterIdentificationResult;
import org.ce.domain.cluster.cvcf.BccA2TModelCvCfTransformations;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.domain.engine.ThermodynamicEngine;
import org.ce.domain.engine.ThermodynamicInput;
import org.ce.domain.result.EquilibriumState;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicRequest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Workflow responsible for thermodynamic calculations.
 *
 * It loads the required scientific data and delegates
 * the computation to thermodynamic engines.
 */
public class ThermodynamicWorkflow {

    private static final Logger LOG = Logger.getLogger(ThermodynamicWorkflow.class.getName());

    private final ClusterDataStore clusterStore;

    private final CECManagementWorkflow cecWorkflow;

    private final ThermodynamicEngine cvmEngine;

    private final ThermodynamicEngine mcsEngine;

    public ThermodynamicWorkflow(
            ClusterDataStore clusterStore,
            CECManagementWorkflow cecWorkflow,
            ThermodynamicEngine cvmEngine,
            ThermodynamicEngine mcsEngine) {

        this.clusterStore = clusterStore;
        this.cecWorkflow = cecWorkflow;
        this.cvmEngine = cvmEngine;
        this.mcsEngine = mcsEngine;
    }

    /**
     * Loads the scientific data required for thermodynamic calculations.
     *
     * @param clusterId     element-agnostic cluster dataset ID, e.g. BCC_A2_T_bin
     * @param hamiltonianId element-specific Hamiltonian ID, e.g. Nb-Ti_BCC_A2_T
     */
    public ThermodynamicData loadThermodynamicData(
            String clusterId, String hamiltonianId) throws Exception {
        return loadThermodynamicData(clusterId, hamiltonianId, null);
    }

    public ThermodynamicData loadThermodynamicData(
            String clusterId, String hamiltonianId, Consumer<String> progressSink) throws Exception {

        LOG.fine("STAGE 3a: Load cluster data (Stages 1-3 topology)...");
        emit(progressSink, "STAGE 3a: Load cluster data (Stages 1-3 topology)");
        AllClusterData clusterData = clusterStore.load(clusterId);
        emit(progressSink, "  Loaded cluster data: " + clusterId);
        LOG.fine("  ✓ Loaded " + clusterId);

        LOG.fine("STAGE 3b: Load Hamiltonian (ECI coefficients)...");
        emit(progressSink, "STAGE 3b: Load Hamiltonian (ECI coefficients)");
        CECEntry cec = cecWorkflow.loadAndValidateCEC(clusterId, hamiltonianId);
        emit(progressSink, "  Loaded Hamiltonian: " + hamiltonianId);
        LOG.fine("  ✓ Loaded " + hamiltonianId);

        // Extract and log CEC data
        if (cec != null) {
            LOG.fine("    elements: " + cec.elements);
            LOG.fine("    structurePhase: " + cec.structurePhase);
            emit(progressSink, "  elements: " + cec.elements + " | structure: " + cec.structurePhase);
            if (cec.cecTerms != null && cec.cecTerms.length > 0) {
                LOG.fine("    ECI terms: " + cec.cecTerms.length + " (format: a + b*T)");
                emit(progressSink, "  ECI terms: " + cec.cecTerms.length + " (format: a + b*T)");
                logFullCecTable(cec, progressSink);
                logFirstNnPairCMatrix(clusterData, progressSink);
                logFullCmatForVerification(clusterData, cec, progressSink);
            } else {
                LOG.fine("    ECI terms: none");
                emit(progressSink, "  ECI terms: none");
            }
        }

        String systemName = cec.elements + "_" + cec.structurePhase;
        LOG.fine("STAGE 3c: Create ThermodynamicData bundle");
        LOG.fine("  systemName: " + systemName);
        emit(progressSink, "STAGE 3c: Create ThermodynamicData bundle: " + systemName);

        return new ThermodynamicData(clusterData, cec, clusterId, systemName);
    }

    /**
     * Prepares the data required for a thermodynamic calculation.
     */
    public ThermodynamicData prepareCalculation(ThermodynamicRequest request) throws Exception {

        return loadThermodynamicData(request.clusterId, request.hamiltonianId, request.progressSink);
    }

    /**
     * Runs a thermodynamic calculation.
     */
    public ThermodynamicResult runCalculation(ThermodynamicRequest request) throws Exception {

        LOG.info("ThermodynamicWorkflow.runCalculation — ENTER");
        LOG.info("  clusterId: " + request.clusterId);
        LOG.info("  hamiltonianId: " + request.hamiltonianId);
        LOG.info("  temperature: " + request.temperature + " K");
        LOG.info("  composition: " + Arrays.toString(request.composition));
        LOG.info("  engineType: " + request.engineType);
        LOG.info("  cvmBasisMode: " + request.cvmBasisMode);

        String resolvedHamiltonianId = resolveHamiltonianIdForEngine(
                request.hamiltonianId, request.engineType, request.cvmBasisMode, request.progressSink);
        ThermodynamicData data = loadThermodynamicData(request.clusterId, resolvedHamiltonianId, request.progressSink);

        ThermodynamicInput input = new ThermodynamicInput(
                data.clusterData,
                data.cec,
                request.temperature,
                request.composition,
                data.systemId,
                data.systemName,
                request.progressSink,
                request.eventSink,
                request.mcsL,
                request.mcsNEquil,
                request.mcsNAvg
        );

        ThermodynamicEngine engine;

        switch (request.engineType) {
            case "CVM":
                engine = cvmEngine;
                break;
            case "MCS":
                engine = mcsEngine;
                break;
            default:
                throw new IllegalArgumentException("Unknown engine: " + request.engineType);
        }

        EquilibriumState state = engine.compute(input);
        ThermodynamicResult result = ThermodynamicResult.from(state);

        LOG.info("ThermodynamicWorkflow.runCalculation — EXIT: G=" + String.format("%.4e", state.getFreeEnergy())
                + " J/mol");
        return result;
    }

    /**
     * For CVM runs prefer CVCF Hamiltonians when available.
     * Examples:
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_T_CVCF (preferred)
     *   Nb-Ti_BCC_A2_T      -> Nb-Ti_BCC_A2_CVCF   (legacy fallback)
     */
    private String resolveHamiltonianIdForEngine(
            String requestedHamiltonianId,
            String engineType,
            String cvmBasisMode,
            Consumer<String> progressSink) {
        if (!"CVM".equalsIgnoreCase(engineType)) {
            return requestedHamiltonianId;
        }
        if (requestedHamiltonianId == null || requestedHamiltonianId.isBlank()) {
            return requestedHamiltonianId;
        }
        if ("ORTHO".equalsIgnoreCase(cvmBasisMode)) {
            emit(progressSink, "CVM mode: using ORTHO Hamiltonian '" + requestedHamiltonianId + "'");
            LOG.info("CVM mode: using ORTHO Hamiltonian '" + requestedHamiltonianId + "'");
            return requestedHamiltonianId;
        }
        if (requestedHamiltonianId.endsWith("_CVCF")) {
            return requestedHamiltonianId;
        }

        String preferredId = requestedHamiltonianId + "_CVCF";
        if (cecWorkflow.hamiltonianExists(preferredId)) {
            emit(progressSink, "CVM mode: using CVCF Hamiltonian '" + preferredId
                    + "' instead of '" + requestedHamiltonianId + "'");
            LOG.info("CVM mode: using CVCF Hamiltonian '" + preferredId
                    + "' instead of '" + requestedHamiltonianId + "'");
            return preferredId;
        }

        int lastUnderscore = requestedHamiltonianId.lastIndexOf('_');
        if (lastUnderscore > 0) {
            String legacyId = requestedHamiltonianId.substring(0, lastUnderscore) + "_CVCF";
            if (cecWorkflow.hamiltonianExists(legacyId)) {
                emit(progressSink, "CVM mode: using CVCF Hamiltonian '" + legacyId
                        + "' instead of '" + requestedHamiltonianId + "'");
                LOG.info("CVM mode: using CVCF Hamiltonian '" + legacyId
                        + "' instead of '" + requestedHamiltonianId + "'");
                return legacyId;
            }
        }

        emit(progressSink, "CVM mode: CVCF Hamiltonian not found (tried '" + preferredId
                + "' and legacy pattern), falling back to '" + requestedHamiltonianId + "'");
        LOG.warning("CVM mode: CVCF Hamiltonian not found (tried '" + preferredId
                + "' and legacy pattern), falling back to '" + requestedHamiltonianId + "'");
        return requestedHamiltonianId;
    }

    /**
     * Logs the full CEC table in Stage 3.
     */
    private void logFullCecTable(CECEntry cec, Consumer<String> progressSink) {
        LOG.fine("  STAGE 3b.1: Full CEC table (name | a | b)");
        emit(progressSink, "  STAGE 3b.1: Full CEC table (name | a | b)");
        for (int i = 0; i < cec.cecTerms.length; i++) {
            CECTerm term = cec.cecTerms[i];
            if (term == null) {
                continue;
            }
            String row = String.format("    [%02d] %-12s a=% .8e  b=% .8e",
                    i, safeName(term), term.a, term.b);
            LOG.fine(row);
            emit(progressSink, row);
        }
    }

    /**
     * Logs the first-NN pair c-matrix block (not CEC matrix) for Stage 3b.2.
     */
    private void logFirstNnPairCMatrix(AllClusterData clusterData, Consumer<String> progressSink) {
        if (clusterData == null || clusterData.getCMatrixResult() == null) {
            emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix: unavailable (missing cMatrixResult)");
            return;
        }

        CMatrixResult cm = clusterData.getCMatrixResult();
        List<List<double[][]>> cmat = cm.getCmat();
        if (cmat == null || cmat.isEmpty()) {
            emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix: unavailable (empty cmat)");
            return;
        }

        int tPair = findFirstPairTypeIndex(clusterData);
        if (tPair < 0 || tPair >= cmat.size()) {
            emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix: no pair cluster type found");
            return;
        }

        List<double[][]> pairTypeBlocks = cmat.get(tPair);
        if (pairTypeBlocks == null || pairTypeBlocks.isEmpty()) {
            emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix: no blocks for pair cluster type t=" + tPair);
            return;
        }

        int blockIndex = -1;
        double[][] block = null;
        for (int j = 0; j < pairTypeBlocks.size(); j++) {
            double[][] candidate = pairTypeBlocks.get(j);
            if (candidate != null && candidate.length > 0 && candidate[0] != null && candidate[0].length > 0) {
                blockIndex = j;
                block = candidate;
                break;
            }
        }
        if (blockIndex < 0 || block == null) {
            emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix: no non-empty block for pair cluster type t=" + tPair);
            return;
        }

        int rows = block.length;
        int cols = block[0].length;
        emit(progressSink, "  STAGE 3b.2: First-NN pair c-matrix block (clusterType=" + tPair
                + ", orbitIndex=" + blockIndex + ")");
        emit(progressSink, "               dims: rows=" + rows + ", cols=" + cols
                + " (last col is constant term)");

        int[] rowWeights = null;
        List<List<int[]>> wcv = cm.getWcv();
        if (wcv != null && tPair < wcv.size() && wcv.get(tPair) != null && blockIndex < wcv.get(tPair).size()) {
            rowWeights = wcv.get(tPair).get(blockIndex);
        }

        int[][] lcv = cm.getLcv();
        if (lcv != null && tPair < lcv.length && lcv[tPair] != null) {
            emit(progressSink, "               lcv[t]=" + Arrays.toString(lcv[tPair]));
        }

        int maxRowsToPrint = 12;
        int maxColsToPrint = 12;
        int rowsToPrint = Math.min(rows, maxRowsToPrint);
        int colsToPrint = Math.min(cols, maxColsToPrint);
        emit(progressSink, "               preview: first " + rowsToPrint + " rows x " + colsToPrint + " cols");

        for (int r = 0; r < rowsToPrint; r++) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("               r=%-3d", r));
            if (rowWeights != null && r < rowWeights.length) {
                sb.append(String.format(" w=%-4d", rowWeights[r]));
            }
            sb.append(" |");
            for (int c = 0; c < colsToPrint; c++) {
                sb.append(String.format(" % .5e", block[r][c]));
            }
            if (cols > colsToPrint) {
                sb.append(" ...");
            }
            emit(progressSink, sb.toString());
        }
        if (rows > rowsToPrint) {
            emit(progressSink, "               ... (" + (rows - rowsToPrint) + " more rows)");
        }
    }

    private int findFirstPairTypeIndex(AllClusterData clusterData) {
        ClusterIdentificationResult disordered = clusterData.getDisorderedClusterResult();
        if (disordered == null) {
            return -1;
        }
        ClusCoordListResult disClusterData = disordered.getDisClusterData();
        if (disClusterData == null || disClusterData.getClusCoordList() == null) {
            return -1;
        }
        List<Cluster> clusTypes = disClusterData.getClusCoordList();
        for (int t = 0; t < clusTypes.size(); t++) {
            Cluster cluster = clusTypes.get(t);
            if (cluster != null && cluster.getAllSites().size() == 2) {
                return t;
            }
        }
        return -1;
    }

    private static String safeName(CECTerm term) {
        return term.name == null || term.name.isBlank() ? "(unnamed)" : term.name;
    }

    private static void emit(Consumer<String> progressSink, String line) {
        if (progressSink != null) {
            progressSink.accept(line);
        }
    }

    /**
     * Debug-only: emits full CVCF c-matrix and reconstructed orthogonal c-matrix
     * for BCC_A2 binary systems so users can compare against Mathematica output.
     */
    private void logFullCmatForVerification(AllClusterData clusterData, CECEntry cec, Consumer<String> progressSink) {
        if (clusterData == null || clusterData.getCMatrixResult() == null || cec == null) {
            return;
        }
        if (!"BCC_A2".equalsIgnoreCase(cec.structurePhase)) {
            return;
        }
        if (cec.elements == null || cec.elements.split("-").length != 2) {
            return;
        }

        CMatrixResult cm = clusterData.getCMatrixResult();
        List<List<double[][]>> cvcf = cm.getCmat();
        if (cvcf == null || cvcf.isEmpty()) {
            return;
        }

        emit(progressSink, "  STAGE 3b.3: Full CVCF c-matrix (all t,j blocks)");
        dumpFullCmat("CVCF", cvcf, progressSink);

        CvCfBasis basis = BccA2TModelCvCfTransformations.binaryBasis();
        if (basis.Tinv == null || basis.Tinv.length == 0) {
            emit(progressSink, "  STAGE 3b.4: Orthogonal c-matrix reconstruction skipped (missing Tinv)");
            return;
        }
        List<List<double[][]>> orth;
        try {
            orth = reconstructOrthCmat(cvcf, basis.Tinv);
        } catch (IllegalArgumentException e) {
            // In CVCF-first mode, orthogonal c-matrix reconstruction is optional (debug only)
            emit(progressSink, "  STAGE 3b.4: Orthogonal c-matrix reconstruction skipped (CVCF-first mode, not needed for Type-2)");
            return;
        }
        emit(progressSink, "  STAGE 3b.4: Full orthogonal c-matrix (reconstructed from CVCF via Tinv)");
        dumpFullCmat("ORTHO", orth, progressSink);
    }

    private void dumpFullCmat(String label, List<List<double[][]>> cmat, Consumer<String> progressSink) {
        for (int t = 0; t < cmat.size(); t++) {
            List<double[][]> groups = cmat.get(t);
            if (groups == null) {
                continue;
            }
            for (int j = 0; j < groups.size(); j++) {
                double[][] block = groups.get(j);
                if (block == null || block.length == 0 || block[0] == null) {
                    continue;
                }
                emit(progressSink, String.format("    %s cmat[t=%d][j=%d] dims=%dx%d",
                        label, t, j, block.length, block[0].length));
                for (int r = 0; r < block.length; r++) {
                    StringBuilder sb = new StringBuilder("      ");
                    sb.append(String.format("r=%-3d |", r));
                    for (int c = 0; c < block[r].length; c++) {
                        sb.append(String.format(" % .8e", block[r][c]));
                    }
                    emit(progressSink, sb.toString());
                }
            }
        }
    }

    private List<List<double[][]>> reconstructOrthCmat(List<List<double[][]>> cvcf, double[][] tInv) {
        int tcfOld = tInv.length;
        int tcfNew = tInv[0].length;
        List<List<double[][]>> result = new ArrayList<>(cvcf.size());

        for (List<double[][]> groups : cvcf) {
            List<double[][]> outGroups = new ArrayList<>(groups.size());
            for (double[][] block : groups) {
                int rows = block.length;
                int cols = block[0].length;
                int cvcfCfCols = cols - 1;
                if (cvcfCfCols != tcfNew) {
                    throw new IllegalArgumentException("Cannot reconstruct orthogonal c-matrix: CVCF CF cols="
                            + cvcfCfCols + ", expected " + tcfNew);
                }
                double[][] orth = new double[rows][tcfOld + 1];
                for (int r = 0; r < rows; r++) {
                    for (int kOld = 0; kOld < tcfOld; kOld++) {
                        double sum = 0.0;
                        for (int kNew = 0; kNew < tcfNew; kNew++) {
                            sum += block[r][kNew] * tInv[kOld][kNew];
                        }
                        orth[r][kOld] = sum;
                    }
                    orth[r][tcfOld] = block[r][cols - 1];
                }
                outGroups.add(orth);
            }
            result.add(outGroups);
        }
        return result;
    }

}
