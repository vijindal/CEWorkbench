package org.ce.workflow;

import static org.ce.domain.cluster.SpaceGroup.SymmetryOperation;
import static org.ce.domain.cluster.ClusterResults.*;

import org.ce.domain.cluster.*;
import org.ce.domain.cluster.cvcf.CvCfBasis;
import org.ce.storage.InputLoader;
import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * Workflow for cluster and correlation function identification.
 *
 * <p>
 * This workflow orchestrates the complete identification pipeline for both
 * disordered (highest-symmetry) and ordered phases:
 * </p>
 * <ol>
 * <li>Load disordered cluster geometry and symmetry (InputLoader)</li>
 * <li>Load ordered cluster geometry and symmetry (InputLoader)</li>
 * <li>Identify clusters in both phases (ClusterIdentifier)</li>
 * <li>Identify correlation functions in both phases (CFIdentifier)</li>
 * <li>Return complete results as AllClusterData</li>
 * </ol>
 *
 * @author CVM Project
 * @version 2.0
 */
public class ClusterIdentificationWorkflow {

        private static final Logger LOG = Logger.getLogger(ClusterIdentificationWorkflow.class.getName());

        private ClusterIdentificationWorkflow() {
        }

        /**
         * Executes the complete cluster and correlation function identification
         * workflow.
         *
         * @param config the identification request with all configuration parameters
         * @return bundle of identification results containing both disordered and
         *         ordered data
         * @throws RuntimeException if resources are not found or identification fails
         */
        public static AllClusterData identify(ClusterIdentificationRequest config) {
                return identify(config, null);
        }

        public static AllClusterData identify(ClusterIdentificationRequest config, Consumer<String> progressSink) {
                LOG.info("ClusterIdentificationWorkflow.identify — START");
                emit(progressSink, "TYPE-1a START");

                // Derive a basic System ID for display
                String compSuffix = config.getNumComponents() == 2 ? "bin"
                                : config.getNumComponents() == 3 ? "tern"
                                                : config.getNumComponents() == 4 ? "quat"
                                                                : "K" + config.getNumComponents();
                String systemId = config.getStructurePhase() + "_" + compSuffix;

                // emit(progressSink, " System ID : " + systemId);
                // emit(progressSink, " Number of Components : " + config.getNumComponents());
                // emit(progressSink, " Input Configuration:");
                // emit(progressSink, " disorderedClusterFile : " +
                // config.getDisorderedClusterFile());
                // emit(progressSink, " disorderedSymmetryGroup: " +
                // config.getDisorderedSymmetryGroup());
                // emit(progressSink, " orderedClusterFile : " +
                // config.getOrderedClusterFile());
                // emit(progressSink, " orderedSymmetryGroup : " +
                // config.getOrderedSymmetryGroup());

                // Validate CVCF basis is supported for this (structure, numComponents)
                String structurePhase = config.getStructurePhase();
                int numComponents = config.getNumComponents();
                if (!CvCfBasis.Registry.INSTANCE.isSupported(structurePhase, numComponents)) {
                        String msg = "CVCF basis not supported for " + structurePhase
                                        + " with " + numComponents + " components.\n"
                                        + CvCfBasis.Registry.INSTANCE.supportedSummary();
                        LOG.severe(msg);
                        emit(progressSink, "ERROR: " + msg);
                        throw new IllegalArgumentException(msg);
                }
                emit(progressSink, "  basis supported: " + structurePhase + ", K=" + numComponents);

                // =====================================================================
                // 1. Load disordered phase (HSP) data
                // =====================================================================
                LOG.fine("Loading disordered phase data...");
                emit(progressSink, "STAGE 1a: Load disordered phase geometry/symmetry");

                List<Cluster> disorderedClusters = InputLoader.parseClusterFile(config.getDisorderedClusterFile());
                SpaceGroup disorderedSpaceGroup = InputLoader.parseSpaceGroup(config.getDisorderedSymmetryGroup());
                List<SymmetryOperation> disorderedSymOps = disorderedSpaceGroup.getOperations();

                LOG.fine("Disordered clusters loaded: " + disorderedClusters.size());
                emit(progressSink, "    [Metadata] disordered max clusters: " + disorderedClusters.size());
                emit(progressSink, "    [Metadata] disordered sym ops     : " + disorderedSymOps.size());

                // =====================================================================
                // 2. Load ordered phase data
                // =====================================================================
                LOG.fine("Loading ordered phase data...");
                emit(progressSink, "STAGE 1b: Load ordered phase geometry/symmetry");

                List<Cluster> orderedClusters = InputLoader.parseClusterFile(config.getOrderedClusterFile());
                SpaceGroup orderedSpaceGroup = InputLoader.parseSpaceGroup(config.getOrderedSymmetryGroup());
                List<SymmetryOperation> orderedSymOps = orderedSpaceGroup.getOperations();

                LOG.fine("Ordered clusters loaded: " + orderedClusters.size());
                emit(progressSink, "    [Metadata] ordered max clusters   : " + orderedClusters.size());
                emit(progressSink, "    [Metadata] ordered sym ops        : " + orderedSymOps.size());

                // =====================================================================
                // 3. Identify clusters (both disordered and ordered phases)
                // =====================================================================
                LOG.fine("Identifying clusters for both phases...");
                emit(progressSink, "STAGE 1c: Cluster identification (HSP and phase analysis)");

                ClusterIdentificationResult clusterResult = ClusterIdentifier.identify(
                                disorderedClusters,
                                disorderedSymOps,
                                orderedClusters,
                                orderedSymOps,
                                config.getTransformationMatrix(),
                                new double[] { config.getTranslationVector().getX(),
                                                config.getTranslationVector().getY(),
                                                config.getTranslationVector().getZ() });

                LOG.fine("Clusters identified: tcdis=" + clusterResult.getTcdis() + ", tc=" + clusterResult.getTc());
                emit(progressSink, "  STAGE 1 diagnostic results:");

                // Stage 1a: HSP (Disordered) Clusters
                emit(progressSink, "    [Stage 1a] HSP Clusters (Disordered Phase):");
                emit(progressSink, "      tcdis  = " + clusterResult.getTcdis());
                emit(progressSink, "      nxcdis = " + clusterResult.getNxcdis());
                List<Double> mhdis = clusterResult.getDisClusterData().getMultiplicities();
                List<List<Integer>> rcdis = clusterResult.getDisClusterData().getRcList();
                List<Cluster> clusdis = clusterResult.getDisClusterData().getClusCoordList();
                for (int t = 0; t < clusterResult.getTcdis(); t++) {
                        int nc = clusdis.get(t).getAllSites().size();
                        emit(progressSink, String.format("        t=%-4d nc=%-2d rc=%s mhdis=%.4f kbdis=%.8f",
                                        t, nc, rcdis.get(t), mhdis.get(t), clusterResult.getKbCoefficients()[t]));
                }

                // Stage 1b: Phase (Ordered) Clusters
                emit(progressSink, "    [Stage 1b] Ordered Phase Clusters:");
                emit(progressSink, "      tc     = " + clusterResult.getTc());
                emit(progressSink, "      nxc    = " + clusterResult.getNxc());
                emit(progressSink, "      lc     = " + Arrays.toString(clusterResult.getLc()));

                ClassifiedClusterResult ordData = clusterResult.getOrdClusterData();
                double[][] mh = clusterResult.getMh();
                for (int t = 0; t < clusterResult.getTcdis(); t++) {
                        for (int j = 0; j < clusterResult.getLc()[t]; j++) {
                                int nc = ordData.getCoordList().get(t).get(j).getAllSites().size();
                                List<Integer> rc = ordData.getRcList().get(t).get(j);
                                emit(progressSink, String.format("        t=%-4d j=%-4d nc=%-2d rc=%s mh=%.4f",
                                                t, j, nc, rc, mh[t][j]));
                        }
                }

                // =====================================================================
                // 4. Identify CFs (both disordered and ordered phases)
                // =====================================================================
                LOG.fine("Identifying correlation functions for both phases...");
                emit(progressSink, "STAGE 2: CF identification");

                CFIdentificationResult cfResult = CFIdentifier.identify(
                                clusterResult,
                                clusterResult.getDisClusterData().getClusCoordList(),
                                disorderedSymOps,
                                orderedClusters,
                                orderedSymOps,
                                orderedSpaceGroup.getRotateMat(),
                                orderedSpaceGroup.getTranslateMat(),
                                config.getNumComponents());

                LOG.fine("CFs identified: tcfdis=" + cfResult.getTcfdis() + ", tcf=" + cfResult.getTcf());

                // Stage 2a: HSP (Disordered) CFs
                emit(progressSink, "    [Stage 2a] HSP CFs (Disordered Phase):");
                emit(progressSink, "      tcfdis = " + cfResult.getTcfdis());
                List<Double> mcfdis = cfResult.getDisCFData().getMultiplicities();
                List<List<Integer>> rcfdis = cfResult.getDisCFData().getRcList();
                List<Cluster> cluscfdis = cfResult.getDisCFData().getClusCoordList();
                for (int i = 0; i < cfResult.getTcfdis(); i++) {
                        String def = formatCfDef(cluscfdis.get(i));
                        emit(progressSink, String.format("        i=%-4d mcfdis=%.4f rcfdis=%-10s def=%s",
                                        i, mcfdis.get(i), rcfdis.get(i), def));
                }

                // Stage 2b: Phase (Ordered) CFs
                emit(progressSink, "    [Stage 2b] Ordered Phase CFs:");
                emit(progressSink, "      tcf    = " + cfResult.getTcf());
                emit(progressSink, "      ncf    = " + cfResult.getNcf());
                emit(progressSink, "      lcf    = " + Arrays.deepToString(cfResult.getLcf()));

                GroupedCFResult groupedCFData = cfResult.getGroupedCFData();
                List<List<List<Double>>> mcf = groupedCFData.getMultiplicityData();
                List<List<List<List<Integer>>>> rcf = groupedCFData.getRcData();
                List<List<List<Cluster>>> cluscf = groupedCFData.getCoordData();

                int tcdis = clusterResult.getTcdis();
                for (int t = 0; t < tcdis; t++) {
                        for (int j = 0; j < cfResult.getLcf()[t].length; j++) {
                                int numCfsInGroup = cfResult.getLcf()[t][j];
                                for (int p = 0; p < numCfsInGroup; p++) {
                                        String def = formatCfDef(cluscf.get(t).get(j).get(p));
                                        emit(progressSink, String.format(
                                                        "        t=%-4d j=%-4d p=%-4d mcf=%.4f rcf=%-10s def=%s",
                                                        t, j, p, mcf.get(t).get(j).get(p), rcf.get(t).get(j).get(p),
                                                        def));
                                }
                        }
                }

                // Stage 2c: CF Names
                emit(progressSink, "    [Stage 2c] Correlation Function Names (uNames):");
                List<String> uNames = cfResult.getUNames();
                for (int i = 0; i < uNames.size(); i++) {
                        emit(progressSink, String.format("      u[%d] = %s", i, uNames.get(i)));
                }

                emit(progressSink, "  EXPECTED c-matrix columns = tcf + 1 = " + (cfResult.getTcf() + 1)
                                + " (last column is constant term)");

                // =====================================================================
                // 5. Build C-Matrix (Stage 3: Orthogonal Basis)
                // =====================================================================
                LOG.fine("Building Orthogonal C-Matrix...");
                emit(progressSink, "STAGE 3: Build C-matrix (orthogonal basis)");

                CMatrix.Result orthMatrix = CMatrix.buildOrthogonal(
                                clusterResult,
                                cfResult,
                                orderedClusters,
                                numComponents);

                emit(progressSink, "  [Orthogonal Basis C-Matrix Data]");

                // R-Matrix and its Inverse (Vandermonde)
                double[][] rMat = ClusterMath.buildRMatrix(numComponents);
                emit(progressSink, "    R-Matrix (site operators <-> compositions):");
                for (int i = 0; i < rMat.length; i++) {
                        emit(progressSink, "      " + Arrays.toString(rMat[i]));
                }

                double[] basisSymbols = ClusterMath.buildBasis(numComponents);
                double[][] rInv = new double[numComponents][numComponents];
                for (int i = 0; i < numComponents; i++) {
                        for (int j = 0; j < numComponents; j++) {
                                rInv[i][j] = (i == 0) ? 1.0 : Math.pow(basisSymbols[j], i);
                        }
                }
                emit(progressSink, "    R-Matrix Inverse (Vandermonde matrix V_ij = s_j^i):");
                for (int i = 0; i < rInv.length; i++) {
                        emit(progressSink, "      " + Arrays.toString(rInv[i]));
                }

                double[] equiatomicX = new double[numComponents];
                Arrays.fill(equiatomicX, 1.0 / numComponents);

                // Site Operator values (point CFs) = V * x
                double[] siteOps = new double[numComponents];
                for (int i = 0; i < numComponents; i++) {
                        for (int j = 0; j < numComponents; j++) {
                                siteOps[i] += rInv[i][j] * equiatomicX[j];
                        }
                }
                emit(progressSink, "    Site Operator values (site-ops/point-CFs) at equiatomic composition:");
                emit(progressSink, "      " + Arrays.toString(siteOps));

                // Random CFs
                int totalCfs = cfResult.getTcf();
                int[][] cfBasisIndices = orthMatrix.getCfBasisIndices();
                double[] fullRandomCFs = new double[totalCfs + 1];

                // First tcf values are products of site operators
                for (int i = 0; i < totalCfs; i++) {
                        int[] indices = cfBasisIndices[i];
                        double val = 1.0;
                        for (int b : indices) {
                                // basis indices are 1-based in Mathematica/CMatrix logic, pointing to s1, s2...
                                // so b=1 corresponds to siteOps[1]
                                val *= siteOps[b];
                        }
                        fullRandomCFs[i] = val;
                }
                // Last value is always identity (1.0)
                fullRandomCFs[totalCfs] = 1.0;

                emit(progressSink, "    Random CF list (orthogonal basis, length=" + (totalCfs + 1) + "):");
                emit(progressSink, "      " + Arrays.toString(fullRandomCFs));

                int[][] lcv = orthMatrix.getLcv();
                List<List<int[]>> wcv = orthMatrix.getWcv();
                List<List<double[][]>> cmat = orthMatrix.getCmat();

                double[][] orthoCvs = new double[cmat.size()][];
                for (int t = 0; t < cmat.size(); t++) {
                        orthoCvs[t] = new double[cmat.get(t).size()];
                        for (int j = 0; j < cmat.get(t).size(); j++) {
                                double[][] block = cmat.get(t).get(j);
                                double commonCv = 0;
                                boolean first = true;
                                boolean allEqual = true;

                                for (int r = 0; r < block.length; r++) {
                                        double cvVal = 0;
                                        for (int c = 0; c < block[r].length; c++) {
                                                cvVal += block[r][c] * fullRandomCFs[c];
                                        }
                                        if (first) {
                                                commonCv = cvVal;
                                                first = false;
                                        } else if (Math.abs(cvVal - commonCv) > 1e-10) {
                                                allEqual = false;
                                        }
                                }
                                orthoCvs[t][j] = commonCv;

                                emit(progressSink, String.format(
                                                "    Group t=%d j=%d: lcv=%d wcv=%s => CV(random) = %.10f %s",
                                                t, j, lcv[t][j], Arrays.toString(wcv.get(t).get(j)), commonCv,
                                                allEqual ? "(Verified: all rows consistent)"
                                                                : "(WARNING: rows inconsistent!)"));
                        }
                }

                // =====================================================================
                // 6. CVCF Transformation (Stage 4)
                // =====================================================================
                LOG.fine("Transforming to CVCF Basis...");
                emit(progressSink, "STAGE 4: CVCF Basis Transformation");

                CvCfBasis cvcfBasis = CvCfBasis.Registry.INSTANCE.get(config.getStructurePhase(), numComponents);

                // emit(progressSink, " [Stage 4a] CVCF Symbol Names and Transformation
                // Rules:");
                List<String> vNames = cvcfBasis.cfNames;
                double[][] T = cvcfBasis.T;
                for (int i = 0; i < uNames.size(); i++) {
                        StringBuilder rule = new StringBuilder();
                        rule.append(uNames.get(i)).append(" = ");
                        boolean firstTerm = true;
                        for (int j = 0; j < vNames.size(); j++) {
                                double coeff = T[i][j];
                                if (Math.abs(coeff) > 1e-10) {
                                        if (!firstTerm && coeff > 0)
                                                rule.append(" + ");
                                        if (coeff < 0)
                                                rule.append(" - ");
                                        double absCoeff = Math.abs(coeff);
                                        if (Math.abs(absCoeff - 1.0) > 1e-10)
                                                rule.append(String.format("%.4f*", absCoeff));
                                        rule.append(vNames.get(j));
                                        firstTerm = false;
                                }
                        }
                        // emit(progressSink, " " + rule.toString());
                }

                emit(progressSink, "\n    [Stage 4b] Transformation Matrix T (Orthogonal -> CVCF):");
                for (int i = 0; i < T.length; i++) {
                        // emit(progressSink, " " + Arrays.toString(T[i]));
                }

                // Random CFs in CVCF basis: v = Tinv * u
                double[][] Tinv = cvcfBasis.Tinv;
                double[] randomCVCFs = new double[cvcfBasis.totalCfs()];
                for (int j = 0; j < randomCVCFs.length; j++) {
                        for (int i = 0; i < fullRandomCFs.length; i++) {
                                randomCVCFs[j] += Tinv[j][i] * fullRandomCFs[i];
                        }
                }
                emit(progressSink, "\n    [Stage 4c] Random CF list (CVCF basis):");
                emit(progressSink, "      " + Arrays.toString(randomCVCFs));

                CMatrix.Result cMatrix = org.ce.domain.cluster.cvcf.CvCfBasisTransformer.transform(orthMatrix,
                                cvcfBasis);

                emit(progressSink, "\n    [Stage 4d] CVCF Basis C-Matrix Data and Verification:");
                int[][] lcv_cvcf = cMatrix.getLcv();
                List<List<int[]>> wcv_cvcf = cMatrix.getWcv();
                List<List<double[][]>> cmat_cvcf = cMatrix.getCmat();

                for (int t = 0; t < cmat_cvcf.size(); t++) {
                        for (int j = 0; j < cmat_cvcf.get(t).size(); j++) {
                                double[][] block = cmat_cvcf.get(t).get(j);
                                double commonCv = 0;
                                boolean first = true;
                                boolean allEqual = true;

                                for (int r = 0; r < block.length; r++) {
                                        double cvVal = 0;
                                        for (int c = 0; c < block[r].length; c++) {
                                                cvVal += block[r][c] * randomCVCFs[c];
                                        }
                                        if (first) {
                                                commonCv = cvVal;
                                                first = false;
                                        } else if (Math.abs(cvVal - commonCv) > 1e-10) {
                                                allEqual = false;
                                        }
                                }

                                boolean basisConsistent = Math.abs(commonCv - orthoCvs[t][j]) < 1e-10;

                                emit(progressSink, String.format(
                                                "    Group t=%d j=%d (CVCF): lcv=%d wcv=%s => CV(random) = %.10f %s",
                                                t, j, lcv_cvcf[t][j], Arrays.toString(wcv_cvcf.get(t).get(j)), commonCv,
                                                (allEqual && basisConsistent) ? "(Verified)"
                                                                : "(WARNING: Inconsistent!)"));
                        }
                }

                // =====================================================================
                // 7. Generate basis-specific symbol lists
                // =====================================================================
                LOG.fine("Retrieving symbolic lists for CFs and CECs from Stage 2 result...");
                CvCfBasis basis = cvcfBasis;

                List<String> uList = cfResult.getUNames();
                List<String> eOList = cfResult.getEONames();
                List<String> vList = basis.cfNames;
                List<String> eList = basis.eciNames;

                // =====================================================================
                // 7. Package and return results
                // =====================================================================
                // Note: ClusterIdentificationResult and CFIdentificationResult each contain
                // analysis of both disordered and ordered phases. We pass them as both
                // the "disordered" and "ordered" parameters to maintain architectural clarity
                // about the 4 distinct inputs (dis cluster file, dis symmetry, ord cluster
                // file,
                // ord symmetry), even though the results derive from a single analysis.
                AllClusterData finalResult = new AllClusterData(
                                clusterResult, // disordered cluster result
                                clusterResult, // ordered cluster result
                                cfResult, // disordered CF result
                                cfResult, // ordered CF result
                                cMatrix // C-Matrix result
                                , uList, vList, eOList, eList);

                LOG.info("ClusterIdentificationWorkflow.identify — EXIT: " + finalResult.getSummary());

                return finalResult;
        }

        private static void dump2DInt(String label, int[][] table, Consumer<String> sink) {
                if (table == null) {
                        emit(sink, label + " = null");
                        return;
                }
                for (int i = 0; i < table.length; i++) {
                        emit(sink, label + "[" + i + "] = " + Arrays.toString(table[i]));
                }
        }

        private static String formatCfDef(Cluster cluster) {
                if (cluster == null || cluster.getAllSites().isEmpty()) {
                        return "1"; // Identity / constant term
                }
                List<String> ops = new ArrayList<>();
                for (ClusterPrimitives.Site site : cluster.getAllSites()) {
                        if (site.getSymbol() != null && !site.getSymbol().isBlank()) {
                                ops.add(site.getSymbol());
                        }
                }
                return ops.isEmpty() ? "1" : String.join("*", ops);
        }

        private static void emit(Consumer<String> sink, String line) {
                if (sink != null) {
                        sink.accept(line);
                }
        }
}
