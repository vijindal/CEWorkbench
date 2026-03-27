package org.ce.domain.cluster.cvcf;

import org.ce.domain.cluster.CMatrixResult;

import java.util.ArrayList;
import java.util.List;

/**
 * Transforms an old orthogonal-basis {@link CMatrixResult} to the CVCF basis.
 *
 * <p>Given the relationship {@code u_old[i] = Σ_j T[i][j] · v_new[j]},
 * the cluster variable evaluation:
 * <pre>
 *   CV[v] = Σ_{k_old} cmat_old[v][k_old] · u_old[k_old] + cmat_old[v][tcf_old]
 * </pre>
 * becomes:
 * <pre>
 *   CV[v] = Σ_{k_new} cmat_new[v][k_new] · v_new[k_new] + cmat_new[v][tcf_new]
 * </pre>
 * where:
 * <pre>
 *   cmat_new[v][k_new] = Σ_{k_old} cmat_old[v][k_old] · T[k_old][k_new]
 *   cmat_new[v][tcf_new] = cmat_old[v][tcf_old]   (constant term unchanged)
 * </pre>
 * i.e., {@code cmat_new[:, :tcf_new] = cmat_old[:, :tcf_old] · T}.
 */
public final class CvCfBasisTransformer {

    private CvCfBasisTransformer() {}

    /**
     * Transforms a CMatrixResult from the old orthogonal basis to the CVCF basis.
     *
     * @param oldResult  the CMatrixResult in the old orthogonal basis
     * @param basis      the CVCF basis description, including T matrix and CF names
     * @return a new CMatrixResult in the CVCF basis
     */
    public static CMatrixResult transform(CMatrixResult oldResult, CvCfBasis basis) {
        List<List<double[][]>> oldCmat = oldResult.getCmat();
        int[][] oldLcv = oldResult.getLcv();
        List<List<int[]>> wcv = oldResult.getWcv();

        int tcdis = oldCmat.size();
        double[][] T = basis.T;
        int tcfOld = T.length;          // number of old CF columns (non-point + point)
        int tcfNew = basis.totalCfs();  // number of new CF columns

        List<List<double[][]>> newCmat = new ArrayList<>(tcdis);

        for (int t = 0; t < tcdis; t++) {
            List<double[][]> oldGroup = oldCmat.get(t);
            int lc = oldGroup.size();
            List<double[][]> newGroup = new ArrayList<>(lc);

            for (int j = 0; j < lc; j++) {
                double[][] cmOld = oldGroup.get(j);  // [nv][numOldCols]
                int nv = cmOld.length;
                int numOldCols = (nv > 0) ? cmOld[0].length : 0;

                // The last column of old C-matrix is the constant term
                // Structure: [CF0, CF1, ..., CFn, constant]
                // where constant column contains the fixed term (usually same value across rows)

                // New cmat: tcfNew CF columns + 1 constant column
                double[][] cmNew = new double[nv][tcfNew + 1];

                for (int v = 0; v < nv; v++) {
                    // Multiply: cmat_new[v][k_new] = Σ_{k_old} cmat_old[v][k_old] · T[k_old][k_new]
                    for (int kNew = 0; kNew < tcfNew; kNew++) {
                        double sum = 0.0;
                        for (int kOld = 0; kOld < tcfOld && kOld < numOldCols - 1; kOld++) {
                            sum += cmOld[v][kOld] * T[kOld][kNew];
                        }
                        cmNew[v][kNew] = sum;
                    }
                    // Preserve constant term from last column of old matrix
                    cmNew[v][tcfNew] = cmOld[v][numOldCols - 1];
                }
                newGroup.add(cmNew);
            }
            newCmat.add(newGroup);
        }

        // cfBasisIndices is not used in the CVCF basis (CF columns identified by name)
        return new CMatrixResult(newCmat, oldLcv, wcv, null);
    }
}
