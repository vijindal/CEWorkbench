package org.ce.debug;

import org.ce.model.cvm.CVMGibbsModel;
import org.ce.model.cvm.CVMGibbsModel.ModelResult;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.cvm.CvCfBasis.RandomStateVectors;
import org.ce.model.cluster.ClusterMath;
import org.ce.model.PhysicsConstants;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;
import org.ce.model.hamiltonian.CECEntry;

import java.util.Arrays;

/**
 * Debug runner: evaluates G/H/S at the random (disordered) CFs without N-R minimization.
 * Prints all intermediate vectors: point basis functions → orthogonal CFs → CVCF CFs → cv values.
 *
 * Usage: randomEval <elements> <structure> <model> <T> <x1> <x2> [x3 ...]
 */
public class RandomStateEvaluator {

    public static void main(String[] args) throws Exception {
        if (args.length < 6) {
            System.err.println("Usage: randomEval <elements> <structure> <model> <T> <x1> <x2> [x3 ...]");
            System.exit(1);
        }

        String elements  = args[0];
        String structure = args[1];
        String model     = args[2];
        double T         = Double.parseDouble(args[3]);

        int ncomp = elements.split("-").length;
        double[] x = new double[ncomp];
        for (int i = 0; i < ncomp; i++) x[i] = Double.parseDouble(args[4 + i]);

        double sum = Arrays.stream(x).sum();
        if (Math.abs(sum - 1.0) > 1e-9) {
            System.err.printf("Composition sums to %.6f, not 1.0%n", sum);
            System.exit(1);
        }

        sep();
        System.out.println("  RANDOM STATE EVALUATION (no N-R minimization)");
        sep();
        System.out.printf("  System     : %s / %s / %s%n", elements, structure, model);
        System.out.printf("  Temperature: %.2f K%n", T);
        System.out.printf("  Composition: %s%n", Arrays.toString(x));
        System.out.println();

        // Load Hamiltonian
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore store = new DataStore.HamiltonianStore(ws);
        String hamiltonianId = elements + "_" + structure + "_" + model + "_CVCF";
        CECEntry cecEntry = store.load(hamiltonianId);
        System.out.println("  Hamiltonian: " + hamiltonianId);

        // Build model
        CVMGibbsModel cvm = new CVMGibbsModel();
        cvm.initialize(elements, structure, model, cecEntry, null);
        int ncf   = cvm.getNcf();
        int tcf   = cvm.getTcf();
        int tcdis = cvm.getTcdis();
        int[] lc  = cvm.getLc();
        int[][] lcv = cvm.getLcv();
        System.out.printf("  ncf=%d  tcf=%d  ncomp=%d  tcdis=%d%n%n", ncf, tcf, ncomp, tcdis);

        CvCfBasis basis = cvm.getBasis();
        int[][] orthIdx = cvm.getOrthCfBasisIndices();

        // ── RMAT: R-matrix (K×K) and basis vector ─────────────────────────
        double[] basisVec = ClusterMath.buildBasis(ncomp);
        double[][] rMat   = ClusterMath.buildRMatrix(ncomp);
        sep();
        System.out.printf("  RMAT — R-matrix (%dx%d)  and basis vector%n", ncomp, ncomp);
        sep();
        System.out.printf("  Basis vector σ: %s%n", Arrays.toString(basisVec));
        System.out.println("  R-matrix  (R = M⁻¹,  where M[i][j] = σ[j]^i):");
        for (int i = 0; i < rMat.length; i++) {
            System.out.printf("    R[%d][] = ", i);
            for (int j = 0; j < rMat[i].length; j++)
                System.out.printf("%+12.6f ", rMat[i][j]);
            System.out.println();
        }
        System.out.println();

        // ── T matrix (CVCF transformation) ────────────────────────────────
        double[][] tMat    = basis.T;
        double[][] tInvMat = (basis.Tinv != null) ? basis.Tinv
                           : org.ce.model.cluster.LinearAlgebra.invert(tMat);
        System.out.printf("  T-matrix (%dx%d)  (CVCF basis transform: u_cvcf = T⁻¹ · u_orth):%n",
                tMat.length, tMat[0].length);
        for (int i = 0; i < tMat.length; i++) {
            System.out.printf("    T[%2d][] = ", i);
            for (int j = 0; j < tMat[i].length; j++)
                System.out.printf("%+10.4f ", tMat[i][j]);
            System.out.println();
        }
        System.out.println();
        System.out.printf("  T⁻¹-matrix (%dx%d):%n", tInvMat.length, tInvMat[0].length);
        for (int i = 0; i < tInvMat.length; i++) {
            System.out.printf("    Tinv[%2d][] = ", i);
            for (int j = 0; j < tInvMat[i].length; j++)
                System.out.printf("%+10.4f ", tInvMat[i][j]);
            System.out.println();
        }
        System.out.println();

        // ── orthCfBasisIndices ─────────────────────────────────────────────
        System.out.printf("  orthCfBasisIndices (%d entries — topology of non-point CFs):%n", orthIdx.length);
        for (int i = 0; i < orthIdx.length; i++)
            System.out.printf("    idx[%2d] = %s%n", i, Arrays.toString(orthIdx[i]));
        System.out.println();

        // ── STEP 1: Point basis functions φₖ(x) and orthogonal random CFs ──
        RandomStateVectors rv = basis.computeRandomStateVectors(x, orthIdx);

        sep();
        System.out.println("  STEP 1 — Point basis functions  φₖ(x)  (K-1 values)");
        sep();
        System.out.println("  Formula: φ[k] = Σᵢ xᵢ · σᵢᵏ   (k = 1 .. K-1)");
        System.out.println();
        System.out.printf("  Basis vector σ (%d values for K=%d):%n", ncomp, ncomp);
        for (int i = 0; i < basisVec.length; i++)
            System.out.printf("    σ[%d] = %+.1f%n", i, basisVec[i]);
        System.out.println();
        System.out.println("  Composition x:");
        for (int i = 0; i < ncomp; i++)
            System.out.printf("    x[%d] = %.6f%n", i, x[i]);
        System.out.println();
        System.out.println("  Calculation term-by-term:");
        for (int k = 0; k < rv.pointCF.length; k++) {
            System.out.printf("    φ[%d]:  ", k + 1);
            double phiSum = 0.0;
            for (int i = 0; i < ncomp; i++) {
                double term = x[i] * Math.pow(basisVec[i], k + 1);
                phiSum += term;
                System.out.printf("x[%d](%.4f)·σ[%d]^%d(%.4f)=%+.6f",
                        i, x[i], i, k + 1, basisVec[i], term);
                if (i < ncomp - 1) System.out.print("  +  ");
            }
            System.out.printf("  =  %+.10f%n", phiSum);
        }
        System.out.println();
        System.out.println("  Result:");
        for (int k = 0; k < rv.pointCF.length; k++)
            System.out.printf("    φ[%d] = %+.10f%n", k + 1, rv.pointCF[k]);
        System.out.println();

        System.out.printf("  STEP 2 — Orthogonal non-point CFs  u_orth[0..%d]%n", rv.uOrthNonPoint.length - 1);
        System.out.println("  (products of φₖ values according to cluster topology)");
        for (int i = 0; i < rv.uOrthNonPoint.length; i++)
            System.out.printf("    u_orth[%2d] = %+.10f%n", i, rv.uOrthNonPoint[i]);
        System.out.println();

        System.out.printf("  STEP 3 — Full orthogonal vector  uOrthFull[0..%d]%n", rv.uOrthFull.length - 1);
        System.out.println("  (= [u_orth_nonpoint | pointCF | (1?)] — input to T⁻¹ transform)");
        for (int i = 0; i < rv.uOrthFull.length; i++)
            System.out.printf("    uOrthFull[%2d] = %+.10f%n", i, rv.uOrthFull[i]);
        System.out.println();

        System.out.printf("  STEP 4 — CVCF non-point CFs  u_cvcf[0..%d]  (= T⁻¹ · uOrthFull)%n", ncf - 1);
        System.out.println("  (these are what the N-R loop minimizes over)");
        for (int i = 0; i < rv.uCvcfNonPoint.length; i++)
            System.out.printf("    u_cvcf[%2d] = %+.10f%n", i, rv.uCvcfNonPoint[i]);
        System.out.println();

        System.out.printf("  STEP 5 — Full CVCF vector  uFull[0..%d]  (non-point + point appended)%n", tcf - 1);
        for (int i = 0; i < rv.uCvcfFull.length; i++) {
            String tag = (i < ncf) ? "non-point" : "point    ";
            System.out.printf("    uFull[%2d] (%s) = %+.10f%n", i, tag, rv.uCvcfFull[i]);
        }
        System.out.println();

        // ── Cluster variables ──────────────────────────────────────────────
        double[][][] cv = cvm.evaluateClusterVariables(rv.uCvcfNonPoint, x);
        System.out.println("  STEP 6 — Cluster variables  cv[t][j][v]:");
        System.out.printf("  (random state: should all be ∈ (0,1); tetra≈%.4f, tri≈%.4f, pair≈%.4f)%n",
                Math.pow(1.0/ncomp, 4), Math.pow(1.0/ncomp, 3), Math.pow(1.0/ncomp, 2));
        for (int t = 0; t < tcdis; t++) {
            for (int j = 0; j < lc[t]; j++) {
                System.out.printf("    t=%d j=%d  (%d vars):%n", t, j, lcv[t][j]);
                for (int v = 0; v < lcv[t][j]; v++) {
                    double cvVal = cv[t][j][v];
                    String flag  = (cvVal <= 0) ? "  *** NEGATIVE — entropy blow-up!" :
                                   (cvVal > 1)  ? "  *** > 1 — unphysical!" : "";
                    System.out.printf("      cv[%d][%d][%2d] = %+.8f%s%n", t, j, v, cvVal, flag);
                }
            }
        }
        System.out.println();

        // ── Physics ────────────────────────────────────────────────────────
        ModelResult r = cvm.evaluate(rv.uCvcfNonPoint, x, T);
        double sMixJmol = 0.0;
        for (double xi : x) if (xi > 0) sMixJmol -= xi * Math.log(xi);
        sMixJmol *= PhysicsConstants.R_GAS;

        sep();
        System.out.println("  PHYSICS SUMMARY");
        sep();
        System.out.printf("  H          = %+20.6f  J/mol%n", r.H);
        System.out.printf("  S (CVM)    = %+20.6f  J/(mol·K)%n", r.S);
        System.out.printf("  S (ideal)  = %+20.6f  J/(mol·K)   [expected]%n", sMixJmol);
        System.out.printf("  Ratio S/S_ideal = %+.4e   [should be ~1]%n", r.S / sMixJmol);
        System.out.printf("  G = H-TS   = %+20.6f  J/mol%n", r.G);
        System.out.printf("  ||∇G||     = %+.4e   [should be ~0 at true random state]%n", norm(r.Gu));
        sep();
    }

    private static void sep() {
        System.out.println("=======================================================");
    }

    private static double norm(double[] v) {
        double s = 0;
        for (double xi : v) s += Math.abs(xi);
        return s;
    }
}
