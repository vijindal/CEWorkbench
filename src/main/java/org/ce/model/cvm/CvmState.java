package org.ce.model.cvm;

import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.List;

/**
 * Immutable evaluation of the CVM free energy at one thermodynamic point
 * {@code (T, x, u)}.
 *
 * <p>Replaces the {@code setT} / {@code setX} / {@code setU} sequence that
 * {@link CVMGibbsModel} uses to stage its mutable fields before each
 * {@code calculateXxx} call. That pattern has three problems this type
 * removes:</p>
 *
 * <ul>
 *   <li><b>Order dependence.</b> {@code setX} syncs the cluster variables only
 *       {@code if (u != null)} and {@code setU} only {@code if (x != null)}, so
 *       whichever is set first does not sync. It works today only because every
 *       call site happens to set both.</li>
 *   <li><b>Temporal coupling.</b> {@code calculateGm()} is meaningless unless
 *       the setters ran first, with the right values, and nothing intervened --
 *       a constraint the compiler cannot express. {@code checkMinimized()} is a
 *       partial runtime guard bolted on for exactly this reason.</li>
 *   <li><b>Shared mutable state.</b> Two callers cannot evaluate at different
 *       points concurrently, which is why {@code PhaseState} gives every phase
 *       its own {@code CVMGibbsModel}. An immutable state is safely shared.</li>
 * </ul>
 *
 * <p>Two derived quantities are computed eagerly in the constructor because
 * nearly every accessor needs them:</p>
 *
 * <pre>
 *   eci[l]        depends on T only        (a + b*T per Hamiltonian term)
 *   cv[t][j][v]   depends on (u, x) only   (via the CVCF C-matrix)
 * </pre>
 *
 * <p>Everything else is derived on demand. Gradients and Hessians are not
 * memoised: the {@code ncf x ncf} Hessian is not always wanted, and the
 * Newton loop asks for each exactly once per iteration.</p>
 *
 * <p><b>Both solvers consume this same object</b>, differing only in which
 * block they read. Fixed-composition Newton-Raphson treats {@code x} as a
 * constraint and reads the {@code ncf}-wide {@link #gmu()}/{@link #gmuu()};
 * the Hillert per-phase step treats composition as an unknown too and reads
 * the {@code (ncf+K)}-wide {@link #guFull()}/{@link #guuFull()}. Same physics,
 * different active set.</p>
 *
 * <p>Every expression here is ported verbatim from {@code CVMGibbsModel}'s
 * corresponding {@code calculateXxx} method, with field reads rebound to
 * {@link CvmGeometry}. No arithmetic was changed.</p>
 */
public final class CvmState {

    /**
     * Cluster variables below this value use a quadratic extension of
     * {@code cv*ln(cv)} rather than the logarithm itself, keeping the entropy
     * and its first two derivatives finite as a cluster probability approaches
     * zero. Matches {@code CVMGibbsModel.ENTROPY_SMOOTH_EPS}.
     */
    private static final double ENTROPY_SMOOTH_EPS = 1.0e-6;

    private final CvmGeometry geo;

    /** Temperature (K). */
    private final double temp;
    /** Mole fractions, length K. */
    private final double[] x;
    /** Non-point CVCF correlation functions, length ncf. */
    private final double[] u;

    /** Effective cluster interactions at this temperature, length ncf. */
    private final double[] eci;
    /** Cluster probabilities at this (u, x). */
    private final double[][][] cv;

    CvmState(CvmGeometry geo, CECEntry cecEntry, double temp, double[] x, double[] u) {
        if (x.length != geo.numComponents) {
            throw new IllegalArgumentException(
                    "x.length=" + x.length + " != numComponents=" + geo.numComponents);
        }
        if (u.length < geo.ncf) {
            throw new IllegalArgumentException(
                    "u.length=" + u.length + " < ncf=" + geo.ncf);
        }
        this.geo = geo;
        this.temp = temp;
        this.x = x.clone();
        this.u = u.clone();
        this.eci = CECEvaluator.evaluate(cecEntry, temp, geo.basis, "CVM", null);
        this.cv = geo.evaluateCVs(u, x);
    }

    // =========================================================================
    // The point itself
    // =========================================================================

    public double temperature() {
        return temp;
    }

    /** Mole fractions (defensive copy). */
    public double[] composition() {
        return x.clone();
    }

    /** Non-point CVCF correlation functions (defensive copy). */
    public double[] u() {
        return u.clone();
    }

    /** Effective cluster interactions at this temperature (defensive copy). */
    public double[] eci() {
        return eci.clone();
    }

    /** Cluster probabilities {@code cv[t][j][v]} at this point. */
    public double[][][] clusterVariables() {
        return cv;
    }

    /** The geometry this state was evaluated against. */
    public CvmGeometry geometry() {
        return geo;
    }

    // =========================================================================
    // Mixing quantities: Gm = Hm - T*Sm
    //
    // Gm is the CVM mixing contribution -- the ECI energy together with the
    // configurational entropy of mixing. It has no absolute reference zero of
    // its own; see g0m()/g() for the pure-element-anchored total.
    // =========================================================================

    /** Mixing enthalpy {@code Hm = sum_l eci[l]*u[l]}. */
    public double hm() {
        double h = 0.0;
        for (int l = 0; l < geo.ncf; l++) {
            h += eci[l] * u[l];
        }
        return h;
    }

    /** {@code dHm/du = eci} -- Hm is linear in u. */
    public double[] hmu() {
        return eci.clone();
    }

    /** {@code d2Hm/du2 = 0} -- Hm is linear in u. */
    public double[][] hmuu() {
        return new double[geo.ncf][geo.ncf];
    }

    /** Configurational entropy of mixing. */
    public double sm() {
        double sval = 0.0;
        for (int t = 0; t < geo.tcdis; t++) {
            double coeffT = geo.kb[t] * geo.mhdis[t];
            for (int j = 0; j < geo.lc[t]; j++) {
                double mhTj = geo.mh[t][j];
                int[] w = geo.wcv.get(t).get(j);
                int nv = geo.lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = cv[t][j][incv];
                    double sContrib;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        sContrib = cvVal * Math.log(cvVal);
                    } else {
                        double logEps = Math.log(ENTROPY_SMOOTH_EPS);
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        sContrib = ENTROPY_SMOOTH_EPS * logEps + (1.0 + logEps) * d
                                + 0.5 / ENTROPY_SMOOTH_EPS * d * d;
                    }
                    sval -= PhysicsConstants.R_GAS * coeffT * mhTj * w[incv] * sContrib;
                }
            }
        }
        return sval;
    }

    /** {@code dSm/du}, length {@code ncf}. */
    public double[] smu() {
        return entropyGradient(geo.ncf);
    }

    /** {@code d2Sm/du2}, {@code ncf x ncf}. */
    public double[][] smuu() {
        return entropyHessian(geo.ncf);
    }

    /** Mixing Gibbs energy {@code Gm = Hm - T*Sm}. */
    public double gm() {
        return hm() - temp * sm();
    }

    /** {@code dGm/du = dHm/du - T*dSm/du}, length {@code ncf}. */
    public double[] gmu() {
        double[] su = smu();
        double[] gu = new double[geo.ncf];
        for (int i = 0; i < geo.ncf; i++) {
            gu[i] = eci[i] - temp * su[i];
        }
        return gu;
    }

    /** {@code d2Gm/du2 = -T*d2Sm/du2}, since Hm is linear in u. */
    public double[][] gmuu() {
        double[][] suu = smuu();
        double[][] guu = new double[geo.ncf][geo.ncf];
        for (int i = 0; i < geo.ncf; i++) {
            for (int j = 0; j < geo.ncf; j++) {
                guu[i][j] = -temp * suu[i][j];
            }
        }
        return guu;
    }

    // =========================================================================
    // Widened quantities over uFull = [u ; x]
    //
    // For the Hillert per-phase step, where composition is an unknown rather
    // than a constraint. The leading ncf entries match the fixed-composition
    // gradient exactly; the trailing K come from the C-matrix's composition
    // columns.
    // =========================================================================

    /** Width of the widened parameter vector, {@code ncf + K}. */
    public int fullWidth() {
        return geo.ncf + geo.numComponents;
    }

    /** {@code dSm/d(uFull)}, length {@code ncf + K}. */
    public double[] smuFull() {
        return entropyGradient(fullWidth());
    }

    /** {@code d2Sm/d(uFull)2}, {@code (ncf+K) x (ncf+K)}. */
    public double[][] smuuFull() {
        return entropyHessian(fullWidth());
    }

    /**
     * {@code dGm/d(uFull)}. H contributes {@code eci} to the leading
     * {@code ncf} entries and zero to the trailing composition block, since
     * {@code Hm = sum_l eci[l]*u[l]} sums over {@code u} only.
     */
    public double[] gmuFull() {
        int width = fullWidth();
        double[] su = smuFull();
        double[] gu = new double[width];
        for (int i = 0; i < width; i++) {
            double hu = (i < geo.ncf) ? eci[i] : 0.0;
            gu[i] = hu - temp * su[i];
        }
        return gu;
    }

    /** {@code d2Gm/d(uFull)2 = -T * d2Sm/d(uFull)2}; H contributes zero everywhere. */
    public double[][] gmuuFull() {
        int width = fullWidth();
        double[][] suu = smuuFull();
        double[][] guu = new double[width][width];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                guu[i][j] = -temp * suu[i][j];
            }
        }
        return guu;
    }

    // =========================================================================
    // Reference energy and absolute totals: G = G0m + Gm
    // =========================================================================

    /**
     * Reference energy of the mechanical mixture of pure elements,
     * {@code G0m = sum_i x_i * G0(element_i, phase, T)}.
     *
     * <p>Pure energy: linear in composition, independent of {@code u}, and
     * carrying no configurational entropy. Hence {@code H0m = G0m},
     * {@code S0m = 0}, and every {@code u}-derivative of an absolute quantity
     * equals the mixing one.</p>
     */
    public double g0m() {
        return org.ce.model.equilibrium.LatticeStability.g0m(
                geo.elementList, geo.structure, x, temp);
    }

    /** {@code H0m = G0m} -- the reference term carries no entropy. */
    public double h0m() {
        return g0m();
    }

    /** {@code S0m = 0} -- a mechanical mixture of pure elements has no entropy of mixing. */
    public double s0m() {
        return 0.0;
    }

    /** Absolute Gibbs energy {@code G = G0m + Gm}. */
    public double g() {
        return g0m() + gm();
    }

    /** Absolute enthalpy {@code H = H0m + Hm}. */
    public double h() {
        return h0m() + hm();
    }

    /** Absolute entropy {@code S = S0m + Sm = Sm}. */
    public double s() {
        return sm();
    }

    /**
     * {@code dG/du = dGm/du}, since {@code G0m} does not depend on {@code u}.
     */
    public double[] gu() {
        return gmu();
    }

    /** {@code d2G/du2 = d2Gm/du2}. */
    public double[][] guu() {
        return gmuu();
    }

    /**
     * {@code dG/d(uFull)}: the mixing gradient plus {@code G0m}'s composition
     * block. {@code G0m} is linear in x with coefficients
     * {@code G0(element_i, phase, T)}, so no differentiation of the SGTE
     * polynomials is needed -- only evaluation.
     */
    public double[] guFull() {
        double[] g = gmuFull();
        for (int i = 0; i < geo.numComponents; i++) {
            g[geo.ncf + i] += org.ce.model.equilibrium.LatticeStability.g0(
                    geo.elementList.get(i), geo.structure, temp);
        }
        return g;
    }

    /**
     * {@code d2G/d(uFull)2 = d2Gm/d(uFull)2}. {@code G0m} is linear in x and
     * independent of u, so its second derivative vanishes identically.
     */
    public double[][] guuFull() {
        return gmuuFull();
    }

    // =========================================================================
    // Correlation functions
    // =========================================================================

    /**
     * Full CF vector {@code [u ; x]} at this point, length {@link #fullWidth()}.
     */
    public double[] cfs() {
        return geo.buildFullVector(u, x);
    }

    // =========================================================================
    // Shared entropy traversal
    //
    // calculateSmu/calculateSmuu and their Full counterparts were four separate
    // copies of one loop in CVMGibbsModel, differing only in the column bound.
    // They are parameterised on that bound here: the combinatorial prefactor
    // kb[t]*mhdis[t]*mh[t][j]*wcv[t][j][v] is then computed in one place, so a
    // wrong weight cannot appear in one copy and not the others.
    // =========================================================================

    private double[] entropyGradient(int width) {
        double[] su = new double[width];
        for (int t = 0; t < geo.tcdis; t++) {
            double coeffT = geo.kb[t] * geo.mhdis[t];
            for (int j = 0; j < geo.lc[t]; j++) {
                double mhTj = geo.mh[t][j];
                double[][] cm = geo.cmat.get(t).get(j);
                int[] w = geo.wcv.get(t).get(j);
                int nv = geo.lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = cv[t][j][incv];
                    double logEff;
                    if (cvVal > ENTROPY_SMOOTH_EPS) {
                        logEff = Math.log(cvVal);
                    } else {
                        double d = cvVal - ENTROPY_SMOOTH_EPS;
                        logEff = Math.log(ENTROPY_SMOOTH_EPS) + d / ENTROPY_SMOOTH_EPS;
                    }
                    double prefix = coeffT * mhTj * w[incv];
                    for (int l = 0; l < width; l++) {
                        double cml = cm[incv][l];
                        if (cml != 0.0) {
                            su[l] -= PhysicsConstants.R_GAS * prefix * cml * logEff;
                        }
                    }
                }
            }
        }
        return su;
    }

    private double[][] entropyHessian(int width) {
        double[][] suu = new double[width][width];
        for (int t = 0; t < geo.tcdis; t++) {
            double coeffT = geo.kb[t] * geo.mhdis[t];
            for (int j = 0; j < geo.lc[t]; j++) {
                double mhTj = geo.mh[t][j];
                double[][] cm = geo.cmat.get(t).get(j);
                int[] w = geo.wcv.get(t).get(j);
                int nv = geo.lcv[t][j];
                for (int incv = 0; incv < nv; incv++) {
                    double cvVal = cv[t][j][incv];
                    double invEff = (cvVal > ENTROPY_SMOOTH_EPS)
                            ? 1.0 / cvVal
                            : 1.0 / ENTROPY_SMOOTH_EPS;
                    double prefix = coeffT * mhTj * w[incv];
                    for (int l1 = 0; l1 < width; l1++) {
                        double cml1 = cm[incv][l1];
                        if (cml1 == 0.0) {
                            continue;
                        }
                        for (int l2 = l1; l2 < width; l2++) {
                            double cml2 = cm[incv][l2];
                            if (cml2 == 0.0) {
                                continue;
                            }
                            double val = -PhysicsConstants.R_GAS * prefix * cml1 * cml2 * invEff;
                            suu[l1][l2] += val;
                            if (l1 != l2) {
                                suu[l2][l1] += val;
                            }
                        }
                    }
                }
            }
        }
        return suu;
    }

    /**
     * True if every cluster variable lies strictly inside {@code (0, 1)} --
     * the physicality check the solvers use before accepting a trial step.
     */
    public boolean isValid() {
        for (int t = 0; t < geo.tcdis; t++) {
            for (int j = 0; j < geo.lc[t]; j++) {
                for (int v = 0; v < geo.lcv[t][j]; v++) {
                    double value = cv[t][j][v];
                    if (!(value > 0.0 && value < 1.0)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    @Override
    public String toString() {
        return String.format("CvmState[T=%.2f, x=%s, ncf=%d]",
                temp, java.util.Arrays.toString(x), geo.ncf);
    }

    /** Element list of the geometry, for callers building reference terms. */
    List<String> elementList() {
        return geo.elementList;
    }
}
