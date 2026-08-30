package org.ce.model.cvm;

import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * The CVM Gibbs model: a <b>pure evaluator</b> mapping
 *
 * <pre>
 *   system parameters  (elements, structure, ECIs)   bound once, here
 *   macro parameters   (T, x)                        per evaluation
 *   micro parameters   (u)                           per evaluation
 *       --&gt;  G, Gm, Hm, Sm, G0m and their derivatives
 * </pre>
 *
 * <p>It holds no per-point state: no current {@code (T, x, u)}, no
 * {@code isMinimized} flag, no setter sequence to get wrong. Every quantity is
 * reached through {@link #at}, which returns an immutable {@link State}
 * carrying everything evaluated at the point you asked for. Two threads may
 * call it concurrently, and a caller may hold states at several points at
 * once.</p>
 *
 * <p><b>Solvers are not part of this class.</b> They hold an instance and drive
 * it from the outside, each owning its own loop:</p>
 *
 * <ul>
 *   <li>{@code CvmNewtonSolver} -- fixed-composition minimisation, where
 *       {@code x} is a constraint and {@code u} the unknown. Reads the
 *       {@code ncf}-wide {@link State#gmu()} and {@link State#gmuu()}.</li>
 *   <li>{@code HillertSolver.PhaseStep} -- the multi-phase per-phase step, where
 *       composition is an unknown too. Reads the {@code (ncf+K)}-wide
 *       {@link State#gmuFull()} and {@link State#gmuuFull()}.</li>
 * </ul>
 *
 * <p>Same physics, same evaluation, different active set -- which is why one
 * model serves both. Both live in {@code org.ce.model.equilibrium}, alongside
 * the multi-phase machinery they serve; this package holds the model they
 * evaluate against.</p>
 *
 * <p>This class previously carried the Stage 1-4 pipeline, seventeen fields of
 * geometry, a mutable point, five duplicated copies of the entropy expression,
 * and the Newton-Raphson loop itself. Setup now lives in {@link CvmGeometry}
 * and the solvers in their own classes; what remains is evaluation, with each
 * expression defined exactly once.</p>
 */
public final class CVMGibbsModel {

    private final CvmGeometry geo;
    private final CECEntry cecEntry;

    /**
     * @param geometry the lattice's cluster algebra, from {@link CvmGeometry#build}
     * @param cecEntry the Hamiltonian supplying ECIs, matched by name against
     *                 {@code geometry.basis} at each evaluation
     */
    public CVMGibbsModel(CvmGeometry geometry, CECEntry cecEntry) {
        if (geometry == null) {
            throw new IllegalArgumentException("geometry must not be null");
        }
        this.geo = geometry;
        this.cecEntry = cecEntry;
    }

    /**
     * Builds a model for one system identity, running the Stage 1-4 pipeline.
     *
     * <p>Expensive. The geometry it produces is Hamiltonian-independent, so a
     * caller evaluating several Hamiltonians on the same lattice should build
     * one {@link CvmGeometry} and share it across constructor calls rather than
     * repeating this.</p>
     */
    public static CVMGibbsModel of(
            String elements, String structure, String model,
            CECEntry cecEntry, Consumer<String> progressSink) {

        if (progressSink != null && cecEntry != null) {
            progressSink.accept(String.format("  > CEC Entry:         %s (%s)",
                    cecEntry.elements, cecEntry.structurePhase));
            if (cecEntry.cecTerms != null) {
                for (CECEntry.CECTerm term : cecEntry.cecTerms) {
                    progressSink.accept(String.format("    - %-10s: a = %10.6f, b = %10.6f",
                            term.name, term.a, term.b));
                }
            }
        }

        CvmGeometry geometry = CvmGeometry.build(elements, structure, model, progressSink);
        geometry.validate();
        return new CVMGibbsModel(geometry, cecEntry);
    }

    // =========================================================================
    // Evaluation
    // =========================================================================

    /**
     * Evaluates at one thermodynamic point.
     *
     * <p>The returned {@link State} belongs to this model and can name it
     * again through {@link State#model()}. That matters when a state outlives
     * the call: a solver evaluates at its converged iterate and returns the
     * state inside its result, so what a caller receives from
     * {@code solver.solve(...).state()} is exactly what this method produced --
     * already evaluated, and needing no repeat.</p>
     *
     * @param temperature temperature in K
     * @param x           mole fractions, length {@code numComponents}
     * @param u           non-point CVCF correlation functions, length {@code >= ncf}
     */
    public State at(double temperature, double[] x, double[] u) {
        return new State(temperature, x, u);
    }

    /**
     * Evaluates from a joint {@code uFull = [u ; x]} vector -- the form the
     * Hillert solver carries, where composition is part of the unknown rather
     * than a separate input.
     */
    public State atFull(double temperature, double[] uFull) {
        int width = geo.ncf + geo.numComponents;
        if (uFull.length != width) {
            throw new IllegalArgumentException(
                    "uFull.length=" + uFull.length + " != ncf+K=" + width);
        }
        double[] u = new double[geo.ncf];
        double[] x = new double[geo.numComponents];
        System.arraycopy(uFull, 0, u, 0, geo.ncf);
        System.arraycopy(uFull, geo.ncf, x, 0, geo.numComponents);
        return at(temperature, x, u);
    }

    // =========================================================================
    // System parameters
    // =========================================================================

    /** The lattice cluster algebra this model evaluates against. */
    public CvmGeometry geometry() {
        return geo;
    }

    /** The Hamiltonian this model evaluates against. */
    public CECEntry cecEntry() {
        return cecEntry;
    }

    /** Number of non-point CFs: the length of {@code u}, and a solver's dimension. */
    public int ncf() {
        return geo.ncf;
    }

    /** Number of components K. */
    public int numComponents() {
        return geo.numComponents;
    }

    /**
     * Correlation functions of the fully disordered state at this composition
     * -- a natural starting iterate for a minimisation. Returns the leading
     * {@code ncf} block.
     */
    public double[] randomStateU(double[] x) {
        double[] full = geo.basis.computeRandomCvcfCFs(x, geo.pipelineResult);
        double[] u = new double[geo.ncf];
        System.arraycopy(full, 0, u, 0, geo.ncf);
        return u;
    }

    /** Full random-state CVCF vector {@code [u ; x]} at this composition. */
    public double[] randomStateFull(double[] x) {
        return geo.basis.computeRandomCvcfCFs(x, geo.pipelineResult);
    }

    /**
     * Cluster probabilities at an arbitrary {@code (u, x)}, without building a
     * full {@link State} -- for callers that want only the cluster variables,
     * such as SRO post-processing and a solver's step limiter.
     */
    public double[][][] clusterVariablesAt(double[] u, double[] x) {
        return geo.evaluateCVs(u, x);
    }

    @Override
    public String toString() {
        return "CVMGibbsModel[" + geo + "]";
    }

    // =========================================================================
    // State: the model evaluated at one point
    // =========================================================================

    /**
     * Immutable evaluation of the CVM free energy at one thermodynamic point
     * {@code (T, x, u)}.
     *
     * <p>Nested inside the model because a state is meaningless without the
     * geometry and Hamiltonian it was evaluated against -- it is the model at a
     * point, not an independent object.</p>
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
     * memoised: the {@code ncf x ncf} Hessian is not always wanted, and a
     * Newton loop asks for each exactly once per iteration.</p>
     */
    /**
     * A neighbour shell whose pair cluster variables carry SRO information.
     *
     * <p>The cluster-type indices are a property of the BCC_A2 tetrahedron
     * layout, verified empirically: types run {@code t=0} tetrahedron,
     * {@code t=1} triangle, {@code t=2} 2nd-NN pair, {@code t=3} 1st-NN pair,
     * {@code t=4} point. Note the naming trap -- the CVCF term {@code v21}
     * (1st-NN) is {@code t=3}, while {@code v22} (2nd-NN) is {@code t=2}.</p>
     *
     * <p>Pair cluster variables are laid out in upper-triangular species order,
     * confirmed by evaluating the random state at a skewed composition, where
     * each entry reproduces {@code x_P * x_R} exactly.</p>
     */
    public enum Shell {
        /** First nearest neighbour; CVCF name {@code v21...}. */
        FIRST(3, "1NN"),
        /** Second nearest neighbour; CVCF name {@code v22...}. */
        SECOND(2, "2NN");

        private final int clusterType;
        private final String label;

        Shell(int clusterType, String label) {
            this.clusterType = clusterType;
            this.label = label;
        }

        /** Index of this shell's pair cluster among the disordered cluster types. */
        public int clusterType() {
            return clusterType;
        }

        /** Short name used when reporting, e.g. {@code "1NN"}. */
        public String label() {
            return label;
        }

        @Override
        public String toString() {
            return label;
        }
    }

    /** One Cowley-Warren pair parameter. */
    public static final class PairSro {
        /** Index of the first species in the system's canonical element order. */
        public final int i;
        /** Index of the second species. */
        public final int j;
        /** Cluster probability p_ij^PR from the CVM. */
        public final double probability;
        /** Random-state reference x_P*x_R. */
        public final double reference;
        /** alpha = 1 - p/(x_P x_R). */
        public final double alpha;

        public PairSro(int i, int j, double probability, double reference, double alpha) {
            this.i = i;
            this.j = j;
            this.probability = probability;
            this.reference = reference;
            this.alpha = alpha;
        }

        @Override
        public String toString() {
            return String.format("PairSro[%d-%d p=%.6f ref=%.6f alpha=%+.6f]",
                    i, j, probability, reference, alpha);
        }
    }

    /**
     * One Cowley-Warren quaternary tetrahedron parameter (Jindal &amp; Lele
     * 2025, Eq. 41), {@code alpha^MPRT = 1 - rho^MPRT/(x_M x_P x_R x_T)}, read
     * off a single {@code product(...)}-type CVCF CF such as {@code v4ABCD1}.
     */
    public static final class TetrahedronSro {
        /** The CVCF CF name this was read from, e.g. {@code "v4ABCD1"}. */
        public final String cfName;
        /** Cluster probability rho^MPRT from the CVM. */
        public final double probability;
        /** Random-state reference x_M*x_P*x_R*x_T. */
        public final double reference;
        /** alpha = 1 - rho/(x_M x_P x_R x_T). */
        public final double alpha;

        public TetrahedronSro(String cfName, double probability, double reference, double alpha) {
            this.cfName = cfName;
            this.probability = probability;
            this.reference = reference;
            this.alpha = alpha;
        }

        @Override
        public String toString() {
            return String.format("TetrahedronSro[%s p=%.6f ref=%.6f alpha=%+.6f]",
                    cfName, probability, reference, alpha);
        }
    }

    public final class State {


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

        State(double temp, double[] x, double[] u) {
            if (x.length != geo.numComponents) {
                throw new IllegalArgumentException(
                        "x.length=" + x.length + " != numComponents=" + geo.numComponents);
            }
            if (u.length < geo.ncf) {
                throw new IllegalArgumentException(
                        "u.length=" + u.length + " < ncf=" + geo.ncf);
            }
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

        /**
         * The model that produced this state.
         *
         * <p>A state is always an evaluation of some model at some point --
         * {@code model.at(T, x, u)} -- and this makes that link visible rather
         * than leaving it implied by the nested type name. It matters where a
         * state has travelled away from its origin: a solver returns one inside
         * its result, and a caller holding several models otherwise has no way
         * to tell which one an arriving state belongs to.</p>
         */
        public CVMGibbsModel model() {
            return CVMGibbsModel.this;
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
                        // cv*ln(cv) is well behaved in double precision down to
                        // ~5e-324, far below any physical cluster probability, so
                        // no regularisation is needed. Only cv == 0 is special:
                        // its limit is 0, but 0*ln(0) evaluates to 0*(-Inf) = NaN.
                        double sContrib = cvVal > 0.0 ? cvVal * Math.log(cvVal) : 0.0;
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

        // =====================================================================
        // Short-range order
        //
        // Cowley-Warren SRO is a thermodynamic property of this point, on the
        // same footing as G/H/S: it is read off the very cluster probabilities
        // the free energy is built from. Jindal & Lele, Calphad 89 (2025)
        // 102825, make the point explicitly -- the CVM already produces these
        // probabilities during minimisation, so no extra machinery is needed:
        //
        //   pair:      alpha_ij^PR   = 1 - p_ij^PR / (x_P x_R)          (Eq. 40)
        //   multisite: alpha^MPRT    = 1 - rho^MPRT / (x_M x_P x_R x_T) (Eq. 41)
        //
        // Sign convention:
        //   alpha < 0   unlike pairs enriched relative to random -- ordering
        //   alpha = 0   random (ideal) solution
        //   alpha > 0   unlike pairs depleted -- clustering / phase separation
        //
        // Only pair SRO is exposed. Extending to triangle/tetrahedron clusters
        // looks straightforward -- their probabilities are equally available --
        // but most CVCF correlation functions for those clusters are not single
        // physical probabilities: CvCfBasis.VSpec defines each CF as either
        // product(...) (one site-atom-pair probability, SRO-eligible) or
        // diff(...)/combo(...) (a signed combination, e.g. v3AB = rho^RPR -
        // rho^PRP). A diff-type CF has no natural [0,1] reference, so
        // 1 - value/reference is not meaningful for it. Adding multisite SRO
        // means splitting those into their constituent product(...) terms
        // first, not computing from the CF as defined.
        // =====================================================================

        /**
         * Cowley-Warren pair SRO parameters for one neighbour shell at this
         * point, in upper-triangular species order -- for K=3 that is
         * {@code [AA, AB, AC, BB, BC, CC]}.
         *
         * <p>Like pairs ({@code i == j}) are included; they carry the
         * complementary information to the unlike ones. For the unlike pairs
         * alone -- the ones normally reported and plotted -- use
         * {@link #unlikePairSro}.</p>
         *
         * @param shell which neighbour shell, {@link Shell#FIRST} or
         *              {@link Shell#SECOND}
         * @throws IllegalArgumentException if the cluster-variable block does
         *         not have the {@code K(K+1)/2} entries a pair cluster must
         *         have, which would mean the layout assumption no longer holds
         *         and the species mapping cannot be trusted
         */
        public List<PairSro> pairSro(Shell shell) {
            int K = geo.numComponents;
            int expected = K * (K + 1) / 2;
            int t = shell.clusterType();

            if (t < 0 || t >= cv.length || cv[t].length == 0) {
                throw new IllegalArgumentException(
                        "No cluster variables for " + shell + " (cluster type t=" + t + ")");
            }

            double[] p = cv[t][0];
            if (p.length != expected) {
                throw new IllegalArgumentException(
                        shell + " (cluster type t=" + t + ") has " + p.length
                                + " cluster variables, expected " + expected + " for K=" + K
                                + " (upper-triangular species pairs). The cluster-variable"
                                + " layout assumed for pair SRO does not hold for this system.");
            }

            List<PairSro> out = new ArrayList<>(expected);
            int v = 0;
            for (int i = 0; i < K; i++) {
                for (int j = i; j < K; j++, v++) {
                    double ref = x[i] * x[j];
                    double alpha = (ref > 0.0) ? 1.0 - p[v] / ref : Double.NaN;
                    out.add(new PairSro(i, j, p[v], ref, alpha));
                }
            }
            return out;
        }

        /** Pair SRO for the unlike pairs ({@code i != j}) only. */
        public List<PairSro> unlikePairSro(Shell shell) {
            List<PairSro> out = new ArrayList<>();
            for (PairSro s : pairSro(shell)) {
                if (s.i != s.j) out.add(s);
            }
            return out;
        }

        /**
         * Pair SRO for every shell this geometry supports, keyed by shell name
         * -- the form the calculation layer reports.
         *
         * <p>A shell whose cluster-variable layout is not the expected
         * upper-triangular pair block is omitted rather than failing: SRO is
         * supplementary, and an unsupported geometry must not take the rest of
         * a calculation down with it.</p>
         */
        public Map<String, List<PairSro>> pairSroByShell() {
            Map<String, List<PairSro>> out = new LinkedHashMap<>();
            for (Shell shell : Shell.values()) {
                try {
                    out.put(shell.label(), pairSro(shell));
                } catch (IllegalArgumentException e) {
                    // Not a pair block in this geometry; skip this shell.
                }
            }
            return out;
        }

        /**
         * Cowley-Warren quaternary tetrahedron SRO (Eq. 41), read directly off
         * the named CVCF CF -- {@code v4ABCD1} for the standard 4-component
         * BCC_A2/T basis. Unlike a pair cluster, the quaternary tetrahedron
         * has no upper-triangular species enumeration to walk: {@code v4ABCD1}
         * already names one specific site-to-species assignment and is
         * exactly {@code rho^ABCD}, a single {@code product(...)}-type
         * probability (see {@link CvCfBasis.VSpec}) with a well-defined
         * {@code x_A x_B x_C x_D} reference -- SRO-eligible by the same rule
         * {@link #pairSro} uses.
         *
         * @param cfName the CVCF name of a quaternary tetrahedron CF, e.g.
         *               {@code "v4ABCD1"}
         * @throws IllegalArgumentException if this geometry is not a
         *         4-component system, or {@code cfName} does not resolve in
         *         its basis
         */
        public TetrahedronSro quaternaryTetrahedronSro(String cfName) {
            if (geo.numComponents != 4) {
                throw new IllegalArgumentException(
                        "Quaternary tetrahedron SRO requires numComponents=4, got " + geo.numComponents);
            }
            int idx = geo.basis.indexOfCf(cfName);
            if (idx < 0) {
                throw new IllegalArgumentException(
                        "'" + cfName + "' is not a CF in this basis.");
            }
            double rho = cfs()[idx];
            double ref = x[0] * x[1] * x[2] * x[3];
            double alpha = (ref > 0.0) ? 1.0 - rho / ref : Double.NaN;
            return new TetrahedronSro(cfName, rho, ref, alpha);
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
                        // No regularisation: see calculateSm. A non-positive cv
                        // gives -Infinity/NaN here, which is intended -- both
                        // solvers reject such a state before it reaches this
                        // point (CvmNewtonSolver step 3, HillertSolver's
                        // backtracking validity check), so a NaN escaping into
                        // a gradient means a solver skipped its own guard.
                        double logEff = Math.log(cvVal);
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
                        // 1/cv is the true curvature of cv*ln(cv); it is large
                        // at a dilute composition because the curvature really
                        // is large there, not because of a numerical artifact.
                        double invEff = 1.0 / cvVal;
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
         * True if every cluster variable lies strictly inside {@code (0, 1)},
         * across <em>all</em> cluster types <b>including the point block</b> --
         * the physicality check the Hillert solver applies before accepting a
         * trial step.
         *
         * <p>Deliberately broader than the non-point-only minimum a
         * fixed-composition minimisation checks ({@code CvmNewtonSolver}
         * excludes the point type, following the reference's {@code findMin}).
         * Including the point block also catches a trial point whose
         * composition itself has drifted to a pure-element boundary
         * ({@code x_i = 0} or {@code 1}), which is exactly the case the Hillert
         * backtracking exists to catch. Do not narrow this to the non-point
         * types; that would silently reduce what it detects.</p>
         */
        public boolean isValidIncludingPoints() {
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
            return String.format("State[T=%.2f, x=%s, ncf=%d]",
                    temp, java.util.Arrays.toString(x), geo.ncf);
        }


    }
}
