package org.ce.model.cvm;

import org.ce.model.cluster.CMatrixPipeline;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.cluster.StructurePhaseRegistry;
import org.ce.model.storage.InputLoader;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * Immutable product of the Stage 1-4 cluster identification pipeline: the
 * complete description of a lattice's cluster combinatorics in the CVCF
 * basis, for one {@code (elements, structure, model)} identity.
 *
 * <p><b>What is and is not in here.</b> This object holds everything fixed by
 * the lattice and the CVM approximation. It deliberately does <em>not</em>
 * hold the Hamiltonian ({@code CECEntry}) or the thermodynamic point
 * ({@code T, x, u}), which factorises the free energy as</p>
 *
 * <pre>
 *   G = f( CvmGeometry , CECEntry , T, x, u )
 *         built once     per Hamiltonian  per point
 * </pre>
 *
 * <p>Extracted from {@code CVMGibbsModel.initialize}, which previously
 * scattered these values across 17 mutable fields of the model. Two
 * consequences of making it a value:</p>
 *
 * <ul>
 *   <li>It is <b>independently testable</b>. {@link #validate()} checks the
 *       geometry's own invariants with no energy computation involved, so a
 *       bad cluster description is caught before any expression that consumes
 *       it is ever audited.</li>
 *   <li>It is <b>Hamiltonian-independent</b>, so it may be cached on
 *       {@code (elements, structure, model)} alone. Two Hamiltonians for the
 *       same system produce bit-identical geometry; a session-keyed cache
 *       pays the full Stage 1-4 cost for each.</li>
 * </ul>
 *
 * <p><b>Scope: cluster algebra only.</b> This class owns the lattice's cluster
 * data and the one linear map defined over it -- correlation functions to
 * cluster probabilities ({@link #evaluateCVs}) -- plus the self-consistency of
 * that data ({@link #validate()}). It deliberately stops there. Free energy,
 * entropy sums, ECI evaluation and any use of temperature belong to the model
 * layer above; a method here taking a temperature or a Hamiltonian would
 * dissolve the factorisation that makes this object independently testable.</p>
 *
 * <p>Array fields are exposed directly rather than defensively copied: this
 * object is shared read-only across every point of a scan, and copying
 * {@code cmat} per access would be wasteful. Treat every field as read-only.</p>
 *
 * <p><b>TODO -- ordered phases are not yet reviewed.</b> Everything here has
 * been built and verified against <em>disordered</em> parent phases only
 * (BCC_A2 at K=2/3/4). Ordered phases (BCC_B2 and friends) are known to differ
 * in ways this class currently only tolerates rather than handles
 * deliberately:</p>
 *
 * <ul>
 *   <li><b>The point set is larger than K.</b> {@code CvCfBasis} registers
 *       BCC_B2 binary with {@code numPointCfs=3} -- xA, xB and a long-range
 *       order parameter eta, all point-like and ECI-free, but eta is not a
 *       composition. So {@code tcf - ncf > K} there, and any code assuming the
 *       trailing block of {@code uFull} is exactly the mole fractions is wrong
 *       for ordered phases. {@link #buildFullVector} makes that assumption
 *       today.</li>
 *   <li><b>Point orbits split across sublattices.</b>
 *       {@code computeRandomCFs} uses {@code nxcf} rather than {@code K-1}
 *       point columns for this reason; at the random-state limit every
 *       sublattice shares one composition, so the columns coincide, but away
 *       from it they do not.</li>
 *   <li><b>{@link #parentStructure} exists but is unused.</b> It is stored so
 *       the pure-element reference can be taken on the parent phase, which is
 *       what {@code LatticeStability.g0}'s switch accepts; nothing reads it
 *       yet.</li>
 * </ul>
 *
 * <p>{@link #validate()} is written to accept {@code P >= K} so it will not
 * spuriously reject an ordered build, but the ordered path has not been
 * exercised. Revisit before relying on any ordered-phase result.</p>
 */
public final class CvmGeometry {

    // --- Identity ---------------------------------------------------------
    /** Hyphen-separated element string, e.g. {@code "Nb-Ti"}. */
    public final String elements;
    /** Elements in canonical order; size equals {@link #numComponents}. */
    public final List<String> elementList;
    /** Structure/phase this geometry was built for, e.g. {@code "BCC_A2"}. */
    public final String structure;
    /**
     * Disordered parent of {@link #structure} (equal to it for a parent phase).
     * Computed by {@code initialize} previously but discarded; retained here
     * because the pure-element reference {@code G0m} is defined on the parent
     * phase, and {@code LatticeStability.g0}'s switch only accepts parent
     * phase names.
     */
    public final String parentStructure;
    /** Number of components K. Drives the whole CF/CV combinatorics. */
    public final int numComponents;

    // --- Cluster combinatorics (Stage 1; the entropy prefactors) ----------
    /** Number of distinct disordered cluster types; outer {@code t} loop bound. */
    public final int tcdis;
    /** Kikuchi-Barker coefficients, {@code kb[t]}; length {@code tcdis}. */
    public final double[] kb;
    /** Disordered multiplicity per cluster type; length {@code tcdis}. */
    public final double[] mhdis;
    /** {@code mh[t][j]} = multiplicity of sub-orbit j within type t. */
    public final double[][] mh;
    /** {@code lc[t]} = number of symmetry-distinct sub-orbits in type t. */
    public final int[] lc;
    /** Full Stage 1-2 result, retained for diagnostics and SRO. */
    public final ClusterCFIdentificationPipeline.PipelineResult pipelineResult;

    // --- CVCF C-matrix (Stage 3-4; the u to cluster-probability map) ------
    /**
     * {@code cmat[t][j]} maps {@code uFull = [u ; x]} to cluster variables:
     * {@code cv[t][j][v] = sum_k cmat[t][j][v][k] * uFull[k]}. Rows =
     * {@code lcv[t][j]}, columns = {@link #tcf} (non-point CFs plus the full
     * point set), matching the width of {@code uFull}.
     *
     * <p>Note this differs from the <em>orthogonal</em> {@code CMatrixData}
     * produced at Stage 3, whose field comment documents {@code totalCFs + 1}
     * columns for a separate constant column. The Stage 4 CVCF matrices held
     * here carry no such column.</p>
     */
    public final List<List<double[][]>> cmat;
    /** {@code lcv[t][j]} = number of distinct CV rows; the {@code v} loop bound. */
    public final int[][] lcv;
    /** {@code wcv[t][j][v]} = how many raw configurations collapse into CV row v. */
    public final List<List<int[]>> wcv;
    /** Column-to-orthogonal-basis index map (Stage 3 provenance). */
    public final int[][] orthCfBasisIndices;

    // --- CVCF basis -------------------------------------------------------
    /** CVCF basis: CF names, ECI names, T/Tinv transforms, CF definitions. */
    public final CvCfBasis basis;
    /**
     * Number of <b>non-point</b> CFs: clusters larger than a point. These are
     * the only CFs that carry ECIs and the only ones that vary during
     * fixed-composition Newton-Raphson minimisation, so this is the length of
     * {@code u} and the solver's dimension.
     */
    public final int ncf;
    /**
     * Total CFs in the CVCF basis = {@code ncf} non-point CFs plus the full
     * point set.
     *
     * <p>Clusters come in three tiers: (i) the empty cluster, (ii) point
     * clusters, (iii) non-point clusters. Only tier (iii) has ECIs. The two
     * bases account for (i) and (ii) differently while reaching the same
     * total width:</p>
     *
     * <pre>
     *   orthogonal:  1 + (K-1) + ncf  =  ncf + K
     *   CVCF:            P    + ncf
     * </pre>
     *
     * <p>The CVCF basis carries no empty-cluster column -- its point CFs are
     * literal mole fractions whose sum is 1, so the constant is implied by the
     * composition constraint rather than stored. For a disordered phase the
     * point count P equals K and both totals agree at {@code ncf + K}; an
     * ordered phase can add point-like non-compositional CFs (BCC_B2 binary
     * registers {@code numPointCfs=3}: xA, xB and the LRO parameter eta), so
     * P &ge; K in general.</p>
     */
    public final int tcf;

    private CvmGeometry(
            String elements, List<String> elementList, String structure, String parentStructure,
            int numComponents, int tcdis, double[] kb, double[] mhdis, double[][] mh, int[] lc,
            ClusterCFIdentificationPipeline.PipelineResult pipelineResult,
            List<List<double[][]>> cmat, int[][] lcv, List<List<int[]>> wcv, int[][] orthCfBasisIndices,
            CvCfBasis basis, int ncf, int tcf) {
        this.elements = elements;
        this.elementList = elementList;
        this.structure = structure;
        this.parentStructure = parentStructure;
        this.numComponents = numComponents;
        this.tcdis = tcdis;
        this.kb = kb;
        this.mhdis = mhdis;
        this.mh = mh;
        this.lc = lc;
        this.pipelineResult = pipelineResult;
        this.cmat = cmat;
        this.lcv = lcv;
        this.wcv = wcv;
        this.orthCfBasisIndices = orthCfBasisIndices;
        this.basis = basis;
        this.ncf = ncf;
        this.tcf = tcf;
    }

    // =========================================================================
    // Construction
    // =========================================================================

    /**
     * Runs the full Stage 0-4 pipeline for one system identity. Expensive
     * (symmetry orbits, orthogonal C-matrix, CVCF transform) -- build once
     * and cache; see the class note on the Hamiltonian-independent cache key.
     *
     * <p>Body moved verbatim from {@code CVMGibbsModel.initialize}, minus the
     * {@code CECEntry} parameter, which that method only stored and logged --
     * across the whole pipeline it never reached
     * {@code ClusterCFIdentificationPipeline.run}, {@code CMatrixPipeline.run},
     * {@code verifyRandomCVs}, or {@code CvCfBasis.generate}, so it cannot
     * influence any field of the result.</p>
     *
     * @param elements     hyphen-separated element string, e.g. {@code "Nb-Ti"}
     * @param structure    structure/phase, e.g. {@code "BCC_A2"}
     * @param model        CVM approximation, e.g. {@code "T"}; a trailing
     *                     {@code "_CVCF"} is tolerated and stripped
     * @param progressSink nullable progress line consumer
     */
    public static CvmGeometry build(
            String elements,
            String structure,
            String model,
            Consumer<String> progressSink) {

        List<String> elementList = List.of(elements.split("-"));
        int numComponents = elementList.size();

        // --- INTERNAL PATH RESOLUTION ---
        String parentStructure = resolveParentStructure(structure);
        String disorderedFile = resolveClusterFile(parentStructure, model);
        String orderedFile = resolveClusterFile(structure, model);
        String disorderedSGName = resolveSymmetryGroup(parentStructure);
        String orderedSGName = resolveSymmetryGroup(structure);

        emit(progressSink, "\n  [CVM-Setup] Stage 0: Loading Inputs...");
        emit(progressSink, String.format("  > Elements:          %s", elements));
        emit(progressSink, String.format("  > Structure (Child): %s", structure));
        emit(progressSink, String.format("  > Structure (Parent):%s", parentStructure));
        emit(progressSink, String.format("  > Model:             %s", model));
        emit(progressSink, String.format("  > Components:        %d", numComponents));
        emit(progressSink, String.format("  > Files (Disord):    [clus: %s, sym: %s]",
                disorderedFile, disorderedSGName));
        emit(progressSink, String.format("  > Files (Order):     [clus: %s, sym: %s]",
                orderedFile, orderedSGName));

        // --- STAGE 0: LOADING ---
        List<Cluster> disorderedClusters = InputLoader.parseClusterFile(disorderedFile);
        disorderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup disorderedSG = InputLoader.parseSpaceGroup(disorderedSGName);

        List<Cluster> orderedClusters = InputLoader.parseClusterFile(orderedFile);
        orderedClusters.replaceAll(Cluster::sorted);
        SpaceGroup orderedSG = InputLoader.parseSpaceGroup(orderedSGName);

        // Extract transformation from the ordered phase mapping
        double[][] transformationMatrix = orderedSG.getRotateMat();
        double[] translationVector = orderedSG.getTranslateMat();

        // --- STAGE 1 & 2: IDENTIFICATION ---
        emit(progressSink, "  [CVM-Setup] Stage 1: Identification...");
        ClusterCFIdentificationPipeline.PipelineResult pr = ClusterCFIdentificationPipeline.run(
                disorderedClusters,
                disorderedSG.getOperations(),
                orderedClusters,
                orderedSG.getOperations(),
                transformationMatrix,
                translationVector,
                numComponents,
                progressSink);

        // --- STAGE 3: C-MATRIX (ORTHOGONAL) ---
        emit(progressSink, "  [CVM-Setup] Stage 2: C-Matrix Generation...");
        CMatrixPipeline.CMatrixData cmatOrth = CMatrixPipeline.run(
                pr.toClusterIdentificationResult(),
                pr.toCFIdentificationResult(),
                orderedClusters,
                numComponents,
                progressSink);

        // 3b. Build and Verify Random State
        double[] x = new double[numComponents];
        Arrays.fill(x, 1.0 / numComponents);
        emit(progressSink, "\n[Stage 3b] Testing CMatrixPipeline.verifyRandomCVs...");
        CMatrixPipeline.verifyRandomCVs(x, pr, cmatOrth, progressSink);

        // --- STAGE 4: CVCF BASIS & FINAL DATA ---
        emit(progressSink, "  [CVM-Setup] Stage 3: Basis Transformation...");
        CvCfBasis basisRef = CvCfBasis.generate(structure, pr, cmatOrth, model, progressSink);

        CMatrixPipeline.CMatrixData cmatCvcf = basisRef.cvcfCMatrixData;

        CvmGeometry geo = new CvmGeometry(
                elements, elementList, structure, parentStructure, numComponents,
                pr.getTcdis(), pr.getKbdis(), pr.getMhdis(), pr.getMh(), pr.getLc(), pr,
                cmatCvcf.getCmat(), cmatCvcf.getLcv(), cmatCvcf.getWcv(), cmatCvcf.getCfBasisIndices(),
                basisRef, basisRef.numNonPointCfs, basisRef.totalCfs());

        emit(progressSink, "  [CVM-Setup] ✓ Initialization complete.");
        return geo;
    }

    private static String resolveParentStructure(String structure) {
        return StructurePhaseRegistry.parentOf(structure);
    }

    private static String resolveClusterFile(String structure, String model) {
        String mod = model != null ? model.replace("_CVCF", "") : "";
        return "clus/" + structure + "-" + mod + ".txt";
    }

    private static String resolveSymmetryGroup(String structure) {
        return structure + "-SG";
    }

    private static void emit(Consumer<String> sink, String msg) {
        if (sink != null) sink.accept(msg);
    }

    // =========================================================================
    // Cluster algebra
    //
    // The C-matrix's whole content is one linear map, from correlation
    // functions to cluster probabilities. These two methods are that map;
    // everything else in this class is the data it is defined over.
    //
    // Nothing here takes a temperature or a Hamiltonian: cluster probabilities
    // are fixed by the lattice and the CF values alone. Energy, entropy and
    // ECI evaluation belong to the model layer above.
    // =========================================================================

    /**
     * Concatenates non-point CFs and composition into the full CVCF vector
     * {@code uFull = [u ; x]} this geometry's C-matrix multiplies against.
     *
     * <p>{@code ncf} is the split point: {@code u} occupies {@code [0, ncf)}
     * and {@code x} occupies {@code [ncf, ncf + K)}, for a total width of
     * {@link #tcf}. Note the CVCF vector carries <em>no</em> trailing constant
     * column, unlike the orthogonal-basis vector built by
     * {@code CMatrixPipeline.buildFullCFVector}, which appends {@code 1.0} for
     * the empty cluster and is therefore one element wider.</p>
     *
     * <p><b>Disordered phases only</b> -- see the class-level ordered-phase
     * TODO. This assumes the trailing block of {@code uFull} is exactly the K
     * mole fractions. An ordered phase's point set is wider than K (BCC_B2
     * binary adds an LRO parameter eta), so the trailing block is not
     * composition alone and this method cannot fill it from {@code x}. The
     * guard below rejects a mismatched {@code x} but does <em>not</em> detect
     * that case; for ordered phases assemble the vector yourself and use
     * {@link #evaluateCVsFull}.</p>
     *
     * @param u non-point CVCF correlation functions, length {@code >= ncf}
     * @param x mole fractions, length {@link #numComponents}
     */
    public double[] buildFullVector(double[] u, double[] x) {
        if (x.length != numComponents) {
            throw new IllegalArgumentException(
                    "buildFullVector: x.length=" + x.length + " != numComponents=" + numComponents);
        }
        if (tcf - ncf != numComponents) {
            throw new IllegalStateException(
                    "buildFullVector: this geometry has " + (tcf - ncf) + " point CFs but only K="
                            + numComponents + " mole fractions to fill them -- an ordered phase's"
                            + " point set includes non-compositional CFs (e.g. the BCC_B2 LRO"
                            + " parameter eta). Assemble uFull explicitly and call evaluateCVsFull."
                            + " See the ordered-phase TODO on this class.");
        }
        return CMatrixPipeline.buildFullCVCFVector(u, x, ncf);
    }

    /**
     * Evaluates every cluster variable (cluster probability) at the given
     * correlation functions and composition:
     * {@code cv[t][j][v] = sum_k cmat[t][j][v][k] * uFull[k]}.
     *
     * <p>Supplying the four geometry arrays ({@code cmat}, {@code lcv},
     * {@code tcdis}, {@code lc}) from one object is not merely convenience:
     * {@code CMatrixPipeline.evaluateCVs} guards at runtime against being
     * handed arrays from mismatched {@code PipelineResult}/{@code CMatrixData}
     * pairs, and routing through this method makes that mismatch
     * unrepresentable.</p>
     *
     * @return freshly allocated {@code cv[t][j][v]}; this object retains no
     *         reference to it and holds no state between calls
     */
    public double[][][] evaluateCVs(double[] u, double[] x) {
        return evaluateCVsFull(buildFullVector(u, x));
    }

    /**
     * Evaluates cluster variables from an already-assembled
     * {@code uFull = [u ; x]} vector -- for callers that build or perturb the
     * full vector themselves (the Hillert solver widens over it directly).
     *
     * @param uFull full CVCF vector, length {@link #tcf}
     */
    public double[][][] evaluateCVsFull(double[] uFull) {
        return CMatrixPipeline.evaluateCVs(uFull, cmat, lcv, tcdis, lc);
    }

    // =========================================================================
    // Validation
    // =========================================================================

    /**
     * Checks the geometry's own structural invariants -- shape agreement and
     * configuration-count conservation. Pure combinatorics: no Hamiltonian,
     * no temperature, no correlation functions, and therefore no physics.
     *
     * <p>These are worth checking separately from any energy because a fault
     * here corrupts the <em>inputs</em> to every free-energy expression. The
     * column-count check in particular guards the failure mode CLAUDE.md
     * records as having silently broken K&ge;3 CVCF energies (a point-CF
     * column-ordering regression) while every internal-consistency check kept
     * passing.</p>
     *
     * @throws IllegalStateException on the first invariant violated, naming
     *                               the offending indices and both values
     */
    public void validate() {
        // --- Outer dimension agreement ---
        require(cmat.size() == tcdis, "cmat.size()=" + cmat.size() + " != tcdis=" + tcdis);
        require(lc.length == tcdis, "lc.length=" + lc.length + " != tcdis=" + tcdis);
        require(lcv.length == tcdis, "lcv.length=" + lcv.length + " != tcdis=" + tcdis);
        require(wcv.size() == tcdis, "wcv.size()=" + wcv.size() + " != tcdis=" + tcdis);
        require(kb.length == tcdis, "kb.length=" + kb.length + " != tcdis=" + tcdis);
        require(mhdis.length == tcdis, "mhdis.length=" + mhdis.length + " != tcdis=" + tcdis);
        require(mh.length == tcdis, "mh.length=" + mh.length + " != tcdis=" + tcdis);
        require(elementList.size() == numComponents,
                "elementList.size()=" + elementList.size() + " != numComponents=" + numComponents);

        // --- CF tier accounting ---
        // Three tiers of cluster: (i) empty, (ii) point, (iii) non-point.
        // Only tier (iii) carries ECIs and varies during fixed-composition
        // minimisation, so it is exactly ncf wide.
        //
        // The two bases reach the same total by different bookkeeping:
        //
        //   orthogonal:  1 + (K-1) + ncf  =  ncf + K   (empty cluster explicit,
        //                                               K-1 independent points)
        //   CVCF:            P    + ncf              (no empty cluster, full
        //                                               point set instead)
        //
        // For a disordered phase P = K and both totals coincide at ncf + K.
        // Asserted as P >= K rather than P == K so an ordered phase's wider
        // point set is not spuriously rejected -- see the ordered-phase TODO
        // on this class; that path is unreviewed, not supported.
        int numPointCfs = tcf - ncf;
        require(numPointCfs >= numComponents,
                "tcf - ncf = " + numPointCfs + " point CFs, fewer than K=" + numComponents
                        + " -- the CVCF basis carries the full point set, not the K-1"
                        + " independent point CFs of the orthogonal basis");
        require(tcf == basis.totalCfs(),
                "tcf=" + tcf + " != basis.totalCfs()=" + basis.totalCfs());
        require(ncf == basis.numNonPointCfs,
                "ncf=" + ncf + " != basis.numNonPointCfs=" + basis.numNonPointCfs);

        for (int t = 0; t < tcdis; t++) {
            require(cmat.get(t).size() == lc[t],
                    "at t=" + t + ": cmat.get(t).size()=" + cmat.get(t).size() + " != lc[t]=" + lc[t]);
            require(lcv[t].length == lc[t],
                    "at t=" + t + ": lcv[t].length=" + lcv[t].length + " != lc[t]=" + lc[t]);
            require(wcv.get(t).size() == lc[t],
                    "at t=" + t + ": wcv.get(t).size()=" + wcv.get(t).size() + " != lc[t]=" + lc[t]);
            require(mh[t].length == lc[t],
                    "at t=" + t + ": mh[t].length=" + mh[t].length + " != lc[t]=" + lc[t]);

            int size = clusterSize(t);
            long expected = ipow(numComponents, size);

            for (int j = 0; j < lc[t]; j++) {
                double[][] block = cmat.get(t).get(j);
                int[] w = wcv.get(t).get(j);
                int nv = lcv[t][j];

                // --- C-matrix block shape ---
                require(block.length == nv,
                        "at (t=" + t + ", j=" + j + "): cmat rows=" + block.length
                                + " != lcv[t][j]=" + nv);
                require(w.length == nv,
                        "at (t=" + t + ", j=" + j + "): wcv length=" + w.length
                                + " != lcv[t][j]=" + nv);
                for (int v = 0; v < nv; v++) {
                    require(block[v].length == tcf,
                            "at (t=" + t + ", j=" + j + ", v=" + v + "): cmat columns="
                                    + block[v].length + " != tcf=" + tcf
                                    + " (= ncf non-point CFs + K point CFs, matching"
                                    + " the width of uFull = [u ; x])");
                }

                // --- Configuration-count conservation ---
                // The CV rows partition all K^(cluster size) configurations,
                // each row carrying its degeneracy, so the weights must sum
                // to that total exactly.
                long wsum = 0;
                for (int v = 0; v < nv; v++) {
                    require(w[v] > 0,
                            "at (t=" + t + ", j=" + j + ", v=" + v + "): wcv=" + w[v]
                                    + " must be positive");
                    wsum += w[v];
                }
                require(wsum == expected,
                        "at (t=" + t + ", j=" + j + "): sum(wcv)=" + wsum + " != K^size="
                                + numComponents + "^" + size + "=" + expected
                                + " -- the CV rows must partition every configuration");
            }
        }
    }

    /**
     * Number of sites in disordered cluster type {@code t}, taken from the
     * Stage 1 cluster data rather than inferred from array shapes, so the
     * configuration-count check in {@link #validate()} is an independent
     * assertion rather than a tautology.
     */
    private int clusterSize(int t) {
        return pipelineResult.getDisClusData().getClusCoordList().get(t).getAllSites().size();
    }

    private static long ipow(int base, int exp) {
        long r = 1;
        for (int i = 0; i < exp; i++) r *= base;
        return r;
    }

    private static void require(boolean condition, String message) {
        if (!condition) {
            throw new IllegalStateException("CvmGeometry invariant violated: " + message);
        }
    }

    @Override
    public String toString() {
        return String.format(
                "CvmGeometry[%s %s (parent %s) K=%d, tcdis=%d, ncf=%d, tcf=%d]",
                elements, structure, parentStructure, numComponents, tcdis, ncf, tcf);
    }
}
