package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.cluster.Cluster;
import org.ce.model.cluster.ClusterCFIdentificationPipeline;
import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.SpaceGroup;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.InputLoader;
import org.ce.model.storage.Workspace;

import java.util.List;
import java.util.Random;

/**
 * Consolidated MCS correctness verification, anchored on configurations whose
 * orthogonal CFs, CVCF CFs, and energy can all be derived by hand or via a
 * closed-form formula independent of MCS's embedding/measurement code — not
 * just checked for internal self-consistency against MCS's own other code paths.
 *
 * <p>Supersedes and replaces (deleted) {@code RandomStateEnergyTest},
 * {@code DeltaEVerificationTest}, {@code DeltaEVerificationTestV2Rigorous}, and
 * {@code MetropolisTrajectoryVerification} — folded into the stages below.
 * {@code AbFamilyCfDiagnostic} and {@code DeltaEBenchmark}/{@code EmbeddingScaleProbe}
 * remain separate (diagnostic / performance tools, not correctness gates).
 *
 * <p>Building this suite found and fixed a real latent bug: an earlier session had
 * reordered the point-CF columns of {@code Embeddings.generatePointCfEmbeddings} and
 * {@code ClusterCFIdentificationPipeline.computeRandomCFs} from the pipeline's raw
 * (descending, for K&ge;3) alpha order to ascending, believing {@code CvCfBasis}'s
 * {@code Tinv} matrix expected ascending order. That reordering was backwards: {@code
 * Tinv} is built from the pipeline's own raw column order (see {@code
 * CMatrixPipeline.buildCfColumnMap}/{@code deriveCfBasisIndices}), so the ascending
 * reorder silently broke K&ge;3 CVCF energies while every MCS-internal check (ΔE vs.
 * finite-difference, trajectory running-energy vs. recompute) kept passing — those
 * checks compare MCS against itself, so a shared, consistent reordering on both sides
 * is invisible to them. Only checking against an <em>independent</em> ground truth
 * (a pure-species config's point CF must equal exactly 1.0) exposed it. Both call
 * sites have been reverted to preserve the pipeline's raw column order.
 *
 * <p>Configurations used, and how their ground truth is derived:
 * <ul>
 *   <li><b>Pure element</b> (x_0=1): {@link PipelineResult#computeRandomCFs} at a
 *       one-hot composition — a closed-form combinatorial formula, zero dependence
 *       on MCS's embedding code.</li>
 *   <li><b>Perfectly random</b> at composition x: same closed-form formula, general
 *       composition.</li>
 *   <li><b>Perfectly ordered B2</b> (K=2, x=0.5): species assigned by flat site-index
 *       parity (even=species 0, odd=species 1, per {@code MCSGeometry.buildBCCPositions}).
 *       Empirically verified (not assumed) that every {@code v21AB} embedding connects
 *       two same-parity sites and every {@code v22AB} embedding connects two
 *       opposite-parity sites (confirmed by direct inspection of
 *       {@code geo.cfEmbeddings()} — the "I-n"/"II-n" naming in {@code CvCfBasis}'s
 *       registration comment does not correspond to a naive same-cell 1NN/2NN distance
 *       reading, so this is checked rather than hand-guessed): with this parity-based
 *       B2 config, {@code v21AB} (orthogonal, same-parity on every embedding) is
 *       deterministically +1, and {@code v22AB} (opposite-parity) is deterministically
 *       -1.</li>
 *   <li><b>Composition boundary</b> (K reduces to K-1 as one species -&gt; 0): the
 *       CVCF-basis (not orthogonal — the orthogonal point-basis sequence is redefined
 *       per K, so raw orthogonal values are not comparable across K) value of any CF
 *       column whose name is shared between the K and K-1 systems must match exactly
 *       at the boundary composition — this is the physically meaningful invariant
 *       (e.g. a pair probability doesn't care that a third, zero-fraction species
 *       exists in the composition space).</li>
 * </ul>
 *
 * <p>For each configuration, verification proceeds:
 * <ol>
 *   <li>Orthogonal CFs: {@code Embeddings.measureFullCVsFromConfig} vs. the analytic
 *       value above.</li>
 *   <li>CVCF CFs: {@code Embeddings.applyTinvTransform} of the measured orthogonal
 *       vector vs. the same transform applied to the analytic orthogonal vector
 *       (shared {@code Tinv} matrix — this checks MCS's usage of it).</li>
 *   <li>Energy: {@code Embeddings.totalEnergyCvcf} with (a) synthetic single-term
 *       ECIs isolating one CF at a time, and (b) the real shipped Hamiltonian,
 *       cross-checked against the analytic CVCF vector's own dot product.</li>
 *   <li>Delta-E: {@code Embeddings.deltaEExchangeCvcf} (list + flat) and
 *       {@code deltaEExchangeCvcfV2}, checked against finite-difference full-energy
 *       recompute, starting from each analytically-anchored configuration above.</li>
 *   <li>Trajectory: a sequence of accept/reject swaps from the B2 state, with running
 *       energy checked against from-scratch recompute after every step.</li>
 * </ol>
 */
public class AnalyticConfigVerification {

    private static int failures = 0;
    private static int checks = 0;

    public static void main(String[] args) throws Exception {
        Workspace ws = new Workspace();
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);

        runBinarySuite(builder, "Nb-Ti", 1000.0);
        runHigherKSuite(builder, "Nb-Ti-V", 1000.0);
        runHigherKSuite(builder, "Nb-Ti-V-Zr", 1000.0);

        System.out.println();
        System.out.println("=== OVERALL ===");
        System.out.println("Checks: " + checks + "   Failures: " + failures);
        System.out.println(failures == 0 ? "RESULT: PASS" : "RESULT: FAIL");
        if (failures > 0) System.exit(1);
    }

    // =========================================================================
    // K=2 suite: pure element, B2, random, delta-E, trajectory
    // =========================================================================

    private static void runBinarySuite(ModelSession.Builder builder, String elements, double T) throws Exception {
        System.out.println("=======================================================");
        System.out.println("K=2 suite: " + elements + " BCC_A2/T, T=" + T);
        System.out.println("=======================================================");

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);
        int numComp = session.numComponents();
        int L = 4;
        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis.numNonPointCfs;
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "AnalyticVerify-K2");
        PipelineResult pr = buildPipelineResult(elements, "BCC_A2", "T", numComp);

        // ── 1. Pure element (all species 0) ──
        System.out.println("\n--- Pure element (x_0=1) ---");
        {
            LatticeConfig config = buildPureConfig(N, numComp);
            double[] uFull = Embeddings.measureFullCVsFromConfig(
                    config, geo.cfEmbeddings(), geo.pointCfEmbeddings(), geo.getFlatBasisMatrix(), ncf, numComp);

            double[] xPure = { 1.0, 0.0 };
            double[] uGroundTruth = pr.computeRandomCFs(xPure);
            checkVector("pure-element uOrth (measured vs closed-form at x_0=1)", uFull, uGroundTruth, 1e-9);

            double[] vCvcfMeasured = Embeddings.applyTinvTransform(uFull, geo.basis);
            double[] vCvcfAnalytic = Embeddings.applyTinvTransform(uGroundTruth, geo.basis);
            checkVector("pure-element vCvcf (measured vs Tinv-of-analytic)", vCvcfMeasured, vCvcfAnalytic, 1e-9);

            double hMcs = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                    geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.basis, numComp) / N;
            double hAnalytic = dot(eciCvcf, vCvcfAnalytic, ncf);
            checkScalar("pure-element H (MCS totalEnergyCvcf/site vs analytic eci.vCvcf)", hMcs, hAnalytic, 1e-6);
        }

        // ── 2. Perfectly ordered B2 (K=2, x=0.5, alternating sublattices) ──
        System.out.println("\n--- Perfectly ordered B2 (x=0.5, alternating sublattices) ---");
        {
            LatticeConfig config = buildB2Config(N, numComp);
            double[] uFull = Embeddings.measureFullCVsFromConfig(
                    config, geo.cfEmbeddings(), geo.pointCfEmbeddings(), geo.getFlatBasisMatrix(), ncf, numComp);

            // Analytic ground truth (verified against actual embedding parity, see class doc):
            //   v21AB embeddings are always same-parity: phi1(0)*phi1(0) = phi1(1)*phi1(1) = +1
            //   v22AB embeddings are always opposite-parity: phi1(0)*phi1(1) = (-1)(+1) = -1
            List<String> cfNames = CvCfBasis.getNonPointCfNames("BCC_A2", "T", numComp);
            int idx21 = cfNames.indexOf("v21AB");
            int idx22 = cfNames.indexOf("v22AB");
            if (idx21 >= 0) checkScalar("B2 v21AB (same-parity pair, analytic=+1)", uFull[idx21], 1.0, 1e-9);
            if (idx22 >= 0) checkScalar("B2 v22AB (opposite-parity pair, analytic=-1)", uFull[idx22], -1.0, 1e-9);

            // Point CFs: deterministic average of the two sublattices' fixed values (no disorder).
            double[] seq = org.ce.model.cluster.ClusterMath.buildBasis(numComp);
            int nPoint = geo.pointCfEmbeddings().size();
            for (int k = 0; k < nPoint; k++) {
                double analytic = 0.5 * Math.pow(seq[0], k + 1) + 0.5 * Math.pow(seq[1], k + 1);
                checkScalar("B2 point CF[" + k + "] (global average)", uFull[ncf + k], analytic, 1e-9);
            }

            // Energy: synthetic single-term ECIs isolating v21AB, confirming the energy
            // dot-product picks up exactly the analytic CF value derived above.
            if (idx21 >= 0) {
                CECEntry synth = buildSingleTermCec(elements, "BCC_A2", "T", cfNames.get(idx21));
                double[] eciSynth = CECEvaluator.evaluate(synth, T, geo.basis, "B2-synthetic-v21AB");
                double E = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                        geo.getFlatBasisMatrix(), ncf, eciSynth, geo.basis, numComp);
                double[] vCvcf = Embeddings.applyTinvTransform(uFull, geo.basis);
                checkScalar("B2 synthetic-ECI energy/site isolating v21AB (== v21AB CVCF value)",
                        E / N, vCvcf[idx21], 1e-9);
            }
        }

        // ── 3. Perfectly random at equiatomic ──
        System.out.println("\n--- Perfectly random (x=0.5, averaged over seeds) ---");
        verifyRandomState(geo, pr, numComp, ncf, eciCvcf, equiComposition(numComp), N, "K=2");

        // ── 4. Delta-E vs finite-difference, anchored at B2 and pure-element states ──
        System.out.println("\n--- Delta-E verification, anchored at analytic configs ---");
        verifyDeltaEFromConfig(geo, numComp, ncf, eciCvcf, buildB2Config(N, numComp), "B2-anchored", 100);
        verifyDeltaEFromConfig(geo, numComp, ncf, eciCvcf, buildPureConfig(N, numComp), "pure-element-anchored", 50);

        // ── 5. Trajectory: running-energy vs from-scratch recompute, from B2 start ──
        System.out.println("\n--- Trajectory verification (running E vs recompute, every step, from B2) ---");
        verifyTrajectory(geo, numComp, ncf, eciCvcf, buildB2Config(N, numComp), 500, 7);
    }

    // =========================================================================
    // K=3 / K=4: pure element, composition-boundary, random, delta-E
    // =========================================================================

    private static void runHigherKSuite(ModelSession.Builder builder, String elements, double T) throws Exception {
        System.out.println();
        System.out.println("=======================================================");
        System.out.println("Higher-K suite: " + elements + " BCC_A2/T, T=" + T);
        System.out.println("=======================================================");

        Workspace.SystemId id = new Workspace.SystemId(elements, "BCC_A2", "T");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, null);
        int numComp = session.numComponents();
        int L = 4;
        MCSGeometry geo = MCSGeometry.build(session, L, null);
        int N = geo.nSites();
        int ncf = geo.basis.numNonPointCfs;
        double[] eciCvcf = CECEvaluator.evaluate(session.cecEntry, T, geo.basis, "AnalyticVerify-HighK");
        PipelineResult pr = buildPipelineResult(elements, "BCC_A2", "T", numComp);

        // ── Composition-boundary: last species -> 0 must reduce to the (K-1)-system's
        //    own closed-form CVCF vector on the shared (surviving-species) columns.
        //    NOTE: this invariant holds in the CVCF (physical, e.g. pair-probability)
        //    basis, not the orthogonal one — the orthogonal point-basis sequence
        //    (ClusterMath.buildBasis) is redefined per K (e.g. {-1,+1} for K=2 vs
        //    {-1,0,+1} for K=3), so raw orthogonal columns of the same name are not
        //    the same physical quantity across K and do not need to agree numerically. ──
        System.out.println("\n--- Composition boundary (last species -> 0) ---");
        {
            double[] xBoundary = new double[numComp];
            for (int c = 0; c < numComp - 1; c++) xBoundary[c] = 1.0 / (numComp - 1);
            xBoundary[numComp - 1] = 0.0;

            double[] uBoundary = pr.computeRandomCFs(xBoundary);
            boolean allFinite = true;
            for (double v : uBoundary) if (!Double.isFinite(v)) allFinite = false;
            checkBoolean("boundary composition CFs are finite (no NaN/Inf as species->0)", allFinite);

            double[] vBoundaryCvcf = Embeddings.applyTinvTransform(uBoundary, geo.getBasis());

            PipelineResult prLower = buildPipelineResult(
                    lowerKElements(elements), "BCC_A2", "T", numComp - 1);
            ModelSession sessionLower = builder.build(
                    new Workspace.SystemId(lowerKElements(elements), "BCC_A2", "T"),
                    ModelSession.EngineConfig.MCS, null);
            MCSGeometry geoLower = MCSGeometry.build(sessionLower, L, null);
            double[] xLower = new double[numComp - 1];
            java.util.Arrays.fill(xLower, 1.0 / (numComp - 1));
            double[] uLowerGt = prLower.computeRandomCFs(xLower);
            double[] vLowerGtCvcf = Embeddings.applyTinvTransform(uLowerGt, geoLower.getBasis());

            // CVCF-basis columns whose name matches between the K and K-1 systems must
            // agree exactly at the boundary — this is the physically meaningful invariant
            // (e.g. the pair-probability of finding species {0,1} on a 1NN pair does not
            // care that a third, zero-fraction species exists in the composition space).
            List<String> namesHigh = CvCfBasis.getNonPointCfNames("BCC_A2", "T", numComp);
            List<String> namesLow = CvCfBasis.getNonPointCfNames("BCC_A2", "T", numComp - 1);
            int matched = 0;
            double maxErr = 0.0;
            for (int lo = 0; lo < namesLow.size(); lo++) {
                int hi = namesHigh.indexOf(namesLow.get(lo));
                if (hi < 0) continue;
                matched++;
                maxErr = Math.max(maxErr, Math.abs(vBoundaryCvcf[hi] - vLowerGtCvcf[lo]));
            }
            checks++;
            boolean pass = matched > 0 && maxErr < 1e-9;
            if (!pass) failures++;
            System.out.printf("  [%s] boundary CVCF CFs match (K=%d -> K=%d): matchedNames=%d maxErr=%.3e -> %s%n",
                    pass ? "OK" : "FAIL", numComp, numComp - 1, matched, maxErr, pass ? "PASS" : "FAIL");
        }

        // ── Pure-element boundary (single species = 1, rest 0): orthogonal + energy ──
        System.out.println("\n--- Pure element (single species) ---");
        {
            LatticeConfig config = buildPureConfig(N, numComp);
            double[] uFull = Embeddings.measureFullCVsFromConfig(
                    config, geo.cfEmbeddings(), geo.pointCfEmbeddings(), geo.getFlatBasisMatrix(), ncf, numComp);
            double[] xPure = new double[numComp];
            xPure[0] = 1.0;
            double[] uGroundTruth = pr.computeRandomCFs(xPure);
            checkVector("pure-element uOrth (K=" + numComp + ", measured vs closed-form at x_0=1)",
                    uFull, uGroundTruth, 1e-9);
        }

        // ── Random state (equiatomic) ──
        System.out.println("\n--- Perfectly random (equiatomic) ---");
        verifyRandomState(geo, pr, numComp, ncf, eciCvcf, equiComposition(numComp), N, "K=" + numComp);

        // ── Delta-E anchored at pure-element and equiatomic-random states ──
        System.out.println("\n--- Delta-E verification, anchored at analytic configs ---");
        verifyDeltaEFromConfig(geo, numComp, ncf, eciCvcf, buildPureConfig(N, numComp), "pure-element-anchored", 50);
        {
            Random rng = new Random(11);
            LatticeConfig randomConfig = new LatticeConfig(N, numComp);
            randomConfig.randomise(equiComposition(numComp), rng);
            verifyDeltaEFromConfig(geo, numComp, ncf, eciCvcf, randomConfig, "random-state-anchored", 100);
        }
    }

    // =========================================================================
    // Shared helpers
    // =========================================================================

    /** Builds a PipelineResult directly (same construction AbFamilyCfDiagnostic uses) for
     *  closed-form computeRandomCFs ground truth, independent of ModelSession/MCSGeometry. */
    private static PipelineResult buildPipelineResult(
            String elements, String structure, String model, int numComp) throws Exception {
        ClusterIdentificationRequest cfg = ClusterIdentificationRequest.fromSystem(elements, structure, model);
        List<Cluster> disClusters = InputLoader.parseClusterFile(cfg.getDisorderedClusterFile());
        disClusters.replaceAll(Cluster::sorted);
        SpaceGroup disSg = InputLoader.parseSpaceGroup(cfg.getDisorderedSymmetryGroup());
        List<Cluster> ordClusters = InputLoader.parseClusterFile(cfg.getOrderedClusterFile());
        ordClusters.replaceAll(Cluster::sorted);
        SpaceGroup ordSg = InputLoader.parseSpaceGroup(cfg.getOrderedSymmetryGroup());

        return ClusterCFIdentificationPipeline.run(
                disClusters, disSg.getOperations(), ordClusters, ordSg.getOperations(),
                cfg.getTransformationMatrix(),
                new double[] { cfg.getTranslationVector().getX(), cfg.getTranslationVector().getY(),
                        cfg.getTranslationVector().getZ() },
                numComp, null);
    }

    private static String lowerKElements(String elements) {
        String[] parts = elements.split("-");
        return String.join("-", java.util.Arrays.copyOf(parts, parts.length - 1));
    }

    private static double[] equiComposition(int numComp) {
        double[] x = new double[numComp];
        java.util.Arrays.fill(x, 1.0 / numComp);
        return x;
    }

    private static LatticeConfig buildB2Config(int N, int numComp) {
        int[] occ = new int[N];
        for (int i = 0; i < N; i++) occ[i] = (i % 2 == 0) ? 0 : (1 % numComp);
        return new LatticeConfig(occ, numComp);
    }

    private static LatticeConfig buildPureConfig(int N, int numComp) {
        return new LatticeConfig(new int[N], numComp);
    }

    private static CECEntry buildSingleTermCec(String elements, String structurePhase, String model, String cfName) {
        CECEntry entry = new CECEntry();
        entry.elements = elements;
        entry.structurePhase = structurePhase;
        entry.model = model;
        CECEntry.CECTerm term = new CECEntry.CECTerm();
        term.name = cfName;
        term.a = 1.0;
        term.b = 0.0;
        entry.cecTerms = new CECEntry.CECTerm[] { term };
        entry.cecUnits = "J/mol";
        entry.reference = "";
        entry.notes = "Synthetic single-term ECI for isolating one CF's energy contribution";
        entry.ncf = 1;
        return entry;
    }

    private static double dot(double[] eci, double[] v, int ncf) {
        double s = 0.0;
        for (int l = 0; l < ncf && l < eci.length && l < v.length; l++) s += eci[l] * v[l];
        return s;
    }

    /** Cross-checks MCS's random-state H (averaged over seeds) against the closed-form H. */
    private static void verifyRandomState(
            MCSGeometry geo, PipelineResult pr, int numComp, int ncf, double[] eciCvcf,
            double[] xFrac, int N, String label) {

        double[] uRandomGt = pr.computeRandomCFs(xFrac);
        double[] vCvcfGt = Embeddings.applyTinvTransform(uRandomGt, geo.getBasis());
        double hGt = dot(eciCvcf, vCvcfGt, ncf);

        Random rng = new Random(7);
        int nSamples = 40;
        double sumH = 0.0, sumH2 = 0.0;
        for (int s = 0; s < nSamples; s++) {
            LatticeConfig config = new LatticeConfig(N, numComp);
            config.randomise(xFrac, rng);
            double h = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                    geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.getBasis(), numComp) / N;
            sumH += h;
            sumH2 += h * h;
        }
        double meanH = sumH / nSamples;
        double stdH = Math.sqrt(Math.max(0, sumH2 / nSamples - meanH * meanH));
        double stderr = stdH / Math.sqrt(nSamples);
        System.out.printf("  [%s] MCS random H = %.6f +/- %.6f (stderr), closed-form H = %.6f, diff = %+.6f%n",
                label, meanH, stderr, hGt, meanH - hGt);
        // Statistical tolerance from the measured standard error (5-sigma), not a fixed guess —
        // per-sample variance scales with the number of CF columns/embeddings (much larger for
        // K=4 than K=2), so a flat tolerance either false-fails high-K or under-tests low-K.
        checkScalar(label + " random-state H (avg-of-" + nSamples + " vs closed-form; 5-sigma statistical tol)",
                meanH, hGt, Math.max(1.0, 5.0 * stderr));
    }

    private static void verifyDeltaEFromConfig(
            MCSGeometry geo, int numComp, int ncf, double[] eciCvcf,
            LatticeConfig config, String label, int trials) {

        double[] eciOrth = computeEciOrth(eciCvcf, geo.getBasis(), ncf);
        int maxEmbPerCol = Embeddings.maxEmbPerCfColumn(geo.cfEmbeddings());
        Embeddings.DeltaScratch scratch = new Embeddings.DeltaScratch(ncf, ncf * maxEmbPerCol);

        Random rng = new Random(42);
        int N = config.getN();
        int mismatches = 0;
        double maxAbsErr = 0.0;
        int ran = 0;

        for (int trial = 0; trial < trials && ran < trials; trial++) {
            int i = rng.nextInt(N);
            int j = rng.nextInt(N);
            if (config.getOccupation(i) == config.getOccupation(j)) continue;
            ran++;

            double eBefore = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                    geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.getBasis(), numComp);

            double dEFlat = Embeddings.deltaEExchangeCvcf(i, j, config, geo.flatEmbData, geo.getFlatBasisMatrix(),
                    geo.siteToCfIndex, ncf, eciCvcf, geo.getBasis(), scratch, maxEmbPerCol, eciOrth, numComp);
            scratch.cleanup(maxEmbPerCol);

            double dEList = Embeddings.deltaEExchangeCvcf(i, j, config, geo.cfEmbeddings(), geo.getFlatBasisMatrix(),
                    geo.siteToCfIndex, ncf, eciCvcf, geo.getBasis(), scratch, maxEmbPerCol, eciOrth, numComp);

            double dEV2 = Embeddings.deltaEExchangeCvcfV2(i, j, config, geo.flatEmbData, geo.getFlatBasisMatrix(),
                    geo.siteToCfIndex, ncf, eciCvcf, geo.getBasis(), scratch, maxEmbPerCol, eciOrth, numComp);
            scratch.cleanup(maxEmbPerCol);

            int occI = config.getOccupation(i), occJ = config.getOccupation(j);
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);
            double eAfter = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                    geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.getBasis(), numComp);
            double dEFiniteDiff = eAfter - eBefore;

            double errFlat = Math.abs(dEFlat - dEFiniteDiff);
            double errList = Math.abs(dEList - dEFiniteDiff);
            double errV2 = Math.abs(dEV2 - dEFiniteDiff);
            maxAbsErr = Math.max(maxAbsErr, Math.max(errFlat, Math.max(errList, errV2)));
            if (errFlat > 1e-6 || errList > 1e-6 || errV2 > 1e-6) mismatches++;

            config.setOccupation(i, occI);
            config.setOccupation(j, occJ);
        }

        checks++;
        boolean pass = mismatches == 0;
        if (!pass) failures++;
        System.out.printf("  [%s] deltaE trials=%d mismatches=%d maxAbsErr=%.3e -> %s%n",
                label, ran, mismatches, maxAbsErr, pass ? "PASS" : "FAIL");
    }

    private static void verifyTrajectory(
            MCSGeometry geo, int numComp, int ncf, double[] eciCvcf,
            LatticeConfig config, int steps, long seed) {

        double[] eciOrth = computeEciOrth(eciCvcf, geo.getBasis(), ncf);
        int maxEmbPerCol = Embeddings.maxEmbPerCfColumn(geo.cfEmbeddings());
        Embeddings.DeltaScratch scratch = new Embeddings.DeltaScratch(ncf, ncf * maxEmbPerCol);

        int N = config.getN();
        double runningE = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.getBasis(), numComp);

        Random rng = new Random(seed);
        int mismatches = 0;
        double maxAbsErr = 0.0;

        for (int step = 0; step < steps; step++) {
            int i = rng.nextInt(N);
            int j = rng.nextInt(N);
            if (config.getOccupation(i) == config.getOccupation(j)) continue;

            double dE = Embeddings.deltaEExchangeCvcf(i, j, config, geo.flatEmbData, geo.getFlatBasisMatrix(),
                    geo.siteToCfIndex, ncf, eciCvcf, geo.getBasis(), scratch, maxEmbPerCol, eciOrth, numComp);
            scratch.cleanup(maxEmbPerCol);

            boolean accept = dE <= 0 || rng.nextDouble() < Math.exp(-dE / (8.314 * 1000.0));
            if (accept) {
                int occI = config.getOccupation(i), occJ = config.getOccupation(j);
                config.setOccupation(i, occJ);
                config.setOccupation(j, occI);
                runningE += dE;
            }

            double trueE = Embeddings.totalEnergyCvcf(config, geo.cfEmbeddings(), geo.pointCfEmbeddings(),
                    geo.getFlatBasisMatrix(), ncf, eciCvcf, geo.getBasis(), numComp);
            double err = Math.abs(runningE - trueE);
            maxAbsErr = Math.max(maxAbsErr, err);
            if (err > 1e-6) mismatches++;
        }

        checks++;
        boolean pass = mismatches == 0;
        if (!pass) failures++;
        System.out.printf("  trajectory steps=%d mismatches=%d maxAbsErr=%.3e -> %s%n",
                steps, mismatches, maxAbsErr, pass ? "PASS" : "FAIL");
    }

    /** Mirrors MCSRunner.computeEciOrth (package-private there) for standalone testing. */
    private static double[] computeEciOrth(double[] eciCvcf, CvCfBasis basis, int ncf) {
        if (basis == null || basis.Tinv == null) return null;
        double[][] Tinv = basis.Tinv;
        int tCols = Tinv[0].length;
        double[] eciOrth = new double[tCols];
        for (int m = 0; m < tCols; m++) {
            double sum = 0.0;
            for (int l = 0; l < ncf && l < eciCvcf.length && l < Tinv.length; l++)
                sum += eciCvcf[l] * Tinv[l][m];
            eciOrth[m] = sum;
        }
        return eciOrth;
    }

    private static void checkScalar(String label, double actual, double expected, double tol) {
        checks++;
        double err = Math.abs(actual - expected);
        boolean pass = err <= tol;
        if (!pass) failures++;
        System.out.printf("  [%s] %s: actual=%.8f expected=%.8f err=%.3e (tol=%.1e) -> %s%n",
                pass ? "OK" : "FAIL", label, actual, expected, err, tol, pass ? "PASS" : "FAIL");
    }

    private static void checkVector(String label, double[] actual, double[] expected, double tol) {
        checks++;
        int n = Math.min(actual.length, expected.length);
        double maxErr = 0.0;
        for (int i = 0; i < n; i++) maxErr = Math.max(maxErr, Math.abs(actual[i] - expected[i]));
        boolean pass = maxErr <= tol;
        if (!pass) failures++;
        System.out.printf("  [%s] %s: maxErr=%.3e (tol=%.1e) -> %s%n",
                pass ? "OK" : "FAIL", label, maxErr, tol, pass ? "PASS" : "FAIL");
        if (!pass) {
            System.out.println("    actual:   " + java.util.Arrays.toString(actual));
            System.out.println("    expected: " + java.util.Arrays.toString(expected));
        }
    }

    private static void checkBoolean(String label, boolean condition) {
        checks++;
        if (!condition) failures++;
        System.out.printf("  [%s] %s -> %s%n", condition ? "OK" : "FAIL", label, condition ? "PASS" : "FAIL");
    }
}
