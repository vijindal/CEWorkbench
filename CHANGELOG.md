# Changelog

All notable changes to this project are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/) — dated entries summarize milestones rather than individual commits, since this project doesn't cut versioned releases.

See `CLAUDE.md` for the pipeline architecture and session contract, `ARCHITECTURE.md` for the layer contracts, `README.md` for project overview.

---

## [Unreleased]

Active and planned work across subsystems.

### 1. CVCF Basis Derivation

Tracks progress on `CvCfBasis.java` registrations (Stage 4 of cluster/CF identification) across structures and component counts.

**Status: ACTIVE — BCC_B2 binary complete, ternary in progress**

**Verified and working** (regression-tested, 0.00e+00 self-test diff)

| Structure | Model | K | Status |
|---|---|---|---|
| BCC_A2 | T | 2, 3, 4 | Verified |
| FCC_A1 | T | 2, 3, 4 | Verified |
| HCP_A3 | T | 2 (binary) | Verified — rebuilt from Jindal & Lele (2025) CALPHAD paper, Appendix 2, eq 59-64 |
| HCP_A3 | T | 3 (ternary) | Verified — rebuilt from paper eq 65+ with corrected site numbering (math site→our site: 1→1, 2→4, 3→2, 4→3, 5→5, 6→6); tetrahedron ABC-mixed terms are literal Mathematica translations, not paper-inferred |
| BCC_B2 | T | 2 (binary) | Verified — derived by hand (see below), rank 9/9, exact 0.00e+00 diff, η=0 at equiatomic |

**Blocked / not started**

| Structure | Model | K | Status |
|---|---|---|---|
| HCP_A3 | T | 4 (quaternary) | Blocked — registration exists but M-matrix is singular (pipeline reports tcf=85, needs 86-wide basis, but our 84-CV registration is short by 1). Table 19 says 83+1=84 should be right; real pipeline disagrees. Unresolved discrepancy, same status as before the HCP rework. |
| BCC_B2 | T | 3 (ternary) | In progress — see below |
| BCC_B2 | T | 4 (quaternary) | Not started |

**BCC_B2 binary — derivation summary**

B2 has two sublattices: α = {p1, p4}, β = {p2, p3} (same 4-site tetrahedron as A2, `BCC_A2-T.txt`'s coordinates, just partitioned per `BCC_B2-T.txt`'s `{{p1,p4}},{{p2,p3}}` grouping).

Orthogonal-CF orbit structure (tcf=8, from pipeline Stage 2b):
- tetrahedron (t=0): unsplit, ncv=9
- triangle (t=1): splits into {1,2,3} and {1,2,4}
- II-n pair (t=2, mh=4 in A2): unsplit as **one** pipeline orbit, but α-α ({1,4}) and β-β ({2,3}) are physically distinct once sublattice occupation matters — picked as **two** independent CVCF variables anyway (see below)
- I-n pair (t=3, mh=3 in A2): two orbit *entries* in Stage 2 output, but same physical pair type — **one** CVCF variable, unchanged from A2
- point (t=4): splits into α and β

Final 9-variable CVCF basis:

| name | formula | ECI shared with A2's... |
|---|---|---|
| V4AB | `p1[A]p2[A]p3[B]p4[B]` | new — self-paired tetrahedron orbit (AABB≡ABAB≡BABA≡BBAA under (1↔4)(2↔3); **no antisymmetric partner exists**, must use plain product) |
| V31AB | `p1[A]p2[B]p3[B] − p1[B]p2[A]p3[A]` | v3AB |
| V32AB | `p2[A]p1[B]p4[B] − p2[B]p1[A]p4[A]` | v3AB (same ECI as V31AB) |
| V221AB | `p1[A]p4[B]` | v22AB |
| V222AB | `p2[A]p3[B]` | v22AB (same ECI as V221AB) |
| V21AB | `p1[A]p2[B]` | v21AB (unchanged) |
| xA, xB | `p1[A]`, `p1[B]` | — |
| eta | `p1[A] − p2[A]` (LRO parameter, no ECI) | — |

**Pitfall found the hard way:** the antisymmetric combination `p1[A]p4[B] − p1[B]p4[A]` for the II-n pair evaluates to an *exact zero row* in the real M-matrix — sites 1 and 4 (and separately 2, 3) have no AB/BA distinction at that orbit under B2's true space group. Only the **plain product** works. Verify any new antisymmetric pick against the real M-matrix (rank check), not just plausibility — see the verification pattern below.

**Shared-code fixes** (apply to all structures, not just B2)

1. `CvCfBasis.Definition` gained an explicit `numPointCfs` field (default = `numComponents`, matching prior behavior for A2/FCC/HCP). B2 registers `numPointCfs=3` (xA, xB, eta) since it has one more point-like quantity than `numComponents`. Use the `register(..., numPointCfs)` overload for any future ordered structure with order parameters.
2. `ClusterCFIdentificationPipeline.computeRandomCFs` previously hardcoded `pointCfCount = K-1`, silently leaving split point-orbit columns at 0 instead of the correct random-state value (this produced a spurious η=0.5 at equiatomic instead of the correct η=0). Fixed to use the tracked `nxcf` field. This affects the "random CF at composition" self-test for any ordered structure with more than `K-1` point orbits.
3. `ClusterCFIdentificationPipeline.runFullWorkflow` now checks `CvCfBasis.isSupported(structurePhase, ...)` directly before falling back to `resolveParentStructure(structurePhase)` — previously an ordered structure with its own registered CVCF definition (like the new `BCC_B2_T_2`) was silently resolved to its disordered parent's definition instead.

**Next session: BCC_B2 ternary (K=3)**

Confirmed via pipeline: tcf=35 (needs 36-variable basis). Per-cluster-type CV counts (from `classified[t][j][k]` Stage-3 breakdown, counting `(j,k)` entries):

| Cluster type | # CVCF variables |
|---|---|
| Tetrahedron | 9 |
| Triangle | 12 |
| II-n pair (t=2, A2 mh=4) | 4 |
| I-n pair (t=3, A2 mh=3) | 6 |
| Point (incl. order parameters) | 4 |
| **Total** | 35 |

Point CFs confirmed via literature research (Sanchez–Ducastelle–Gratias generalized cluster description): **xA, xB, xC, ηB, ηC** where ηB = xB(α)−xB(β), ηC = xC(α)−xC(β), and ηA = −(ηB+ηC) is dependent (traceless constraint). This matches the pipeline's nxcf=4 (independent point DOF) + xA as a conventional redundant 5th.

**Open question, not yet resolved:** does the ternary I-n pair (t=3) split into 6 CVCF variables (2 per species-pair AB/AC/BC), or does it stay unsplit like binary's I-n pair (which needed only 1 variable, `V21AB`, despite Stage 2 showing two `j`-orbit entries)? The `classified[3][j][k]` data shows k=0,1 for all three `j` (species-pairs), suggesting a split — but binary's analogous `classified[3][0][k=0,1]` *also* showed a k-split and yet the correct answer (verified by rank + self-test) was **no split needed**. This "classified k-split ≠ CVCF split" trap means per-cluster-type CV counts from `classified[]` should be treated as an upper bound / sanity check, not a direct prescription — always verify via M-matrix rank, same as the II-n antisymmetric-cancels-to-zero pitfall above.

Next step: propose I-n pair as 3 unsplit variables (AB, AC, BC, matching A2 exactly) first — cheaper to test — and only fall back to 6 if rank comes up short.

**Verification pattern** (use for any new registration)

1. Get real orbit structure: `./gradlew run --args="type1a <elements> <structure> T --verbose"`, grep `t=X j=Y:` lines for `ncv`/`wcv`, and `[RESULT] lcv:` for orbit counts.
2. Propose CVCF picks, wire into `CvCfBasis.java` via `register(...)`.
3. Compile, run the same command — watch for "Matrix is singular" errors.
4. If singular: temporarily add M-matrix row dump (`for i: emit("MROW "+i+" ("+name+"): "+row)` gated behind an env var, e.g. `DEBUG_MMATRIX`), pull into Python/numpy, check `np.linalg.matrix_rank` and SVD null vector to find which CV/orbit direction is missing or redundant. Remove the debug block once resolved — don't leave it in.
5. Once rank is full: check `[SELF-TEST] CV VERIFICATION` output — every cluster type should show `Diff: 0.00e+00` (or ~1e-16 floating noise) at equiatomic composition. Any order parameter (like η) should read exactly 0 at equiatomic/random state.
6. Regression-check all previously-verified structures before calling it done — the shared-code fixes above can silently affect unrelated structures.

### 2. AlloyMC Monte Carlo Engine

Implementation of a high-performance Metropolis Monte Carlo engine for Cluster Expansion models.

**Status: ACTIVE**

> **Ternary Mapping Resolved** — the previously identified -2.0 energy discrepancy in the Nb-Ti-V system has been resolved. The root cause was a basis function indexing mismatch between the engine's state vector population and the CVCF transformation matrix.

Key resolutions:
- **Dynamic Basis Mapping**: refactored `AlloyMC` to use `MCSGeometry.getCfBasisIndices()` metadata for point CF column assignments.
- **Full Vector Transformation**: updated `MCSGeometry` to transform the entire measured orthogonal vector ($T^{-1}u$) rather than hard-appending compositions.
- **Symbolic Verification**: integrated symbolic definitions ($p[A]*p[B]$) directly into the `CvCfBasis` and `AlloyMCTest` diagnostic suite.

Roadmap:
1. **Finalize Hamiltonian Integration** — integrate the verified CVCF state vector into the `AlloyMC` energy calculator; verify mixing energy for the Nb-Ti-V equiatomic disordered state.
2. **Metropolis Sampling Loop** — implement the standard Metropolis-Hastings acceptance criterion; optimize site-swap logic for local energy updates.
3. **Verification & Benchmarking** — validate against reference CVM results for simple BCC alloys; benchmark site-swap performance for large supercells ($L > 20$).

Verification plan:
- *Automated*: `AlloyMCTest` must produce 0.0 J for pure Nb in the ternary Nb-Ti-V system; verify randomized ternary energy is plausible.
- *Manual*: trace the $T^{-1}$ matrix construction for a ternary tetrahedron to ensure linear independence and correct endpoint centering.

### 3. MCS Pipeline Refactor

Aligning the MCS (Monte Carlo Sampling) calculation pipeline with the strict three-layer contract (`ui` → `calculation` → `model`). Two phases, in dependency order — Phase 1 is a prerequisite for Phase 2's ownership/queueing work.

**Status: PLANNED — neither phase implemented yet**

**Phase 1: Move model building into the model layer** (`MCSRunner` accepts `ModelSession`)

*Context:* the calculation layer should pass model specifications to the model layer, which takes full responsibility for building the model (cluster algebra, CEC reading, initial config, ECI evaluation and transform). Steps 2 and 3 (Hamiltonian evaluation, ECI transform) should happen together with Step 4 (initialization) inside the model layer — not scattered across `ThermodynamicWorkflow` and `MCSRunner`.

*Current problem:*
- `ThermodynamicWorkflow.runMcs()` does Step 2 (`CECEvaluator.evaluate` → `eci[]`) and validates the C-matrix before handing off to `MCSRunner` — model-layer concerns leaked into the calculation layer.
- `MCSRunner.Builder` requires the caller to decompose `ModelSession` into 6 raw fields (`clusterData`, `eci`, `basis`, `matrixData`, `lcf`, `numComp`) — the session was already available.
- `MCSRunner.run()` does Step 3 (ECI transform via `Tinv`) — fine, stays there.

*Target:*
- `ThermodynamicWorkflow.runMcs()` becomes a thin dispatcher — passes `ModelSession` directly.
- `MCSRunner.Builder` gains `session(ModelSession, double temperature)` — extracts all model fields internally, performs `CECEvaluator.evaluate()`, validates the C-matrix, builds ECIs.
- Steps 2, 3, 4 all visible and sequenced inside the model layer.

Layer responsibilities after the change:
```
UI Layer
  → ModelSpecifications (elements, structure, model, engine)
  → CalculationSpecifications (T, x, L, nEquil, nAvg)
  → CalculationService.execute()

Calculation Layer (ThermodynamicWorkflow.runMcs)
  → retrieves pre-built ModelSession                (unchanged)
  → wires progress callbacks                        (unchanged)
  → calls MCSRunner.builder().session(session, T)   ← NEW: passes session directly
  → calls builder.build().run()                     (unchanged)
  → runs MCSStatisticsProcessor on MCResult          (unchanged)

Model Layer (MCSRunner.Builder)
  Step 2: CECEvaluator.evaluate(session.cecEntry, T, session.cvcfBasis) → eci (CVCF basis)
  Step 2b: validate C-matrix dimensions            ← moved from ThermodynamicWorkflow
  Step 3: ECI transform via Tinv (stays in MCSRunner.run() as today)
  Step 4: build lattice, embeddings, LatticeConfig (stays in MCSRunner.run() as today)
```

Precise changes:

1. `MCSRunner.Builder` — add `session()` method, remove `eci()` requirement.

   New import in `MCSRunner.java`:
   ```java
   import org.ce.model.ModelSession;
   import org.ce.model.hamiltonian.CECEvaluator;
   import org.ce.model.PhysicsConstants;
   ```

   New field in `Builder`:
   ```java
   private ModelSession session = null;   // set by session() convenience method
   ```

   New `Builder` method:
   ```java
   public Builder session(ModelSession s, double temperature) {
       this.session = s;
       this.T       = temperature;
       return this;
   }
   ```

   Modified `build()` in `Builder` — add block that auto-fills from session if provided:
   ```java
   public MCSRunner build() {
       if (session != null) {
           // Step 2: Evaluate Hamiltonian — ECI in CVCF basis at temperature T
           var matData = session.clusterData.getMatrixData();
           int cmatCols = /* extract from matData */;
           if (cmatCols != session.cvcfBasis.totalCfs())
               throw new IllegalStateException("C-matrix dimension mismatch ...");
           this.eci       = CECEvaluator.evaluate(session.cecEntry, T, session.cvcfBasis, "MCS");
           // Fill remaining fields from session
           this.clusterData = session.clusterData.getDisorderedClusterResult().getDisClusterData();
           this.numComp     = session.numComponents();
           this.basis       = session.cvcfBasis;
           this.matrixData  = matData;
           this.lcf         = session.clusterData.getLcf();
           if (R <= 0) this.R = PhysicsConstants.R_GAS;
       }
       // existing validation follows unchanged
       if (clusterData == null) throw new IllegalStateException("clusterData required");
       if (eci == null)         throw new IllegalStateException("eci required");
       if (T <= 0)              throw new IllegalStateException("T must be > 0");
       ...
   }
   ```

2. `ThermodynamicWorkflow.runMcs()` — remove Steps 2 & 3, use `session()` builder method.

   Remove from `runMcs()`:
   - `CECEvaluator.evaluate(...)` call (Step 2 — moves to `MCSRunner.Builder.build()`)
   - C-matrix dimension validation block (Step 2b — moves to `MCSRunner.Builder.build()`)
   - All 13 individual `.clusterData()`, `.eci()`, `.numComp()`, `.basis()`, `.matrixData()`, `.lcf()` builder calls

   Replace with:
   ```java
   MCSRunner.Builder builder = MCSRunner.builder()
       .session(session, request.temperature)   // model building now in model layer
       .composition(request.composition)
       .nEquil(nEquil)
       .nAvg(nAvg)
       .L(L)
       .seed(System.currentTimeMillis())
       .cancellationCheck(Thread.currentThread()::isInterrupted);
   ```

   Remove import from `ThermodynamicWorkflow.java`:
   ```java
   // remove:  import org.ce.model.hamiltonian.CECEvaluator;
   // remove:  import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
   ```

3. No other file changes — `ModelSession`, `MetropolisMC`, `CalculationService`, `ThermodynamicRequest`, `LocalEnergyCalcTest` all unchanged.

Critical files:
- [MCSRunner.java](src/main/java/org/ce/model/mcs/MCSRunner.java) — Builder gains `session()` + `build()` fills fields from session
- [ThermodynamicWorkflow.java](src/main/java/org/ce/calculation/workflow/thermo/ThermodynamicWorkflow.java) — `runMcs()` simplified to ~10 builder lines
- [ModelSession.java](src/main/java/org/ce/model/ModelSession.java) — read-only, no changes
- [MetropolisMC.java](src/main/java/org/ce/model/mcs/MetropolisMC.java) — no changes
- [LocalEnergyCalcTest.java](src/test/java/org/ce/model/mcs/LocalEnergyCalcTest.java) — no changes

Verification:
```bash
./gradlew compileJava compileTestJava   # must BUILD SUCCESSFUL
./gradlew test                          # LocalEnergyCalcTest must pass (3 tests)
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 0.5 --verbose"
```

**Phase 2: Ownership, queueing, and result streaming**

*Context & objectives:* simplify the MCS calculation pipeline further by:
1. **Strictly classifying parameters** — moving physical parameters (like supercell size `L`) that require a model rebuild into `ModelSpecifications`.
2. **Clear ownership** — `CalculationService` becomes the authoritative owner of the persistent model instance.
3. **Queueing** — implementing a synchronized execution queue in the calculation layer.
4. **Streaming results** — providing point-by-point thermodynamic feedback to the UI while deferred plotting occurs at the end of a scan.

Specification redesign (`org.ce.calculation`):

[MODIFY] `CalculationDescriptor.java`
- Relocate `Parameter.MCS_L` from `JobSpecifications` to `ModelSpecifications`.
- Update `Registry.getRequirements` logic to ensure `L` is treated as a construction-time requirement for MCS models.

Model layer autonomy (`org.ce.model`):

[NEW] `ThermodynamicModel.java` — a unified domain interface for physics engines.
- `void validate(ModelSpecifications specs)`: internal validation for physical consistency (e.g. $L \ge 1$).
- `void calculate(T, x, JobSpecifications specs, Consumer<ThermodynamicResult> resultSink)`: performs the actual physics and streams results.

[NEW] `McsModel.java` & `CvmModel.java`
- Concrete implementations of `ThermodynamicModel`.
- `McsModel` encapsulates `MCSGeometry` and ensures it persists for the life of the engine instance.

Authorized ownership & queueing (`org.ce.calculation.workflow`):

[MODIFY] `CalculationService.java`
- **Active model management**: hold `private ThermodynamicModel activeModel`.
- **Job queueing**: implement a serialized runner (e.g. `singleThreadExecutor`) to prevent overlapping calculations on a single persistent engine.

Result streaming & plotting (`org.ce.calculation.workflow.thermo`):

[MODIFY] `ThermodynamicWorkflow.java`
- Update implementation to stream results point-by-point via the `resultSink`.
- Append a completion signal to notify the UI to finalize plots and statistical summaries.

Verification plan:
- *Automated*: queueing verification (concurrent request testing to ensure serial execution); contract verification (unit tests for `ModelSpecification` validation within the engines).
- *Manual*: GUI feedback (monitoring the "Building..." vs "Running..." states in the terminal); scan visibility (verifying that chart points appear incrementally during a temperature scan).

---

## 2026-08-06 — HCP/B2 CVCF basis rebuild, documentation consolidation

- Fixed `ClusterCFIdentificationPipeline.runFullWorkflow` passing the wrong `maxClusters` (disordered parent instead of ordered phase) to `CMatrixPipeline.run` — the root cause of `BCC_B2` Stage 3 identification failures.
- Rebuilt `HCP_A3` binary and ternary CVCF registrations from Jindal & Lele (2025, CALPHAD 89) Appendix 2, with corrected site numbering; both verified to exact 0.00e+00 self-test diff.
- Derived and registered `BCC_B2` binary (K=2) CVCF basis by hand — see [Unreleased] for details.
- Added `CvCfBasis.Definition.numPointCfs` and fixed `computeRandomCFs` to support ordered structures with order parameters.
- Moved remaining `clus`/`sym` data files from `src/main/resources` into `data/CEWorkbench/inputs`, the canonical `Workspace` location.
- Consolidated 8 root-level markdown files into 4 (`README.md`, `CLAUDE.md`, `ARCHITECTURE.md`, `STATUS.md` → later folded into this changelog); removed stale debug scratch files.

## 2026-05-01 — AlloyMC performance optimization

- Optimized AlloyMC energy sampling from ~14.3 s/composition to ~2.3 s/composition across eleven incremental optimization passes (per-move allocation elimination, precomputed orthogonal ECI, CSR-layout basis indexing, parallel spatial-decomposition MCS, incremental CF updates, thread-local scratch buffers, flat embedding arrays).
- Implemented the `AlloyMC` stateful engine with CVCF basis transformation and ternary verification; resolved the ternary CF mapping/energy discrepancy.

## 2026-04-20 — 2026-04-29 — MCS engine consolidation and geometry refactor

- Consolidated the MCS package from 10 files to 5; renamed `MCEngine` → `MetropolisMC`.
- Refactored `MCSRunner` as a persistent model with per-step progress output.
- Reorganized `Embeddings.java` and the broader MCS folder structure by algorithm step.
- Refactored MCS geometry and embeddings; added verification scripts.

## 2026-04-10 — 2026-04-17 — Three-layer architecture established

- Organized the codebase into Model, Calculation, and UI layers with strict `ui → calculation → model` dependency direction; eliminated upward-dependency violations (`ModelSession`, `ProgressEvent`, `ThermodynamicResult`/`EquilibriumState`).
- Moved CVM and MCS optimizers into the model layer; separated MCS statistical post-processing (`MCSStatisticsProcessor`) from simulation.
- Refactored the calculation layer to be unified and specification-driven — CLI and GUI consolidated under a common `CalculationService` dispatcher using `CalculationSpecifications`/`CalculationDescriptor`.
- Implemented metadata-driven dynamic GUI generation and multicomponent (N-component) thermodynamic scanning.

## 2026-03-25 — 2026-04-08 — CVCF basis pipeline introduced

- Completed CVCF basis transformation matrices and registry for `BCC_A2`; added the CVCF basis integration adapter and test suite.
- Implemented CEC basis transformation (binary → ternary) and cross-order CEC usage projection.
- Implemented dynamic CVCF matrix generation, replacing static per-structure transformation classes.
- Stabilized the CVM ternary solver (fixed Jackson deserialization for `nijTable`, aligned enthalpy formula with the Mathematica reference).

## 2026-03-16 — Project foundation

- Established the initial layered architecture for the CVM identification pipeline.
- Added the thermodynamic calculation framework and CEC (Hamiltonian) management.
- Initial README with architecture overview and quick start guide.
