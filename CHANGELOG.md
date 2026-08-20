# Changelog

All notable changes to this project are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/) — dated entries summarize milestones rather than individual commits, since this project doesn't cut versioned releases.

See `CLAUDE.md` for the architecture, layer contracts, and session/API contracts; `README.md` for project overview and API usage.

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

Aligning the MCS calculation pipeline with the strict three-layer contract
(`ui` → `calculation` → `model`).

**Status: Phase 1 ACHIEVED (by a different route than planned); Phase 2 still open**

**Phase 1 — model building belongs to the model layer.** The original plan was to add
`MCSRunner.Builder.session(ModelSession, double)`. That builder no longer exists; the
same goal was reached with a static factory instead:

```java
MCSGeometry geo = MCSGeometry.build(session, L, sink);          // geometry, cached per (session, L)
MCSRunner runner = MCSRunner.forTemperature(geo, session, T, sink);  // ECI eval + transform
```

`ThermodynamicWorkflow.runMcs` now passes the `ModelSession` straight through and no
longer decomposes it into raw fields or calls `CECEvaluator` itself — which is what
Phase 1 was for. No further work needed here; the detailed builder specification that
used to occupy this section has been removed as obsolete.

*Known gap:* `MCSStatisticsProcessor` (calculation layer, τ_int / block averaging /
jackknife Cv) exists but is **not** wired into `runMcs`. MCS results currently carry
`Double.NaN` for `stdEnthalpy` and `heatCapacity`. Wiring it up is unfinished Phase-1
adjacent work.


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

## 2026-08-19 — 2026-08-20 — CVM evaluator/solver split, SGTE reference energy, Hillert solver

Undertaken to make the `G` expressions auditable: `CVMGibbsModel` was both evaluator and optimizer, which is what made them hard to review in isolation.

- **`CVMGibbsModel` reduced to a pure evaluator.** It answers for G/H/S, their derivatives, and SRO at given system parameters (elements, structure, ECIs), macro parameters (T, x), and micro parameters (u) — and owns no iteration or convergence logic. Results are read through a nested `State` obtained from `model.at(T, x, u)`; scalars compute on demand. Solvers now hold a model and drive it from outside.
- **Split `model/cvm/` (evaluates) from `model/equilibrium/` (solves).** The Newton–Raphson loop moved out of `CVMGibbsModel` into `CvmNewtonSolver`, its nine stages documented inline; two subtleties are noted where they were previously implicit — the gradient test uses an L1 norm (matching the reference, not Euclidean), and the final step test reads the raw Newton step, not the clamped one.
- **Extracted `CvmGeometry`** — the immutable Stage 1–4 pipeline product, cluster algebra only, independent of any Hamiltonian, so geometry can be checked on its own inputs. Carries a class-level TODO flagging ordered phases (point set wider than K, sublattice orbit splitting) as unreviewed.
- **Replaced the hardcoded lattice-stability tables with `SgteDatabase`**, parsing SGTE Unary v4.4 from `inputs/unary.dat`, and reduced `LatticeStability` to a façade over it (867 → 117 lines). This supplies the pure-element reference `G0m`, without which chemical potentials are not comparable across phases. Verified by 176 cross-check agreements against the previous hardcoded values, analytic T-derivatives against finite differences, and external reference values.
- **Added `HillertSolver`** — multi-phase equilibrium in one file, with `Phase`, `Result`, `PhaseResult`, `PhaseStep`, and `EquilibriumMatrix` nested inside it. Each had exactly one caller and is meaningless outside the outer loop.
- **Merged `SroCalculator` into `CVMGibbsModel.State`** — SRO is a thermodynamic property of a state, like G/H/S. Also removed `CvmEvaluator`/`CvmState` after they had served as the transitional pair.
- Net effect on the two packages: fourteen classes down to seven.
- Verified throughout by parity gates against frozen copies of the pre-refactor implementations, so each step compared against real prior behaviour rather than against itself. The three `CLAUDE.md` reference values are bit-identical.

**Known open, deferred:** the Hillert path is validated only for a single phase; a two-phase reference trace is still needed. It also does not converge at one near-edge composition (Mo-Nb-Ta, 1273 K, x=[0.05, 0.475, 0.475]) where the reference does in 7 outer iterations — the outer loop is a faithful port, so the discrepancy is upstream in the Hessian. Separately, `dGm/du` disagrees with the Mathematica reference on all 18 components for Mo-Nb-Ta while agreeing with finite differences of our own `Gm`; this is pre-existing and unresolved.

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
