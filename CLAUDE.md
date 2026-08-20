# CLAUDE.md — CE Thermodynamics Workbench

Agent context file. Read this before touching any code.

---

## What this project is

A Java/Swing desktop application for **Cluster Expansion (CE) thermodynamic calculations** on alloy systems. It identifies cluster basis functions, manages effective cluster interaction (ECI) Hamiltonians, and computes free energy (G/H/S) using either the Cluster Variation Method (CVM) or Monte Carlo (MCS).

Three types of work:
- **Type-1a** — Cluster identification (4 stages: geometry → CF basis → C-matrix → CVCF transform)
- **Type-1b** — Scaffold an empty Hamiltonian JSON from cluster identification results
- **Type-2** — Thermodynamic equilibrium calculation (CVM Newton–Raphson or MCS)

---

## Build and run

```bash
# GUI
./gradlew runGui

# CLI — full pipeline with defaults
./gradlew run --args="all Nb-Ti BCC_A2 T"

# CLI — single-point CVM calculation
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5 --verbose"

# CLI — temperature scan
./gradlew run --args="type2 Nb-Ti BCC_A2 T --verbose"

# Compile only
./gradlew compileJava
```

All Gradle tasks: `run` (CLI), `runGui`, `runGuiDebug`, `runCli`, `runDebugCli`,
`runScratch` (`-PscratchClass=...`), `build`, `installDist` (builds the launcher the
JSON API is invoked through).

---

## Package layout

```
org.ce
├─ CEWorkbench.java          GUI main()
├─ CEWorkbenchContext.java   App-level wiring; shared by GUI and CLI
│
├─ model/                    Physics evaluators, optimizers, persistent state, disk I/O
│   ├─ ModelSession.java     Immutable session (see below — most important class)
│   ├─ ThermodynamicResult   Immutable result DTO (G/H/S, CFs, SRO, convergence)
│   ├─ ProgressEvent         Structured progress events for charts
│   ├─ PhysicsConstants      R_GAS etc.
│   ├─ cluster/              Cluster geometry, CF basis, C-matrix, CVCF pipeline
│   │   ├─ ClusterCFIdentificationPipeline  Stages 1–4, produces PipelineResult
│   │   ├─ CMatrixPipeline   C-matrix build + evaluateCVs (cluster probabilities)
│   │   └─ StructurePhaseRegistry, SpaceGroup, ClusterMath, LinearAlgebra, …
│   ├─ cvm/                  Pure evaluators — no solver loops (see note below)
│   │   ├─ CVMGibbsModel     Evaluator; nested State/Shell/PairSro
│   │   ├─ CvmGeometry       Immutable Stage 1–4 product; cluster algebra only
│   │   └─ CvCfBasis         CVCF basis registry + transformation matrices
│   ├─ equilibrium/          Solvers + pure-element reference energy
│   │   ├─ CvmNewtonSolver   Fixed-composition Newton–Raphson (nine stages)
│   │   ├─ HillertSolver     Multi-phase; nested Phase/Result/PhaseStep/EquilibriumMatrix
│   │   ├─ LatticeStability  G0m façade
│   │   └─ SgteDatabase      SGTE Unary v4.4 parser (inputs/unary.dat)
│   ├─ mcs/                  MCS supercell state, geometry, Metropolis engine
│   │   ├─ MCSGeometry       Expensive per-(session, L) geometry; built once
│   │   ├─ MCSRunner         forTemperature(): ECI evaluation + transform
│   │   ├─ MetropolisMC      Sweep loop, accept/reject, raw observables
│   │   ├─ LatticeConfig     Atomic occupation array; mutable supercell state
│   │   ├─ Embeddings, LatticeDecomposer, McsSuggester
│   │   └─ AlloyMC           Standalone stateful MC engine (parallel chains)
│   ├─ hamiltonian/          CECEntry (+CECTerm), CECEvaluator, NumericalCECTransformer
│   └─ storage/              Workspace (paths, SystemId), InputLoader, DataStore
│
├─ calculation/              Public API, vocabulary, discovery, dispatch
│   ├─ Conditions/Range/ConditionsScan  Named conditions (see API contract below)
│   ├─ CalculationDescriptor Vocabulary: Property, Mode, Parameter, Registry
│   ├─ EciValidator          Strict ECI name/coverage validation
│   ├─ CalculationResult, QuantityDescriptor, ResultFormatter
│   └─ workflow/
│       ├─ CalculationService        calculate/calculateScan + execute
│       ├─ CECManagementWorkflow     Hamiltonian scaffold/load/save
│       └─ thermo/ ThermodynamicWorkflow, MCSStatisticsProcessor
│
└─ ui/
    ├─ cli/   Main.java, ApiCommand.java (JSON API)
    └─ gui/   MainWindow, DynamicCalculationPanel, WorkbenchContext, …
```

**Dependency rule:** `ui` → `calculation` → `model`. Never reverse. `model` has no upward deps.

## Layer roles

**`model/`** — Physics evaluators AND optimizers. Both belong here, but in
**separate packages**: `model/cvm/` evaluates, `model/equilibrium/` solves.

`CVMGibbsModel` is a **pure evaluator**. It answers "what are G/H/S and their
derivatives at these system parameters (elements, structure, ECIs), macro
parameters (T, x), and micro parameters (u)?" — nothing more. It owns no
iteration and no convergence logic.

Solvers hold a model and drive it from outside. Both do this the same way, and
differ only in which unknowns they solve for:

| | `CvmNewtonSolver` | `HillertSolver` |
|---|---|---|
| composition | fixed constraint | an unknown |
| reads | `gmu` / `gmuu` (`ncf`) | `gmuFull` / `gmuuFull` (`ncf+K`) |
| phases | one | N |

**Do not put a solver loop back on `CVMGibbsModel`.** It was both evaluator and
optimizer until this split; the mixture is what made the G expressions hard to
audit.

**`calculation/`** — **Public API, discovery, and dispatch.**
1. **API**: `CalculationService.calculate`/`calculateScan` is the single named entry
   point (see below). CLI, GUI, and external JVM callers all bottom out here.
2. **Discovery**: `CalculationDescriptor.Registry` tells the GUI what can be
   calculated and which parameters a form needs.
3. **Statistics**: post-processing that must not live in the model layer
   (`MCSStatisticsProcessor` — τ_int, block averaging, jackknife Cv). Note this
   class exists but is not currently wired into `ThermodynamicWorkflow.runMcs`.

**`ui/`** — Collects inputs, dispatches, renders. No physics, no statistics.

---

## The calculation API contract

One entry point, used identically by CLI, GUI, and external callers — modelled on
pycalphad's `equilibrium(dbf, comps, phases, conditions)`:

```java
ThermodynamicResult r = service.calculate(session, conditions, property, sink, eventSink);
List<ThermodynamicResult> pts = service.calculateScan(session, scan, property, sink, eventSink);
```

- **Composition is named, never positional.** `Conditions(T, Map.of("Ti", 0.5))`.
  Any one element may be omitted and is derived to sum to 1. There is no
  "index 0 is the dependent species" convention any more — it caused a real
  silently-wrong-answer bug and was removed.
- **`ConditionsScan`** folds T and composition ranges into one object; at most one
  axis may vary per scan (enforced in the constructor).
- **MCS algorithm parameters** (`L`, `nEquil`, `nAvg`) go in `McsParams`, not in
  `Conditions` — they are not physical conditions.
- `execute(modelSpecs, jobSpecs, …)` still exists for the GUI's metadata-driven
  form, but delegates to `calculateScan`, so there is one code path underneath.
- **External (non-JVM) callers** use the `api` subcommand: JSON on stdin/stdout.
  See README for the schema.

**Two hazards this API guards against — do not weaken them:**

1. **ECI names are matched by string.** `CECEvaluator` silently leaves an unmatched
   term's interaction at `0.0` and still reports success. `EciValidator` rejects
   both unmatched names and incomplete coverage before any calculation runs, reusing
   `CECEvaluator`'s own alias rules so the two cannot diverge.
2. **CVM minimization can fail to converge** and still return plausible numbers
   (it hits a 20-iteration cap). `ThermodynamicResult.converged` carries the flag;
   the API reports it per point. Always check it.

---

## The session contract — read this carefully

`ModelSession` is the central object. It is **immutable** and holds pre-computed
state for one (elements, structure, model, engine) identity:

| Field | Type | Content |
|-------|------|---------|
| `systemId` | `Workspace.SystemId` | elements / structure / model |
| `cecEntry` | `CECEntry` | loaded Hamiltonian (ECI terms) |
| `resolvedHamiltonianId` | `String` | actual Hamiltonian ID used, or `<inline>` |
| `engineConfig` | `EngineConfig` | `CVM` or `MCS` |

**Built by** `ModelSession.Builder.build(systemId, engineConfig, progressSink)`:
1. Resolves the Hamiltonian ID — always the `_CVCF`-suffixed form
2. Loads and validates it from the Hamiltonian store

An overload `build(systemId, engineConfig, cecEntry, progressSink)` takes a
caller-supplied `CECEntry` and skips the store — used by the JSON API when an
external tool supplies its own ECIs. Sessions built that way must **not** go through
`CalculationService.getOrBuildSession`, whose cache is keyed only on
(systemId, engineConfig) and would collide with a stored-Hamiltonian session.

**What `build()` does NOT do:** it does not run cluster identification. Despite the
name, the expensive Stage 1–4 pipeline runs lazily inside
`CVMGibbsModel.initialize(...)` (called from `ThermodynamicWorkflow.runCvm`) and
`MCSGeometry.build(...)` for MCS. Caching is per-engine:
`ThermodynamicWorkflow` holds a `CvmCache` keyed on the session and an `McsCache`
keyed on (session, L, T).

**Passed as first arg** to all `CalculationService` and `ThermodynamicWorkflow`
calculation methods. Never null.

---

## Dataflow

**Type-1a — cluster identification** (`ClusterCFIdentificationPipeline.runFullWorkflow`):
Stage 1 geometric symmetry → Stage 2 algebraic CF orbits → Stage 3 orthogonal
C-matrix → Stage 4 CVCF transformation. Produces `PipelineResult`, held in memory.

**Type-2 CVM** — `ThermodynamicWorkflow.runCvm`:
```
CVMGibbsModel.of(...)              once per session (builds CvmGeometry, cached)
  → new CvmNewtonSolver(model)
      .solve(T, x, ...)            Newton–Raphson; returns Result.converged
       → model.at(T, x, u)         a State per trial point
          → CECEvaluator.evaluate()  eci[i] = a + b·T, matched by name
  → Result.state()                 the State at equilibrium u
     → st.gm() / st.hm() / st.sm() / st.g()   computed on demand
     → st.pairSroByShell()                    Cowley-Warren α
```

**Type-2 MCS** — `ThermodynamicWorkflow.runMcs`:
```
MCSGeometry.build(session, L)             cached per (session, L)
  → MCSRunner.forTemperature(geo, session, T)   ECI eval + transform, cached per T
     → MetropolisMC                             equilibration + averaging sweeps
```

---

## ID conventions

| ID | Formula | Example |
|----|---------|---------|
| `hamiltonianId` | `{elements}_{structure}_{model}` | `Nb-Ti_BCC_A2_T` |
| CVCF Hamiltonian | `{hamiltonianId}_CVCF` | `Nb-Ti_BCC_A2_T_CVCF` |
| `clusterId` | `{structure}_{model}_{ncomp}` | `BCC_A2_T_bin` |
| ncomp suffix | 2→`bin`, 3→`tern`, 4→`quat` | — |

`SystemId` derives these. Use `SystemId.hamiltonianId()` and `SystemId.clusterId()` — do not construct them by hand.

---

## CVM means CVCF only

`EngineConfig("CVM")` always uses the CVCF basis. There is no "ORTHO" mode for CVM. The string `"CVM"` in `EngineConfig.engineType` means CVCF-basis CVM.

---

## Input file naming convention

Cluster files: `inputs/clus/<structure>-<model>.txt` → e.g. `clus/BCC_A2-T.txt`
Symmetry files: `inputs/sym/<structure>-SG.txt` → e.g. `sym/BCC_A2-SG.txt`

`ModelSession.Builder` derives these paths automatically from `SystemId.structure` and `SystemId.model`. Do not hard-code file paths.

---

## GUI session lifecycle

1. `DynamicCalculationPanel` starts with defaults pre-filled (`Nb-Ti / BCC_A2 / T`).
2. `DocumentListener` on each identity field calls `context.setSystem(...)` on every
   keystroke → `WorkbenchContext` updated. A change to the elements combo also
   rebuilds the parameter form, so composition spinners can't go stale.
3. **Run Calculation** → `startExecution()` on a `SwingWorker`. There is no separate
   "Rebuild Session" step: `service.execute(...)` calls `getOrBuildSession()`
   internally, so the session is built (or reused from cache) on demand.
4. On completion, `done()` publishes the session via `context.setActiveSession(...)`
   and the result to the output panel.

**Thread rule:** session building and calculation are blocking and slow (disk I/O,
cluster identification, minimization). Always run on a `SwingWorker`, never on the
EDT. Only touch Swing state in `done()` or via `publish`/`process`.

---

## Workspace location

```
~/CEWorkbench/              (default, or ./data/CEWorkbench/ if it exists locally)
 ├─ inputs/clus/            cluster coordinate files
 ├─ inputs/sym/             symmetry group files
 └─ hamiltonians/<id>/hamiltonian.json
```

`new Workspace()` picks `./data/CEWorkbench/` if it exists, otherwise `~/CEWorkbench/`. Use `appCtx.getWorkspace()` to get paths — do not construct paths manually.

---

## Progress streaming pattern

All long-running operations accept `Consumer<String> progressSink` (text lines → `OutputPanel` log) and `Consumer<ProgressEvent> eventSink` (structured events → `ResultChartPanel` chart). Both may be null — always null-check before calling. Use the `emit(sink, msg)` helper pattern already present in each class.

```java
// correct
if (sink != null) sink.accept(msg);

// or the helper already in most classes
private static void emit(Consumer<String> sink, String msg) {
    if (sink != null) sink.accept(msg);
}
```

---

## What not to do

- Do not specify composition positionally. Use `Conditions` with element names —
  the positional convention caused a real silently-wrong-answer bug and was removed.
- Do not weaken `EciValidator` or duplicate `CECEvaluator`'s alias-matching rules.
  A second, divergent copy would let a name validate and then fail to map, silently
  zeroing that interaction.
- Do not report a CVM result without checking `converged` — a non-converged run
  returns plausible-looking numbers.
- Do not run cluster identification eagerly in `ModelSession.Builder` — it is
  deliberately lazy, inside `CVMGibbsModel.initialize` / `MCSGeometry.build`.
- Do not store mutable state in `ModelSession` — it is shared read-only across all
  scan points.
- Do not add calculation-layer imports to `ModelSession.Builder` or anything else in
  `model/` — the dependency rule is one-way.
- Do not route a caller-supplied-ECI session through
  `CalculationService.getOrBuildSession` — its cache key would collide with the
  stored-Hamiltonian session for the same system.
- Do not use `JScrollPane` as the top-level wrapper for `GridBagLayout` forms in
  Nimbus dark theme — the viewport background renders incorrectly. Use
  `add(buildForm(), BorderLayout.NORTH)` instead.
- Do not call `SwingWorker.get()` on the EDT outside of `done()`.
- Do not bypass `context.setSystem()` when changing system identity in a GUI panel —
  it is the only way to propagate changes and invalidate the session.

---

## Ternary isothermal-section plotting — why Python renders, Java computes

`TernaryGridScan` (Java, in-process, session-cached) sweeps a 2-D composition
grid over a 3-component system at fixed temperature and computes it directly
against `CalculationService.calculate` — not through `ConditionsScan`, which
supports only one varying axis at a time (temperature XOR one composition
element; see its class doc) and cannot express a full ternary sweep.

Rendering the result as a ternary contour is delegated to
`scripts/isothermal_section.py` (mpltern) rather than drawn in Java. This was
a deliberate choice, not a shortcut: no maintained Java library offers
ternary contour plotting — JFreeChart has no ternary axis support, and the
one abandoned point-plotting library found (`jTernaryPlot`, last updated
2013, unclear license) has no fill/contour capability. Plotly was also
checked and rejected: neither plotly.js nor Python plotly has a native
ternary contour trace either — Python's `figure_factory.create_ternary_contour`
is a convenience wrapper that does the same barycentric-to-Cartesian
transform-and-mask mpltern does internally, and that wrapper doesn't exist in
JS at all. mpltern is the only option that provides real ternary contour
support without hand-rolling the triangle geometry.

Java remains the sole source of truth for the physics: `TernaryPlotRenderer`
writes the already-computed grid to a temp JSON file and shells out to the
script's `--from-json` mode, which only turns numbers into pixels — it never
recomputes anything. This split (Java computes / Python renders) is shared
by both the GUI (`TernaryPlotPanel` → `OutputPanel`'s ternary card) and the
JSON API (`TernaryGridCommand`'s `"render":true`), so there is one rendering
code path, not two.

### Four composition regions, each handled differently

Every grid point falls into exactly one region, classified by how many of
its three mole fractions are exactly zero:

- **Interior** (all three &gt; 0) — a genuine ternary CVM solve via the
  session passed into `TernaryGridScan.run`. If it fails to converge, the
  point is skipped — no interpolation or synthetic data is substituted for
  it. The ternary CVM solver is known to be numerically fragile in the
  composition band adjacent to a binary edge (a near-edge Newton-Raphson
  instability, not a bug in the scan itself), so some interior points near
  an edge will legitimately be missing from the result; increasing `n`
  recovers coverage density since the failure is per-point, not systemic.
- **Edge** (exactly one = 0) — never evaluated on the ternary Hamiltonian at
  all. Swept as its own 1-D scan on a genuine **binary** CVM session.
  `BinarySubsystemExtractor` extracts the binary-pair CVCF terms
  (point/pair/triangle/tetrahedron terms whose CVCF name refers only to the
  two edge elements — e.g. `e21AB`, `e22AB`, `e3AB`, `e4AB`) from the ternary
  Hamiltonian's `cecTerms`, renamed to the binary basis's own A/B letter
  assignment, and builds a real 2-component `ModelSession` from them (cached
  per pair for the life of the scan). This is sound because binary-cluster
  CECs are inherited unchanged into any higher-order system containing that
  pair (Eq. 30 of the paper) — no transformation, only extraction. Evaluating
  the ternary Hamiltonian at a zero-composition edge was tried first and
  found to still be numerically fragile exactly where it's needed most (a
  large fraction of edge points failing to converge even there, for some
  systems); the true binary problem has a much smaller configuration space
  and solves far more reliably.
- **Corner** (two = 0, i.e. a pure element) — no calculation at all. G/H/S
  are trivially 0 by definition (nothing to mix with a single component).
  Pair SRO's Eq. 40 reference `x_P*x_R` is 0 there, making alpha
  mathematically undefined, so no point is emitted for SRO at a corner.

An earlier version of this scan tried to paper over interior near-edge
failures by linearly interpolating between the nearest edge value and a
converged interior point further along the same composition ray. This was
removed: for a near-**corner** interior point (two mole fractions both
small, but nonzero), "nearest edge" is ambiguous between two edges, and the
single-axis interpolation could snap the result almost entirely onto one
edge's value — producing a visibly wrong, discontinuous-looking point right
next to the correct interior contour. Skipping non-convergent points outright
is simpler and never wrong; only real solved values (interior, edge, or the
analytic corner) are ever reported.

### SRO in the ternary grid — why only pair SRO, and why some CFs are skipped

`TernaryGridScan.Quantity` is either a `PropertyQuantity` (G/H/S) or a
`PairSroQuantity` (1st-neighbour Cowley-Warren pair SRO for one unlike
species pair). Both route through the same single-point `calculate` call —
every CVM calculation already computes SRO as a side effect
(`ThermodynamicWorkflow.computeSro`), so `PairSroQuantity` just extracts a
different field from the same `ThermodynamicResult`, not a separate
calculation path.

Only pair SRO (1NN) is exposed. Extending this to triangle/tetrahedron
multi-site SRO looks straightforward at first — cluster probabilities for
*every* cluster type are already available via
`CVMGibbsModel.State`'s cluster variables — but most CVCF correlation
functions for triangle/tetrahedron clusters are not single physical
probabilities. `CvCfBasis.VSpec` makes this explicit: each CF is defined as
either `product(...)` (a single site-atom-pair probability, directly
SRO-eligible via `alpha = 1 − ρ/reference`) or `diff(...)`/`combo(...)` (a
signed linear combination of multiple probabilities, e.g. the ternary
binary-triangle CFs `v3AB = ρ^RPR − ρ^PRP`). A `diff`-type CF has no natural
`[0,1]` reference, so `1 − value/reference` isn't meaningful for it.

The literature confirms there's no shortcut around this: Goff, Li, Sinnott,
Dabo (PRB 104, 054109, 2021) — one of the two papers the Eq. 41 SRO code
cites — define one SRO-like parameter per distinct, symmetry-labeled
*occupation probability*, never on a signed difference of probabilities.
Nor does the "orthogonal" (Chebyshev/Inden-polynomial) basis used internally
by `CMatrixPipeline`/`ClusterCFIdentificationPipeline` help: those orthogonal
CFs are themselves signed polynomial moments of occupation, related to the
CVCF `diff`-type CFs by an invertible linear transform (`CvCfBasis`'s
`T`/`Tinv`) — swapping basis just trades one signed combination for another,
it doesn't produce a probability. If multi-site SRO is added later, the
`diff`-type CFs need to be split into their constituent `product(...)` terms
(each a real probability with a real reference) rather than computed
directly from the CF as defined.

**A pair-SRO scan only attempts that pair's own edge.** `PairSroQuantity`
(A,B) is undefined on the *other* two binary edges of a ternary system —
those edges each exclude one of A/B, so the pair being asked about doesn't
even exist there. `TernaryGridScan.run` checks `isPairEdge` before running
the edge's binary solve for an SRO request and, if it's the wrong edge,
skips the point without attempting a calculation or counting it in
`skipped` — exactly like a corner. This matters for interpreting `skipped`:
before this check existed, an SRO scan's `skipped` count included every
point on the two irrelevant edges (a large, systematic inflation that made
SRO look far less reliable than G/H/S on the *same* underlying CVM solves,
when in fact the interior/edge success rate is identical — verified by
diffing the exact composition sets that succeeded for G vs. SRO on the same
grid). After the fix, `skipped` means the same thing for every quantity:
a real computation was attempted and failed.

---

## Key files for context

| File | Why |
|------|-----|
| `model/ModelSession.java` | Session contract + Builder — most important class |
| `model/storage/Workspace.java` | `SystemId`, ID derivation, path layout, data-root resolution |
| `model/cvm/CVMGibbsModel.java` | Pure CVM evaluator; nested `State` (G/H/S, derivatives, SRO) |
| `model/cvm/CvmGeometry.java` | Immutable Stage 1–4 product; cluster algebra only |
| `model/cvm/CvCfBasis.java` | CVCF basis registry; expected ECI names per system |
| `model/equilibrium/CvmNewtonSolver.java` | Fixed-composition Newton–Raphson over a model |
| `model/equilibrium/HillertSolver.java` | Multi-phase equilibrium; whole method in one file |
| `model/equilibrium/SgteDatabase.java` | SGTE Unary v4.4 parser; G0m and its T-derivatives |
| `model/hamiltonian/CECEvaluator.java` | Name→basis ECI mapping and its alias rules |
| `model/mcs/MCSGeometry.java` · `MCSRunner.java` · `MetropolisMC.java` | MCS geometry, ECI setup, sweep loop |
| `calculation/Conditions.java` | Named-composition API entry type |
| `calculation/EciValidator.java` | Strict ECI validation |
| `calculation/workflow/CalculationService.java` | Public API (`calculate`/`calculateScan`) |
| `calculation/workflow/thermo/ThermodynamicWorkflow.java` | CVM/MCS dispatch, caching, SRO/convergence wiring |
| `ui/cli/ApiCommand.java` | JSON API for external callers |
| `ui/cli/TernaryGridCommand.java` | JSON API for ternary composition-grid scans (`ternary_grid` subcommand); optional `"render":true` returns a base64 PNG |
| `calculation/workflow/TernaryGridScan.java` | In-process ternary composition-grid sweep, near-edge interpolation via binary-edge fallback |
| `calculation/workflow/BinarySubsystemExtractor.java` | Extracts a binary CVCF sub-Hamiltonian for one pair from a ternary+ Hamiltonian's inherited terms |
| `calculation/workflow/TernaryPlotRenderer.java` | Shells out to `scripts/isothermal_section.py` (mpltern) to render a grid result to PNG |
| `ui/gui/DynamicCalculationPanel.java` | GUI calculation entry point |
| `ui/gui/TernaryPlotPanel.java` | GUI ternary isothermal-section panel (explorer column controls) |
| `ui/gui/WorkbenchContext.java` | GUI session state, listeners |
| `CEWorkbenchContext.java` | App wiring — how layers connect |

---

## Testing

There is **no test source tree** (`src/test/java` does not exist; only stale build
artifacts). Verification is by CLI invocation against known values:

```bash
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"                    # G = -3480.5209063901
./gradlew run --args="calc_min Nb-Ti-V BCC_A2 T 1000 Ti=0.33 V=0.34"          # G = -7051.1257304632
./gradlew run --args="calc_min Nb-Ti-V-Zr BCC_A2 T 1273 Ti=0.25 V=0.25 Zr=0.25 S"  # S = 11.0812146249
```

Re-run these after any change to the calculation path. If you add tests, the layer
split makes model-layer classes directly unit-testable with no mocks.

For MCS-specific changes (`model/mcs/**`), also run:

```bash
./gradlew runScratch -PscratchClass=org.ce.model.mcs.AnalyticConfigVerification
```

Cross-checks orthogonal CFs, CVCF CFs, energy, ΔE, and Metropolis trajectory against
configurations with independently-derivable ground truth (pure element, perfectly
ordered B2, perfectly random, composition boundaries) for K=2/3/4 — not just
internal self-consistency between MCS's own code paths. Expect `RESULT: PASS`.
This suite is what caught a real point-CF-column-ordering bug that silently broke
K≥3 CVCF energies while every internal-consistency check kept passing (see the
class doc for the full story) — internal-consistency-only checks are not sufficient
for this codebase's MCS correctness gate.

For CVM solver or evaluator changes (`model/cvm/**`, `model/equilibrium/**`), also run:

```bash
./gradlew runScratch -PscratchClass=org.ce.scratch.TernaryReferenceValidation
```

Validates both minimisers against a Mathematica `phaseq` reference point
(Mo-Nb-Ta / BCC_A2 / T, 1000 K, x=[0.33, 0.33, 0.34]): G/Gm/Hm/Sm to 1e-5
relative, and all 18 converged CFs to 1e-4. Expect `RESULT: PASS`.

This is an external check, not self-consistency: `CvmNewtonSolver` (fixed
composition, `ncf` system) and `HillertSolver` (composition as unknown, widened
`ncf+K` system with a Lagrange constraint and outer μ loop) share only the
evaluator, so agreement between them and with the reference is evidence about
`CVMGibbsModel` itself.

**The reference is stored as (ECI, CF) pairs, and checked as pairs.** The two
codes differ in cluster-algebra conventions — labels may be exchanged, and block
order differs (the reference emits the 1NN pair block before 2NN; we emit 2NN
first). Each code is internally consistent, but neither a position nor a label
identifies a cluster across them on its own. What must correspond is the cluster:
the ECI and the equilibrium CF have to be the same one. The gate resolves a name,
checks the ECI at that slot, and only then compares the CF — a CF checked alone
could pass while attached to the wrong cluster, if a label were exchanged
consistently in both the ECI input and the CF output.

Two conventions are reconciled, and they differ in kind: pair CFs are spelled
shell-last in the reference (`v2AB1`) and shell-first here (`v21AB`), which is a
spelling difference, not a reordering — the `v21`/`v22` blocks are **not**
transposed, though a positional read of an unlabelled vector makes them look as
though they are. Separately, `v3ABC1`/`v3ABC3` genuinely are exchanged; both carry
ECI 0, so the ECI cannot discriminate there and that one mapping rests on the CF
cross-comparison alone.

μ is not compared entry-by-entry: for a single phase it is underdetermined
(one Gibbs-Duhem equation, K unknowns), so the gate asserts `Σμᵢxᵢ = G` instead.
`CvmNewtonSolver` has no μ at all — composition is fixed there.

This covers one interior composition. The near-edge band (e.g. x=[0.05, 0.05, 0.90]
on the same system) is where **both** solvers still fail — a separate open item,
not covered here.
