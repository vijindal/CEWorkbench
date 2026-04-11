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
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 0.5 --verbose"

# CLI — temperature scan
./gradlew run --args="type2 Nb-Ti BCC_A2 T --verbose"

# Compile only
./gradlew compileJava
```

All Gradle tasks: `runGui`, `runGuiDebug`, `run` (CLI), `runDebugCli`, `build`.

---

## Package layout

```
org.ce
├─ CEWorkbench.java          GUI main()
├─ CEWorkbenchContext.java   App-level wiring; shared by GUI and CLI
│
├─ model/                    Physics models, persistent state, disk I/O
│   ├─ ModelSession.java     Immutable session (see below — most important class)
│   ├─ EngineConfig.java     "CVM" or "MCS"
│   ├─ cluster/              Cluster geometry, CF basis, C-matrix, CVCF
│   │   └─ cvcf/             CVCF basis and transformation
│   ├─ cvm/                  CVM Gibbs functional + Newton-Raphson minimizer
│   │   ├─ CVMGibbsModel     Physics evaluator; evaluates G/H/S/gradients given (T, x, CFs)
│   │   └─ CVMSolver         Algorithm driver; owns NR loop, convergence logic
│   ├─ mcs/                  MCS supercell state + physics evaluators + algorithm drivers
│   │   ├─ LatticeConfig     Atomic occupation array; mutable supercell state
│   │   ├─ EmbeddingData     Pre-computed site→embedding map for supercell
│   │   ├─ Embedding         One cluster instance mapped onto lattice sites
│   │   ├─ EmbeddingGenerator Builds EmbeddingData from cluster data + lattice positions
│   │   ├─ SiteOperatorBasis  Evaluates orthogonal basis functions φ_α(σ)
│   │   ├─ CvCfEvaluator     Measures CVCF cluster variables from a LatticeConfig
│   │   ├─ Vector3D          3D vector for lattice position arithmetic
│   │   ├─ MCEngine          Algorithm driver; owns equilibration+averaging sweep loop
│   │   ├─ MCSampler         Algorithm driver; owns tau_int, block averaging, jackknife Cv
│   │   ├─ MCSRunner         Algorithm driver; owns BCC geometry, ECI transform
│   │   ├─ ExchangeStep      Algorithm driver; owns Metropolis accept/reject
│   │   ├─ LocalEnergyCalc   Physics evaluator; computes ΔE (static utility)
│   │   └─ MCResult          Immutable result DTO
│   ├─ hamiltonian/          CECEntry, CECTerm, CECEvaluator
│   ├─ result/               ThermodynamicResult, EquilibriumState
│   └─ storage/              Workspace, InputLoader, HamiltonianStore, DataStore
│
├─ calculation/              Metadata providers and unified calculation dispatcher
│   ├─ CalculationDescriptor.java   Vocabulary: Property, Mode, Parameter
│   ├─ CalculationSpecifications.java Value object for calculation parameters
│   ├─ CalculationRegistry.java     Discovery: provides available options per engine
│   └─ workflow/             Unified CalculationService + specialized scans
│       ├─ CalculationService.java  Unified entry point (execute/executeScan)
│       └─ thermo/           ThermodynamicWorkflow, LineScanWorkflow, etc.
│
└─ ui/
    ├─ cli/   Main.java      Redesigned to use unified service entry points
    └─ gui/   MainWindow, DynamicCalculationPanel (unified), OutputPanel...
```

**Dependency rule:** `ui` → `calculation` → `model`. Never reverse. `model` has no upward deps.

## Layer roles

**`model/`** — Physics evaluators AND optimizers. Evaluators (e.g., `CVMGibbsModel`) are queried for properties. Optimizers (e.g., `CVMSolver`, `MCEngine`) own algorithm loops and convergence logic. Both belong in the model layer.

**`calculation/`** — **Discovery and Dispatch Layer.**
1.  **Discovery**: Provides `CalculationRegistry` which the UI uses to "discover" what can be calculated and what parameters are needed.
2.  **Dispatch**: `CalculationService` acts as the unified facade. It manages `ModelSession` construction ("Model Construction Role") and routes specifications to model-layer optimizers ("Execution Role").

**`ui/`** — Specification-driven. Collects inputs into `ModelSpecifications` and `CalculationSpecifications`, then dispatches via the unified service.

---

## The session contract — read this carefully

`ModelSession` is the central object. It is **immutable** and holds everything pre-computed for one (elements, structure, model) identity:

| Field | Type | Content |
|-------|------|---------|
| `systemId` | `Workspace.SystemId` | elements / structure / model |
| `clusterData` | `AllClusterData` | clusters, CF basis, C-matrix |
| `cecEntry` | `CECEntry` | loaded Hamiltonian (ECI terms) |
| `resolvedHamiltonianId` | `String` | actual Hamiltonian ID used |
| `cvcfBasis` | `CvCfBasis` | CVCF basis for this system |
| `engineConfig` | `EngineConfig` | "CVM" or "MCS" |

**Built once** by `ModelSession.Builder.build(systemId, engineConfig, progressSink)`:
1. Runs `ClusterIdentificationWorkflow.identify()` — cluster files derived automatically as `clus/<structure>-<model>.txt` and `<structure>-SG`
2. Resolves Hamiltonian ID — for CVM prefers `<id>_CVCF`, falls back to base ID
3. Loads Hamiltonian via `CECManagementWorkflow.loadAndValidateCEC()`
4. Looks up `CvCfBasis` from registry

**Passed as first arg** to all `CalculationService` and `ThermodynamicWorkflow` methods. Never null — `ThermodynamicInput` constructor throws if `clusterData` is null.

**Do not** re-run cluster identification inside engines. `CVMEngine` constructs `CVMGibbsModel` directly from `input.clusterData` — it does NOT call `ClusterIdentificationWorkflow`.

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

1. `CalculationPanel` starts with defaults pre-filled (`Nb-Ti / BCC_A2 / T`).
2. `DocumentListener` on each field calls `context.setSystem(...)` on every keystroke → `WorkbenchContext` updated → session invalidated.
3. **Rebuild Session** → `ModelSession.Builder.build()` on a `SwingWorker` background thread → `context.setActiveSession(session)` on EDT.
4. `Calculate` button enabled only when `context.hasActiveSession()`.
5. Any change to system identity in any panel calls `context.setSystem()` → session invalidated → Calculate disabled until rebuilt.

**Thread rule:** `ModelSession.Builder.build()` is blocking and slow (disk I/O + cluster identification). Always run it on a `SwingWorker`, never on the EDT.

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

- Do not call `AllClusterData.identify()` inside optimizers or workflow classes — it belongs only in `ModelSession.Builder.build()`.
- Do not add fields to `ThermodynamicRequest` for system identity — it holds only calculation parameters (T, composition, MCS params, sinks). Session state lives in `ModelSession`.
- Do not store mutable state in `ModelSession` — it is shared read-only across all scan points.
- Do not add calculation-layer classes to `ModelSession.Builder` — it should only import from model layer (after move: no calculation-layer imports).
- Do not use `JScrollPane` as the top-level wrapper for `GridBagLayout` forms in Nimbus dark theme — the viewport background renders incorrectly. Use `add(buildForm(), BorderLayout.NORTH)` instead.
- Do not call `SwingWorker.get()` on the EDT outside of `done()`.
- Do not bypass `context.setSystem()` when changing system identity in a GUI panel — it is the only way to propagate changes and invalidate the session.

---

## Key files for context

| File | Why |
|------|-----|
| `model/ModelSession.java` | Session contract + Builder — most important class |
| `model/storage/Workspace.java` | ID derivation, path layout |
| `model/cvm/CVMGibbsModel.java` | CVM physics evaluator — evaluates G/H/S/gradients given (T, x, CFs) |
| `model/cvm/CVMSolver.java` | CVM optimizer — owns Newton–Raphson loop, convergence logic |
| `model/mcs/MCEngine.java` | MCS optimizer — owns equilibration+averaging sweep loop |
| `model/mcs/MCSampler.java` | MCS statistics — owns tau_int, block averaging, jackknife Cv |
| `calculation/workflow/thermo/ThermodynamicWorkflow.java` | Thin dispatcher — inlines CVM/MCS dispatch to model-layer optimizers |
| `calculation/workflow/CalculationService.java` | Public API used by GUI and CLI |
| `ui/gui/WorkbenchContext.java` | GUI session state, listeners |
| `ui/gui/CalculationPanel.java` | GUI calculation entry point |
| `CEWorkbenchContext.java` | App wiring — how layers connect |
