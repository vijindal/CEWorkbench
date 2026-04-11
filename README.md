# CE Thermodynamics Workbench

A scientific software framework for **Cluster Expansion (CE) based thermodynamic calculations** and **cluster identification**. The workbench provides tools to identify cluster basis functions, manage effective cluster interactions (ECI), and compute thermodynamic equilibrium states for alloy systems.

## Overview

Three classes of work are supported:

| Type | Name | Description |
|------|------|-------------|
| **1a** | Cluster Identification | Load ordered/disordered cluster files and symmetry groups. Four stages: geometric identification (1–2), orthogonal basis / C-matrix construction (3), CVCF transformation (4). |
| **1b** | Hamiltonian Scaffold | Auto-generate an empty ECI (Hamiltonian) JSON file from cluster identification results, ready for editing. |
| **2** | Thermodynamic Equilibrium | Minimize free energy with CVM (Newton–Raphson, CVCF basis) or Monte Carlo (MCS) to produce G / H / S at given T and composition. |

---

## Quick Start

### Prerequisites

- Java 21 or later
- Gradle 9.3+

### Launch the GUI

```bash
./gradlew runGui
```

Opens the VS Code-style dark workbench. Use the activity bar on the left to switch between the three panels.

### Run the CLI

```bash
# Full pipeline — explicit system
./gradlew run --args="all Nb-Ti BCC_A2 T"

# Individual stages
./gradlew run --args="type1a Nb-Ti BCC_A2 T"
./gradlew run --args="type1b Nb-Ti BCC_A2 T"
./gradlew run --args="type2  Nb-Ti BCC_A2 T"

# Single-point CVM with minimisation
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 0.5"

# Evaluate at fixed correlation functions (no minimisation)
./gradlew run --args="calc_fixed Nb-Ti BCC_A2 T 1000 0.5 0.5 0.1 0.2 0.3"

# View loaded Hamiltonian
./gradlew run --args="view Nb-Ti BCC_A2 T"

# Verbose output for any mode
./gradlew run --args="type2 Nb-Ti BCC_A2 T --verbose"
```

**Mode reference:**

| Mode | Description |
|------|-------------|
| `type1a` | Cluster identification only |
| `type1b` | Hamiltonian scaffold only |
| `type2` | Thermodynamic calculation (temperature scan) |
| `all` | Runs type1a → type1b → type2 in sequence |
| `calc_min` | Single-point CVM with Newton–Raphson minimisation |
| `calc_fixed` | Evaluate G/H/S at a user-supplied CF vector |
| `view` | Print Hamiltonian ECI table to stdout |

### Build

```bash
./gradlew build
```

---

## Architecture

Three explicit layers with strict one-way dependencies, plus a root-level application context:

```
org.ce
 ├─ CEWorkbench.java          GUI entry point
 ├─ CEWorkbenchContext.java   Wires all layers; shared by GUI and CLI
 │
 ├─ model/       Physics models, cluster data, Hamiltonian, storage I/O
 ├─ calculation/ Engine algorithms and workflow orchestration
 └─ ui/          GUI (Swing) and CLI interfaces
```

### Model Layer — `org.ce.model`

Static, persistent physics evaluators. Built once from system info (elements, structure, model, ECI) and live for the duration of a calculation. The calculation layer calls into them with (T, composition, CFs) and gets back thermodynamic quantities. No dependency on `calculation` or `ui`.

| Package | Key Classes | Role |
|---------|-------------|------|
| `model` | `ModelSession`, `EngineConfig` | Immutable session holding pre-computed cluster data, Hamiltonian, and CVCF basis for one system identity |
| `model.cluster` | `ClusterIdentifier`, `CFIdentifier`, `CMatrix`, `AllClusterData`, `ClusterVariableEvaluator`, `ClusterPrimitives` | Cluster geometry, CF basis construction, C-matrix, correlation function evaluation |
| `model.cluster.cvcf` | `CvCfBasis`, `CvCfBasisTransformer`, `CvCfMatrixGenerator`, `CvCfDefinition` | CVCF basis and transformation machinery |
| `model.cvm` | `CVMGibbsModel`, `CVMSolver` | CVM Gibbs functional evaluator + Newton–Raphson minimizer (NR loop, convergence logic) |
| `model.mcs` | `LatticeConfig`, `EmbeddingData`, `Embedding`, `EmbeddingGenerator`, `SiteOperatorBasis`, `CvCfEvaluator`, `Vector3D`, `MCEngine`, `MCSampler`, `MCSRunner`, `ExchangeStep`, `LocalEnergyCalc`, `MCResult` | MCS supercell state, physics evaluators, algorithm drivers (equilibration+averaging loops, statistics, Metropolis step) |
| `model.hamiltonian` | `CECEntry`, `CECTerm`, `CECEvaluator` | ECI model objects; `CECEvaluator` computes `eci[i] = a + b·T` |
| `model` | `ThermodynamicInput` | Bundled input to CVM/MCS optimizers |
| `model.result` | `ThermodynamicResult` | Output value object (unified; no EquilibriumState wrapper) |
| `model.storage` | `Workspace`, `InputLoader`, `HamiltonianStore`, `DataStore` | All disk I/O — path resolution, cluster file parsing, JSON persistence |

**`ModelSession` is the session contract.** It is built once per (elements, structure, model) identity by `ModelSession.Builder`, which runs cluster identification, resolves the Hamiltonian, and looks up the CVCF basis. All calculations within a session reuse this pre-computed state — no redundant disk reads or re-identification.

### Calculation Layer — `org.ce.calculation`

Active algorithm drivers. Accepts a `ModelSession`; drives algorithms using model-layer physics evaluators. Performs no cluster identification or Hamiltonian loading.

| Package | Key Classes | Role |
|---------|-------------|------|
| `calculation` | `CalculationDescriptor`, `CalculationSpecifications`, `CalculationRegistry` | **Discovery-based Metadata**: Defines Property, Mode, and Parameter vocabularies; provides available options for a given engine. |
| `calculation.workflow` | `CalculationService` | **Unified Entry Point**: Orchestrates model construction (`execute`) and result dispatching (`executeScan`). |
| `calculation.workflow.cec` | `CECManagementWorkflow` | Type-1b: scaffold, load, validate, save Hamiltonian. |
| `calculation.workflow.thermo` | `ThermodynamicWorkflow`, `LineScanWorkflow`, `GridScanWorkflow`, `FiniteSizeScanWorkflow` | Type-2 internal implementation: single-point, temperature/composition/grid scans, finite-size scaling. |

### UI Layer — `org.ce.ui`

| Package | Key Classes | Role |
|---------|-------------|------|
| `ui.cli` | `Main` | CLI entry point; argument parsing for all modes |
| `ui.gui` | `MainWindow`, `WorkbenchContext` | Window frame; shared observable session state |
| `ui.gui` | `CalculationPanel` | System identity fields, session rebuild, CVM/MCS calculation |
| `ui.gui` | `DataPreparationPanel` | Type-1a inputs: cluster files, symmetry groups, component count |
| `ui.gui` | `CECManagementPanel` | Type-1b: scaffold / load / edit / save Hamiltonian ECI table |
| `ui.gui` | `OutputPanel`, `ResultsPanel` | Always-visible results display and scrollable log |
| `ui.gui` | `ActivityBar`, `ExplorerPanel`, `HeaderBar`, `StatusBar` | Chrome and navigation |

---

## GUI Layout

```
JFrame (BorderLayout)
 ├─ NORTH  → HeaderBar        — app name + active system identity
 ├─ CENTER → JSplitPane
 │            ├─ LEFT  → ActivityBar + ExplorerPanel  — CardLayout; one parameter panel at a time
 │            └─ RIGHT → OutputPanel                  — RESULTS (top) / LOG (bottom)
 └─ SOUTH  → StatusBar        — one-line operation status
```

**Session lifecycle in the GUI:**

1. Open the app → Calculation panel is shown with pre-filled defaults (`Nb-Ti / BCC_A2 / T`).
2. Edit elements / structure / model as needed → fields sync to `WorkbenchContext`.
3. Click **Rebuild Session** → `ModelSession.Builder.build()` runs on a background thread (cluster identification + Hamiltonian load).
4. Session status turns teal → **Calculate** button becomes enabled.
5. Set temperature and composition → click **Calculate** → results stream to the output panel in real time.

Alternatively, complete Type-1a (Data Prep panel) or load a Hamiltonian (Hamiltonian panel) — those panels also trigger a session build on success, so the Calculation panel becomes ready automatically.

---

## Calculation Workflows

### Type-1a — Cluster Identification

```
Input files (clus/*.txt + sym/*.txt)
  ↓
ClusterIdentificationWorkflow.identify(request)
  Stage 1 — ClusterIdentifier     (geometric symmetry)
  Stage 2 — CFIdentifier          (algebraic CF basis)
  Stage 3 — CMatrix               (orthogonal basis / C-matrix)
  Stage 4 — CvCfBasisTransformer  (CVCF transformation)
  ↓
AllClusterData  (held in ModelSession.clusterData)
```

### Type-1b — Hamiltonian Scaffold

```
CECManagementWorkflow.scaffoldFromClusterData(hamiltonianId, ...)
  ↓
hamiltonian.json  (all ECI terms a=0, b=0 — edit and save)
~/CEWorkbench/hamiltonians/<hamiltonianId>/hamiltonian.json
```

### Type-2 — Thermodynamic Calculation (CVM)

```
ModelSession (pre-built: AllClusterData + CECEntry + CvCfBasis)
  ↓
ThermodynamicWorkflow.runCalculation(session, request)
  ↓
CVMEngine.compute(input)                          ← calculation layer
  ├─ CECEvaluator   (eci[i] = a + b·T)
  ├─ CVMGibbsModel  (model layer: holds cluster geometry + ECI)
  │   └─ .evaluate(T, x, CFs) → G, H, S, ∇G, ∇²G
  └─ CVMSolver      (Newton–Raphson, tolerance 1e-5, max 400 iter)
  ↓
ThermodynamicResult { gibbsEnergy, enthalpy, entropy, temperature, composition }
```

---

## Package Structure

```
src/main/java/org/ce
├─ CEWorkbench.java
├─ CEWorkbenchContext.java
│
├─ model/
│  ├─ ThermodynamicInput.java    (bundled input to optimizers)
│  ├─ EngineConfig.java
│  ├─ ModelSession.java          (+ nested Builder)
│  ├─ cluster/
│  │  ├─ AllClusterData.java
│  │  ├─ ClusterIdentifier.java
│  │  ├─ CFIdentifier.java
│  │  ├─ CMatrix.java
│  │  ├─ ClusterVariableEvaluator.java
│  │  ├─ ClusterPrimitives.java
│  │  └─ (15+ supporting classes)
│  │  └─ cvcf/
│  │     ├─ CvCfBasis.java
│  │     ├─ CvCfBasisTransformer.java
│  │     ├─ CvCfMatrixGenerator.java
│  │     └─ CvCfDefinition.java
│  ├─ cvm/
│  │  ├─ CVMGibbsModel.java      (Gibbs functional; evaluates G/H/S/gradients)
│  │  └─ CVMSolver.java          (Newton–Raphson loop; NR minimizer)
│  ├─ mcs/
│  │  ├─ LatticeConfig.java      (mutable supercell occupation array)
│  │  ├─ EmbeddingData.java      (pre-computed site→embedding map)
│  │  ├─ Embedding.java          (one cluster instance on lattice sites)
│  │  ├─ EmbeddingGenerator.java (builds EmbeddingData from cluster data)
│  │  ├─ SiteOperatorBasis.java  (evaluates φ_α(σ) basis functions)
│  │  ├─ CvCfEvaluator.java      (measures CVCF variables from LatticeConfig)
│  │  ├─ Vector3D.java           (3D vector for lattice arithmetic)
│  │  ├─ MCEngine.java           (equilibration+averaging sweep loop)
│  │  ├─ MCSampler.java          (accumulates statistics; tau_int, block avg, jackknife)
│  │  ├─ MCSRunner.java          (geometry setup; ECI transform)
│  │  ├─ ExchangeStep.java       (Metropolis exchange trial)
│  │  ├─ LocalEnergyCalc.java    (ΔE computation)
│  │  └─ MCResult.java           (immutable result from completed MCS run)
│  ├─ hamiltonian/
│  │  ├─ CECEntry.java
│  │  ├─ CECTerm.java
│  │  ├─ CECEvaluator.java
│  │  └─ Hamiltonian.java
│  ├─ result/
│  │  ├─ ThermodynamicResult.java
│  │  └─ EquilibriumState.java
│  └─ storage/
│     ├─ Workspace.java
│     ├─ InputLoader.java
│     ├─ HamiltonianStore.java
│     └─ DataStore.java
│
├─ calculation/
│  ├─ CalculationDescriptor.java
│  ├─ CalculationSpecifications.java
│  ├─ CalculationRegistry.java
│  └─ workflow/
│     ├─ CalculationService.java
│     ├─ CECManagementWorkflow.java (Hamiltonian scaffold/load/save)
│     └─ thermo/
│        ├─ ThermodynamicWorkflow.java
│        ├─ LineScanWorkflow.java
│        ├─ GridScanWorkflow.java
│        └─ FiniteSizeScanWorkflow.java
│
└─ ui/
   ├─ cli/
   │  └─ Main.java
   └─ gui/
      ├─ MainWindow.java
      ├─ WorkbenchContext.java
      ├─ ActivityBar.java
      ├─ ExplorerPanel.java
      ├─ OutputPanel.java
      ├─ ResultsPanel.java
      ├─ HeaderBar.java
      ├─ StatusBar.java
      ├─ DataPreparationPanel.java
      ├─ CECManagementPanel.java
      └─ CalculationPanel.java
```

---

## Input Data Files

Located in `~/CEWorkbench/inputs/` (runtime workspace). Defaults are bundled in `src/main/resources/` and copied on first run.

| File | Description |
|------|-------------|
| `clus/BCC_A2-T.txt` | BCC A2 (disordered) cluster coordinates |
| `clus/BCC_B2-T.txt` | BCC B2 (ordered) cluster coordinates |
| `clus/FCC_A1-TO.txt` | FCC A1 cluster coordinates |
| `sym/BCC_A2-SG.txt` | BCC A2 space group symmetry operations |
| `sym/BCC_B2-SG.txt` | BCC B2 space group symmetry operations |
| `sym/FCC_A1-SG.txt` | FCC A1 space group symmetry operations |
| `sym/HCP_A3-SG.txt` | HCP A3 space group symmetry operations |

Cluster files follow the naming convention `<structure>-<model>.txt`; symmetry files follow `<structure>-SG.txt`. The `ModelSession.Builder` derives file paths automatically from the system identity.

---

## Workspace Layout

All persistent data is stored under `~/CEWorkbench/`:

```
~/CEWorkbench/
 ├─ inputs/
 │   ├─ clus/      cluster coordinate files (*.txt)
 │   └─ sym/       symmetry group files (*.txt)
 └─ hamiltonians/
     └─ <hamiltonianId>/
         └─ hamiltonian.json    edit a and b ECI values here
```

`SystemId` derives IDs deterministically from elements / structure / model:
- `hamiltonianId` = `<elements>_<structure>_<model>` e.g. `Nb-Ti_BCC_A2_T`
- For CVM, the builder prefers a `_CVCF`-suffixed Hamiltonian (`Nb-Ti_BCC_A2_T_CVCF`) when available.

---

## Project Status

**Complete:**
- Three-layer architecture (Model / Calculation / UI) with one-way dependencies
- `ModelSession` — immutable session object; cluster identification + Hamiltonian load run once per system identity, reused across all scan points
- Cluster identification — Stages 1–4 (cluster, CF, C-matrix, CVCF transformation)
- CVM thermodynamic engine with Newton–Raphson minimisation
- MCS thermodynamic engine with finite-size scaling
- Hamiltonian scaffold, load, edit, save workflow
- VS Code dark GUI — self-contained Calculation panel, Activity Bar, Explorer, Output panel
- Bidirectional system identity sync across all GUI panels via `WorkbenchContext`
- CLI with all modes: `type1a`, `type1b`, `type2`, `all`, `calc_min`, `calc_fixed`, `view`

**Planned:**
- Additional symmetry groups and crystal structures
- Phase diagram grid scan visualization
- Export results to CSV / JSON

---

## License

[Add your license here]
