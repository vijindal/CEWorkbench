# CE Thermodynamics Workbench

A scientific software framework for **Cluster Expansion (CE) based thermodynamic calculations** and **cluster identification**. This workbench provides tools to generate reusable scientific data and compute thermodynamic equilibrium states for alloy systems.

## Overview

The CE Thermodynamics Workbench implements two fundamental classes of work:

1. **Type-1a — Cluster identification:** Load ordered/disordered cluster files and symmetry groups, identify correlation functions, and build the C-matrix.
2. **Type-1b — Hamiltonian scaffold:** Auto-generate an empty ECI (Hamiltonian) JSON file from saved cluster data, ready for editing.
3. **Type-2 — Thermodynamic equilibrium:** Minimize free energy with CVM or run Monte Carlo simulations (MCS) to compute G, H, S as a function of temperature and composition.

---

## Quick Start

### Prerequisites

- Java 25 or later
- Gradle 9.3+

### Launch the GUI

```bash
./gradlew runGui
```

Opens the VS Code-style dark workbench with three panels:

| Panel | Purpose |
|-------|---------|
| **Activity bar** (left strip) | Switch between Data Prep / Hamiltonian / Calculate |
| **Explorer** (side panel) | Parameters for the selected task |
| **Output** (main area) | Results display + log output |

### Run the CLI

```bash
# Full pipeline — built-in defaults (A-B / BCC_B2 / T)
./gradlew run

# Full pipeline — explicit system
./gradlew run --args="all Nb-Ti BCC_A2 T"

# Single mode
./gradlew run --args="type1a"
./gradlew run --args="type1b Nb-Ti BCC_A2 T"
./gradlew run --args="type2  Nb-Ti BCC_A2 T"
```

**Argument signature:** `[mode] [elements] [structure] [model]`

| Argument | Values | Default |
|----------|--------|---------|
| mode | `type1a` \| `type1b` \| `type2` \| `all` | `all` |
| elements | e.g. `Nb-Ti`, `A-B` | `A-B` |
| structure | e.g. `BCC_A2`, `BCC_B2` | `BCC_B2` |
| model | e.g. `T` | `T` |

### Build & Test

```bash
./gradlew build
./gradlew test
```

---

## Architecture

Four clean layers with strict one-way dependencies:

```
org.ce
 ├─ domain      → physics models and algorithms  (no upward deps)
 ├─ storage     → disk IO and persistence         (→ domain)
 ├─ workflow    → calculation orchestration       (→ domain, storage)
 └─ ui          → GUI and CLI interfaces          (→ workflow)
```

### Domain Layer

Pure physics and algorithms. No dependency on any other layer.

| Package | Contents |
|---------|----------|
| `domain/cluster` | `ClusterIdentifier`, `CFIdentifier`, `CMatrixBuilder`, `AllClusterData`, `SpaceGroup`, `Vector3D`, `Position`, `Site`, `Sublattice`, + 30 supporting classes |
| `domain/engine` | `ThermodynamicEngine` interface, `CVMEngine` (Cluster Variation Method), `MCSEngine` (Monte Carlo) |
| `domain/hamiltonian` | `Hamiltonian`, `CECEntry`, `CECTerm` — ECI model objects |
| `domain/result` | `ThermodynamicResult`, `EquilibriumState` |

### Storage Layer

All disk IO and data persistence.

| Class | Role |
|-------|------|
| `Workspace` | Resolves paths under `~/CEWorkbench/` |
| `InputLoader` | Parses `.txt` cluster and symmetry-group files |
| `ClusterDataStore` | Save/load `AllClusterData` as JSON |
| `HamiltonianStore` | Save/load `CECEntry` (Hamiltonian) as JSON |
| `SystemId` | Derives canonical `clusterId` and `hamiltonianId` from elements/structure/model |

### Workflow Layer

Orchestrates domain operations without implementing physics.

| Class | Role |
|-------|------|
| `ClusterIdentificationWorkflow` | Type-1a: cluster + CF identification |
| `ClusterIdentificationRequest` | Builder-pattern config for Type-1a |
| `CECManagementWorkflow` | Type-1b: scaffold, load, validate, update Hamiltonian |
| `ThermodynamicWorkflow` | Type-2: single-point, temperature scan, composition scan |
| `CalculationService` | High-level façade used by both GUI and CLI |
| `LineScanWorkflow`, `GridScanWorkflow` | Batched parameter sweeps |

### UI Layer

| Package | Contents |
|---------|----------|
| `ui/cli` | `Main.java` — entry point with positional argument parsing |
| `ui/gui` | VS Code-style dark workbench (see below) |

---

## GUI Components

```
JFrame (BorderLayout)
 ├─ NORTH  → HeaderBar          — app name + active system (elements · structure · model)
 ├─ WEST   → ActivityBar        — icon strip; switches ExplorerPanel cards
 ├─ CENTER → JSplitPane
 │            ├─ LEFT  → ExplorerPanel   — CardLayout; one parameter panel at a time
 │            └─ RIGHT → OutputPanel     — JSplitPane; RESULTS top / LOG bottom
 └─ SOUTH  → StatusBar          — one-line operation status
```

| Component | Description |
|-----------|-------------|
| `ActivityBar` | 52 px strip, 3 items (Data Prep / Hamiltonian / Calculate), blue left-accent on active |
| `ExplorerPanel` | CardLayout holder for the 3 parameter panels |
| `DataPreparationPanel` | Type-1a inputs: ordered + disordered cluster files, symmetry groups, component count |
| `CECManagementPanel` | Type-1b: scaffold / load / edit / save Hamiltonian ECI table |
| `CalculationPanel` | Type-2: elements, structure, model, T, x_B, engine selector, Calculate button |
| `OutputPanel` | Always-visible split pane: RESULTS rows (G/H/S/T/x_B/engine) + scrollable log |
| `HeaderBar` | Dark top bar; auto-updates from `WorkbenchContext` |
| `StatusBar` | Blue footer bar; one-line status from panels |
| `WorkbenchContext` | Shared observable session state — system identity propagated across all panels |

---

## Calculation Workflows

### Type-1a — Cluster Identification

```
Input files (clus/*.txt + sym/*.txt)
  ↓
ClusterIdentificationWorkflow.identify(config)
  ↓  Stage 1 — ClusterIdentifier
  ↓  Stage 2 — CFIdentifier
  ↓  Stage 3 — CMatrixBuilder
  ↓
AllClusterData  →  ClusterDataStore.save(clusterId, data)
                   ~/CEWorkbench/cluster-data/<clusterId>/cluster_data.json
```

### Type-1b — Hamiltonian Scaffold

```
ClusterDataStore.load(clusterId)  ← uses ncf from disordered CF result
  ↓
CECManagementWorkflow.scaffoldFromClusterData(...)
  ↓
hamiltonian.json  (all ECI terms set to 0 — edit a and b values)
~/CEWorkbench/hamiltonians/<hamiltonianId>/hamiltonian.json
```

### Type-2 — Thermodynamic Calculation

```
ClusterDataStore.load(clusterId)
HamiltonianStore.load(hamiltonianId)
  ↓
ThermodynamicEngine (CVMEngine or MCSEngine)
  ↓
ThermodynamicResult { gibbsEnergy, enthalpy, entropy, temperature, composition }
```

---

## Package Structure

```
app/src/main/java/org/ce
├─ domain
│  ├─ cluster
│  │  ├─ ClusterIdentifier.java
│  │  ├─ CFIdentifier.java
│  │  ├─ CMatrixBuilder.java
│  │  ├─ AllClusterData.java
│  │  ├─ Vector3D.java
│  │  └─ (30+ supporting classes)
│  ├─ engine
│  │  ├─ ThermodynamicEngine.java
│  │  ├─ CVMEngine.java
│  │  ├─ MCSEngine.java
│  │  └─ cvm/  (CVMFreeEnergy, CVMPhaseModel, NewtonRaphsonSolverSimple)
│  ├─ hamiltonian
│  │  ├─ Hamiltonian.java
│  │  ├─ CECEntry.java
│  │  └─ CECTerm.java
│  └─ result
│     ├─ ThermodynamicResult.java
│     └─ EquilibriumState.java
│
├─ storage
│  ├─ SystemId.java
│  ├─ Workspace.java
│  ├─ InputLoader.java
│  ├─ ClusterDataStore.java
│  └─ HamiltonianStore.java
│
├─ workflow
│  ├─ CalculationService.java
│  ├─ ClusterIdentificationWorkflow.java
│  ├─ ClusterIdentificationRequest.java
│  ├─ cec/
│  │  ├─ CECManagementWorkflow.java
│  │  └─ CFMetadata.java
│  └─ thermo/
│     ├─ ThermodynamicWorkflow.java
│     ├─ LineScanWorkflow.java
│     ├─ GridScanWorkflow.java
│     ├─ ThermodynamicRequest.java
│     └─ ThermodynamicData.java
│
└─ ui
   ├─ cli/
   │  └─ Main.java
   └─ gui/
      ├─ MainWindow.java
      ├─ WorkbenchContext.java
      ├─ ActivityBar.java
      ├─ ExplorerPanel.java
      ├─ OutputPanel.java
      ├─ HeaderBar.java
      ├─ StatusBar.java
      ├─ DataPreparationPanel.java
      ├─ CECManagementPanel.java
      └─ CalculationPanel.java
```

---

## Input Data Files

Located in `data/CEWorkbench/inputs/` (runtime workspace) and mirrored in `app/src/main/resources/` (bundled defaults):

| File | Description |
|------|-------------|
| `clus/BCC_A2-T.txt` | BCC A2 (disordered) cluster coordinates |
| `clus/BCC_B2-T.txt` | BCC B2 (ordered) cluster coordinates |
| `clus/FCC_A1-TO.txt` | FCC A1 cluster coordinates |
| `sym/BCC_A2-SG.txt` | BCC A2 space group symmetry operations |
| `sym/BCC_B2-SG.txt` | BCC B2 space group symmetry operations |
| `sym/FCC_A1-SG.txt` | FCC A1 space group symmetry operations |
| `sym/HCP_A3-SG.txt` | HCP A3 space group symmetry operations |

---

## Workspace Layout

All persistent data is stored under `~/CEWorkbench/`:

```
~/CEWorkbench/
 ├─ inputs/
 │   ├─ clus/   ← cluster coordinate files
 │   └─ sym/    ← symmetry group files
 ├─ cluster-data/
 │   └─ <clusterId>/
 │       └─ cluster_data.json
 └─ hamiltonians/
     └─ <hamiltonianId>/
         └─ hamiltonian.json    ← edit a and b ECI values here
```

`SystemId` derives IDs deterministically:
- `clusterId`      = `<elements>_<structure>_<model>`  e.g. `Nb-Ti_BCC_A2_T`
- `hamiltonianId`  = same format, e.g. `Nb-Ti_BCC_A2_T`

---

## Project Status

**Complete:**
- Layered architecture with clean dependency rules
- Cluster identification — Stages 1, 2, 3 (cluster, CF, C-matrix)
- CVM thermodynamic engine
- MCS thermodynamic engine
- Hamiltonian scaffold, load, edit, save workflow
- VS Code dark GUI — Activity Bar, Explorer, Output panel, Header, Status bar
- Shared `WorkbenchContext` — system identity propagated across all panels
- CLI with positional system-identity arguments

**In Progress / Planned:**
- Additional symmetry groups and crystal structures
- Phase diagram grid scan visualization
- Export results to CSV / JSON

---

## Design Principles

1. **Physics-first** — domain layer contains pure physics; no IO, no UI
2. **Minimal abstractions** — no unnecessary factories, repositories, or adapters
3. **Clean separation** — each layer has a single responsibility
4. **Scientific clarity** — structure mirrors the actual thermodynamic workflow

---

## License

[Add your license here]

## Authors

- Project: CE Thermodynamics Workbench
