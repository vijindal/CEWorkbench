# Simplified Final Architecture

## CE Thermodynamics Workbench

This document defines the **final simplified architecture** for the CE
Thermodynamics Workbench.\
The goal is to keep the structure **minimal, easy to understand, and
aligned with the scientific workflow**.

The architecture mirrors the two fundamental tasks of the software:

1.  **Generate reusable scientific data** (Type‑1)
2.  **Compute thermodynamic equilibrium** (Type‑2)

------------------------------------------------------------------------

# 1. Core Design Principle

The program fundamentally performs two classes of work:

    1. Generate reusable scientific data   (Type 1)
    2. Compute thermodynamic equilibrium   (Type 2)

Therefore the architecture should revolve around these tasks rather than
complex enterprise abstractions.

------------------------------------------------------------------------

# 2. High‑Level Architecture

The system is divided into three layers with strict one-way dependencies:

    model       → physics evaluators, persistent state, disk I/O
    calculation → algorithm drivers, workflow orchestration
    ui          → GUI / CLI

Top‑level package:

    org.ce
     ├─ model
     ├─ calculation
     └─ ui

Dependency rule: ui → calculation → model. Never reverse.

------------------------------------------------------------------------

# 3. Domain Layer (Physics)

The **domain layer contains only scientific concepts and algorithms**.

    org.ce.domain

Subpackages:

    domain
     ├─ cluster
     ├─ hamiltonian
     ├─ engine
     └─ result

The domain layer **must never depend on storage, workflow, or UI**.

------------------------------------------------------------------------

## 3.1 cluster

Represents structural cluster information.

Typical classes:

    ClusterData
    Cluster
    CorrelationBasis
    CMatrix

Responsibilities:

    represent cluster topology
    provide correlation basis
    store C‑matrix

------------------------------------------------------------------------

## 3.2 hamiltonian

Energy models.

    Hamiltonian
    CEHamiltonian

Responsibilities:

    evaluate energy from correlation functions
    store ECI parameters

------------------------------------------------------------------------

## 3.3 engine

**In the actual implementation, both physics evaluators AND algorithm drivers are in the model layer.**

Physics evaluators (model layer):

    CVMGibbsModel   — CVM Gibbs functional; evaluates G/H/S/gradients given (T, x, CFs)
    LocalEnergyCalc — MCS cluster-expansion energy calculator (static utility)
    LatticeConfig   — MCS supercell occupation state
    EmbeddingData   — pre-computed site→cluster-instance map
    SiteOperatorBasis — evaluates φ_α(σ) basis functions
    CvCfEvaluator   — measures CVCF variables from a LatticeConfig

Algorithm drivers (model layer):

    CVMSolver       — Newton–Raphson loop; calls CVMGibbsModel.evaluate()
    MCEngine        — Metropolis sweep loop (equilibration + averaging); calls LocalEnergyCalc
    MCSampler       — accumulates statistics (tau_int, block averaging, jackknife Cv)
    MCSRunner       — BCC supercell geometry setup; ECI transformation
    ExchangeStep    — Metropolis accept/reject for one exchange trial

Responsibilities:

    Physics evaluators: stateless; respond to queries with thermodynamic quantities
    Algorithm drivers: own loops and convergence decisions; coordinate physics evaluators
    Both persist for the duration of a session, reused across all calculation points

------------------------------------------------------------------------

## 3.4 result

Objects representing thermodynamic results.

    EquilibriumState
    EngineMetrics

Responsibilities:

    store thermodynamic properties
    represent equilibrium results
    remain immutable

------------------------------------------------------------------------

# 4. Calculation Layer (Metadata & Dispatch)

The calculation layer **provides metadata for UI discovery and dispatches execution requests**.

    org.ce.calculation

Classes:

    CalculationRegistry       — discovery: what can be calculated and required parameters
    CalculationDescriptor     — vocabulary for Properties, Modes, and Parameters
    CalculationService        — unified facade for model construction and calculation dispatch
    ThermodynamicWorkflow     — internal driver for CVM and MCS optimizers

Responsibilities:

    Define the "Calculation Vocabulary" shared by UI and CLI.
    Provide the "Source of Truth" for valid engine/property/mode combinations.
    Manage the lifecycle of ModelSession (Model Construction).
    Route specifications to model-layer optimizers (CVMSolver, MCSRunner).

------------------------------------------------------------------------

# 5. Storage Layer (Part of Model)

Handles **all disk IO and data persistence** as part of the model layer.

    org.ce.model.storage

Classes:

    Workspace           — file path resolution, workspace layout
    InputLoader         — cluster file parsing
    HamiltonianStore    — JSON persistence for ECIs
    DataStore           — unified storage interface

Responsibilities:

    read cluster files
    read/write JSON files
    resolve file paths
    manage caches

Storage is part of the model layer because the model objects own their persistence.

------------------------------------------------------------------------

# 6. UI Layer

Contains the user interfaces.

    org.ce.ui

Subpackages:

    ui
     ├─ gui
     └─ cli

Typical GUI components:

    CalculationPanel
    ResultsPanel
    CECManagementPanel
    DataPreparationPanel

Responsibilities:

    collect user inputs
    display results
    trigger workflows

The UI interacts with the system through **CalculationService**.

------------------------------------------------------------------------

# 7. Final Package Structure

    org.ce
    ├─ model                         (physics evaluators, optimizers, persistent state, disk I/O)
    │  ├─ cluster                    (cluster geometry, CF basis, C-matrix, CVCF)
    │  │  └─ cvcf
    │  ├─ cvm
    │  │  ├─ CVMGibbsModel           (Gibbs functional evaluator)
    │  │  └─ CVMSolver               (Newton–Raphson minimizer — algorithm driver)
    │  ├─ mcs
    │  │  ├─ LatticeConfig           (supercell occupation state)
    │  │  ├─ EmbeddingData, Embedding, EmbeddingGenerator
    │  │  ├─ SiteOperatorBasis       (basis function evaluator)
    │  │  ├─ CvCfEvaluator           (correlation function measurer)
    │  │  ├─ Vector3D                (lattice position arithmetic)
    │  │  ├─ MCEngine                (sweep loop driver)
    │  │  ├─ MCSampler               (statistics accumulator)
    │  │  ├─ MCSRunner               (geometry + ECI transform)
    │  │  ├─ ExchangeStep            (Metropolis step)
    │  │  ├─ LocalEnergyCalc         (energy evaluator)
    │  │  └─ MCResult                (result DTO)
    │  ├─ hamiltonian                (CECEntry, CECTerm, CECEvaluator)
    │  ├─ result                     (ThermodynamicResult)
    │  └─ storage                    (Workspace, InputLoader, HamiltonianStore, DataStore)
    │
    ├─ calculation                   (thin dispatcher — no intermediate engine classes)
    │  ├─ engine
    │  │  └─ ProgressEvent           (UI-facing event type)
    │  └─ workflow
    │     ├─ CalculationService      (public API)
    │     ├─ CECManagementWorkflow   (Hamiltonian management)
    │     └─ thermo
    │        ├─ ThermodynamicWorkflow (inlines CVM/MCS dispatch)
    │        ├─ LineScanWorkflow, GridScanWorkflow, FiniteSizeScanWorkflow
    │        └─ ThermodynamicRequest, ThermodynamicData
    │
    └─ ui
       ├─ gui
       └─ cli

------------------------------------------------------------------------

# 8. Calculation Execution Flow

When the user runs a thermodynamic calculation:

    ┌─ UI layer ────────────────────────────────────────────────────────┐
    │  GUI / CLI                                                        │
    │   provides: ModelSpecifications + CalculationSpecifications       │
    │   receives: ThermodynamicResult → display / store results          │
    └──────────────┬───────────────────────────────▲────────────────────┘
                   ↓                               │ ThermodynamicResult
    ┌─ calculation layer (CalculationService) ─────┼────────────────────┐
    │   1. getOrBuildSession(ModelSpecifications)  │                    │
    │   2. execute(CalculationSpecifications)       │                    │
    │      ↓                                       │                    │
    │   internal workflow (CVM/MCS)  ──────────────┘                    │
    │      │ query (T, x, CFs)        ▲ Result data (G, H, S...)        │
    └──┬───┼──────────────────────────┼─────────────────────────────────┘
       ↓   │                          │
    ┌─ model layer ─────────────────────────────────────────────────────┐
    │  CVMGibbsModel / LatticeConfig / EmbeddingData etc.               │
    │   (pre-built from modelSpecs; persists for the full execution)    │
    └───────────────────────────────────────────────────────────────────┘

    Rules:
    - Model layer never communicates with UI directly.
    - All results flow: model → calculation → UI.
    - ModelSession is not a step in the flow — it is the pre-built
      container of model-layer objects passed as an argument to the
      calculation layer. Built once on session rebuild; reused
      across all scan points without re-identification.

------------------------------------------------------------------------

# 9. Cluster Identification Workflow (Type‑1)

For generating cluster data:

    GUI
     ↓
    ClusterIdentificationWorkflow
     ↓
    ClusterIdentifier (Stage 1)
    CFIdentifier (Stage 2)
    CMatrixBuilder (Stage 3 - Orthogonal Basis)
    CvCfBasisTransformer (Stage 4 - CVCF Basis)
     ↓
    ClusterDataStore.save()

------------------------------------------------------------------------

# 10. Advantages of This Architecture

This simplified architecture removes unnecessary enterprise abstractions
such as:

    repositories
    factories
    registries
    ports
    adapters

Instead it keeps:

    physics clean
    orchestration simple
    storage isolated
    UI independent

------------------------------------------------------------------------

# 11. Estimated Class Distribution

    domain      ~60 classes
    workflow    ~10 classes
    storage     ~15 classes
    ui          ~30 classes

Total \~115 classes, but clearly organized.

------------------------------------------------------------------------

# 12. Golden Rule

The most important architectural rule:

    domain must not depend on workflow
    domain must not depend on storage
    domain must not depend on ui

Everything depends on the domain --- never the reverse.

------------------------------------------------------------------------

# 13. Summary

This architecture keeps the CE Thermodynamics Workbench:

    simple
    scientifically oriented
    easy to maintain
    easy to extend

It aligns directly with the physics workflow of cluster‑expansion‑based
thermodynamic calculations.
