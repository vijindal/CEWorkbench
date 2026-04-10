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

**Note:** In the actual implementation, the domain (model) layer is split into physics evaluators and the calculation layer holds the algorithm drivers. See section 3.5.

Physics evaluators (model layer):

    CVMGibbsModel   — CVM Gibbs functional; evaluates G/H/S/gradients given (T, x, CFs)
    LatticeConfig   — MCS supercell occupation state
    EmbeddingData   — pre-computed site→cluster-instance map
    SiteOperatorBasis — evaluates φ_α(σ) basis functions
    CvCfEvaluator   — measures CVCF variables from a LatticeConfig

Algorithm drivers (calculation layer):

    ThermodynamicEngine  — interface
    CVMEngine / CVMSolver — Newton–Raphson loop; calls CVMGibbsModel.evaluate()
    MCSEngine / MCEngine  — Metropolis sweep loop; calls model-layer objects

Responsibilities of physics evaluators:

    hold persistent state (cluster geometry, ECI, lattice topology)
    respond to queries with thermodynamic quantities
    do not drive loops or decide convergence

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

# 4. Workflow Layer

The workflow layer **coordinates domain operations**.

    org.ce.workflow

Classes:

    CalculationService
    ClusterIdentificationWorkflow
    ThermodynamicWorkflow

Responsibilities:

    load required data
    construct engines
    run calculations
    coordinate domain objects

No physics is implemented here.

------------------------------------------------------------------------

# 5. Storage Layer

Handles **all disk IO and data persistence**.

    org.ce.storage

Classes:

    ClusterDataStore
    HamiltonianStore
    Workspace

Responsibilities:

    read JSON files
    write JSON files
    resolve file paths
    manage caches

This layer replaces repositories, registries, and adapters.

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
    ├─ model                         (physics evaluators + persistent state + disk I/O)
    │  ├─ cluster                    (cluster geometry, CF basis, C-matrix, CVCF)
    │  │  └─ cvcf
    │  ├─ cvm                        (CVMGibbsModel — Gibbs functional evaluator)
    │  ├─ mcs                        (LatticeConfig, EmbeddingData, Embedding,
    │  │                              EmbeddingGenerator, SiteOperatorBasis,
    │  │                              CvCfEvaluator, Vector3D)
    │  ├─ hamiltonian                (CECEntry, CECTerm, CECEvaluator)
    │  ├─ result                     (EquilibriumState, ThermodynamicResult)
    │  └─ storage                    (Workspace, InputLoader, HamiltonianStore, DataStore)
    │
    ├─ calculation                   (algorithm drivers — use model, produce results)
    │  ├─ engine
    │  │  ├─ cvm                     (CVMEngine, CVMSolver)
    │  │  └─ mcs                     (MCSEngine, MCSRunner, MCEngine, ExchangeStep,
    │  │                              MCSampler, LocalEnergyCalc, MCResult, MCSUpdate)
    │  └─ workflow
    │     ├─ CalculationService
    │     ├─ ClusterIdentificationWorkflow
    │     ├─ cec                     (CECManagementWorkflow)
    │     └─ thermo                  (ThermodynamicWorkflow, LineScanWorkflow, ...)
    │
    └─ ui
       ├─ gui
       └─ cli

------------------------------------------------------------------------

# 8. Calculation Execution Flow

When the user runs a thermodynamic calculation:

    ┌─ UI layer ──────────────────────────────────────────────────┐
    │  GUI / CLI                                                  │
    │   provides: system info + calculation parameters            │
    │   receives: EquilibriumState → display / store results      │
    └──────────────┬──────────────────────────▲───────────────────┘
                   ↓                          │ EquilibriumState
    ┌─ calculation layer ─────────────────────┼───────────────────┐
    │  CalculationService                     │                   │
    │   ↓                                     │                   │
    │  ThermodynamicWorkflow                  │                   │
    │   ↓                                     │                   │
    │  CVMEngine / MCSEngine  ────────────────┘                   │
    │   │ query (T, x, CFs)        ▲ G, ∇G, ∇²G / sweep result  │
    └───┼──────────────────────────┼─────────────────────────────-┘
        ↓                          │
    ┌─ model layer ────────────────────────────────────────────── ┐
    │  CVMGibbsModel / LatticeConfig / EmbeddingData etc.        │
    │   (pre-built from system info; persist for full session)    │
    └─────────────────────────────────────────────────────────────┘

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
