# CE Thermodynamics Workbench: Architecture & Contracts

This document defines the **unified final architecture** for the CE Thermodynamics Workbench, merging the simplified architectural guidelines with the established clean layer contracts. The goal is to keep the structure **minimal, easy to understand, testable, and strictly aligned with the scientific workflow**.

## Executive Summary

The program fundamentally performs two classes of work:
1. **Generate reusable scientific data** (Type-1: Cluster Identification & Hamiltonian Scaffold)
2. **Compute thermodynamic equilibrium** (Type-2: CVM Minimization & Monte Carlo Simulation)

The codebase revolves around these tasks using **clean, testable contracts at layer boundaries**:
- **Model Layer Contract:** Given `(T, x, ECI)`, return equilibrium state with raw thermodynamic values.
- **Calculation Layer Contract:** Given model outputs, compute scientific quantities (e.g. error bars, derived properties, statistical post-processing).
- **UI Layer Contract:** Present results and accept user input.

All layer boundaries are strictly unidirectional: `ui → calculation → model`. Never reverse.

---

## 1. High-Level Architecture & Dependency

The system is divided into three primary layers, plus a root-level storage inclusion:

```text
org.ce
 ├─ model       → physics evaluators, optimizers, persistent state, disk I/O
 ├─ calculation → algorithm dispatch, workflow orchestration, statistical post-processing
 └─ ui          → GUI / CLI interfaces
```

**Golden Architectural Rule:** 
`model` must never depend on `calculation` or `ui`. Everything depends on the domain model.

### Dependency Flow

```text
┌─ Model Layer (org.ce.model) ──────────────────────────┐
│  Physics Evaluators      Algorithm Drivers            │
│  • CVMGibbsModel         • CVMSolver                  │
│  • LocalEnergyCalc       • MCEngine / MCSRunner       │
│  • EmbeddingData         • MCSampler                  │
│                                                       │
│  Produces: ProgressEvent                              │
└────────────────────────┬──────────────────────────────┘
                         ↑ (Dispatches to optimizers)
┌─ Calculation Layer (org.ce.calculation) ──────────────┐
│  Orchestrators & Statistics Processors                │
│  • ThermodynamicWorkflow                              │
│  • MCSStatisticsProcessor                             │
│  • CalculationService (Unified dispatch API)          │
│                                                       │
│  Consumes: ProgressEvent for feedback                 │
└────────────────────────┬──────────────────────────────┘
                         ↑ (Accepts ThermodynamicResult)
┌─ UI Layer (org.ce.ui) ────────────────────────────────┐
│  Input & Display                                      │
│  • DynamicCalculationPanel                            │
│  • OutputPanel (consumes ProgressEvent)               │
└───────────────────────────────────────────────────────┘
```

---

## 2. Model Layer (Domain Physics & Algorithms)

The **domain layer contains only scientific concepts and algorithms**. It responds to queries with raw thermodynamic quantities.

### Packages & Responsibilities
- **`cluster`**: `ClusterData`, `CorrelationBasis`, `CMatrix`. Represents cluster topology and CVCF transformations. (The unified `ClusterCFIdentificationPipeline` handles all geometric/algebraic basis construction).
- **`hamiltonian`**: Energy models (`CECEntry`, `CECEvaluator`). Evaluates energy from CFs and stores ECI parameters.
- **`storage`**: Disk I/O handles persistence *as part of the model*. It owns reading cluster logic and saving ECIs (`Workspace`, `DataStore`).
- **`result`**: Immutable objects (`ThermodynamicResult`) representing equilibrium results.

### Physics Evaluators (Stateless)
Evaluators are pure functions; deterministic and no side effects.
- **`CVMGibbsModel.evaluate(...)`**: Evaluates Gibbs free energy, enthalpy, entropy, and thermodynamic gradients `(Gu, Guu, etc)` given correlation variables, `T`, `x`, and ECIs.
- **`MCSampler.sample(...)`**: Accumulates **raw** time series for simulation parameters (`Hmix`, `E`, `CF`). No statistics processing.

### Algorithm Drivers (Stateful Loops)
Optimizers own their convergence and loops but retain the model contract.
- **`CVMSolver.minimize(...)`**: Newton-Raphson loop. Inputs conditions (`T`, `x`, tolerance) and outputs an equilibrium state + convergence status.
- **`MCEngine.run(...)`**: Metropolis sweep loop. Runs equilibration and averaging sweeps. Outputs resulting `MCResult` containing *only raw averages*. Statistical processing (like Jackknife CV or standard error) is explicitly rejected from the model layer.

---

## 3. Calculation Layer (Metadata & Dispatch)

The calculation layer **provides metadata for UI discovery, coordinates model state, and runs higher-order statistical analysis**.

### Responsibilities
- Define "Calculation Vocabulary" via `CalculationRegistry` and `CalculationSpecifications`.
- Route calculations using `CalculationService`.
- Act as the source of truth for post-processing and statistics.

### Orchestrators
- **`ThermodynamicWorkflow.runCvm(...)`**: Executes a standard point scan. Feeds the `ModelSession` and request parameters to `CVMSolver`. Extracts scalar free energy points immediately without post-processing.
- **`ThermodynamicWorkflow.runMcs(...)`**: Coordinates the MC run by passing inputs to `MCSRunner`. It then applies the *Calculation-Layer Statistical processor*.

### Statistical Post-Processor
- **`MCSStatisticsProcessor`**: Completely separated from the simulation model. Takes raw time series supplied by `MCSampler` and performs:
  - Sokal automatic windowing for correlated integrated times (`τ_int`).
  - Strict Block Averaging with standard error of the mean (SEM).
  - Unbiased Jackknife resampling for heat capacity (`Cv`).

---

## 4. UI Layer

The UI Layer is entirely independent of system algorithms and handles only display and user dispatch logic.

### Components
- **`DynamicCalculationPanel` & `ParameterFieldFactory`**: Constructs responsive UI input components based purely on the `CalculationSpecifications` metadata provided by the layer below.
- **`DataPreparationPanel` & `CECManagementPanel`**: Specific Type-1 forms.
- **Contract**: Collect inputs, fire `CalculationService.executeScan()`, listen to `ProgressEvent` from the model (passed blindly through calculation), and paint `ThermodynamicResult` objects. **No physics or statistics logic.**

---

## 5. Type-1a: Cluster Identification Dataflow

This is the **Structural Identification Pipeline** (Stages 1-4). It transforms crystal geometry into mathematical structures (C-Matrix, CF orbits) required for computation. It is an **Always-Fresh** pipeline triggered automatically during `ModelSession` building, so that the basis strictly matches the target system and model.

**Execution Flow**:
```text
┌────────────────────────────────────────────────────────────────────────────────┐
│ UI / CalculationService / ModelSession.Builder                                 │
│                                                                                │
│ 1. Build ClusterIdentificationRequest (paths to clus/ and sym/ txt files)      │
│ 2. Call ClusterCFIdentificationPipeline.runFullWorkflow(request)               │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ ClusterCFIdentificationPipeline                                                │
│                                                                                │
│ [STAGE 1&1c]: Geometric Symmetry (ClusterIdentifier)                           │
│  - Target: Extract cluster topological symmetries.                             │
│  - Logic:  Analyze HSP/Phase orbits, calculate Kikuchi-Baker coefficients.     │
│                                                                                │
│ [STAGE 2]: Algebraic CFs (CFIdentifier)                                        │
│  - Target: Build site cluster representations.                                 │
│  - Logic:  Enumerate Correlation Function (CF) orbits.                         │
│                                                                                │
│ [STAGE 3]: Orthogonal C-Matrix (CMatrix)                                       │
│  - Target: Base mathematical foundation.                                       │
│  - Logic:  Build C-matrix using site-occupancy polynomials and configurations. │
│                                                                                │
│ [STAGE 4]: CVCF Transformation (CvCfBasis)                                     │
│  - Target: Chemical model abstraction.                                         │
│  - Logic:  Apply geometric transformation matrix T to the orthogonal system.   │
│                                                                                │
│ OUTPUT: PipelineResult (held purely in-memory in the ModelSession)             │
└────────────────────────────────────────────────────────────────────────────────┘
```

---

## 6. Type-2: CVM Minimization Dataflow

The Cluster Variation Method computes equilibrium by finding root derivatives of the Gibbs Free Energy function. 

```text
┌────────────────────────────────────────────────────────────────────────────────┐
│ Execution Landing: ThermodynamicWorkflow.runCalculation()                      │
│ - Unpacks ThermodynamicRequest (T, x)                                          │
│ - Evaluates ECI at T via CECEvaluator: eci[i] = a + b*T                        │
│ - Initializes CVMGibbsModel using ModelSession's PipelineResult & Hamiltonian. │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Optimization Loop: CVMSolver.minimize()                                        │
│                                                                                │
│ For each iteration (up to Newton-Raphson max-iter):                            │
│  1. Call CVMGibbsModel.evaluate(u, x, T, eci)                                  │
│  2. Compute physical gradient & Hessian of Gibbs free energy (G).              │
│  3. Calculate Newton step Δu = - Hessian^-1(G) * Gradient(G).                  │
│  4. Enforce domain constraints and apply line-search logic.                    │
│  5. Check convergence (||∇G|| < tolerance).                                    │
│                                                                                │
│ OUTPUT: EquilibriumResult (converged u vector, scalar G, H, S).                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Finalization: Validation and Packing                                           │
│ - ThermodynamicWorkflow maps EquilibriumResult to ThermodynamicResult.         │
│ - CalculationService dispatches result to the UI layer OutputPanel.            │
└────────────────────────────────────────────────────────────────────────────────┘
```

---

## 7. Type-2: Monte Carlo (MCS) Dataflow

Monte Carlo computes equilibrium stochastically by swapping atomic configurations across a lattice using Metropolis algorithms. 

```text
┌────────────────────────────────────────────────────────────────────────────────┐
│ Execution Landing: ThermodynamicWorkflow.runCalculation()                      │
│ - Unpacks ThermodynamicRequest (T, x, Box L, MCS blocks/sweeps).               │
│ - Evaluates ECI at T via CECEvaluator.                                         │
│ - Allocates LatticeConfig and evaluates EmbeddingData from session geometry.   │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Simulation Loop: MCSRunner & MCEngine.run()                                    │
│                                                                                │
│  1. Equilibration Phase (N sweeps):                                            │
│     - Iterates lattice exchanges via ExchangeStep.                             │
│     - Accepts/rejects based on Metropolis-Hastings and LocalEnergyCalc.        │
│                                                                                │
│  2. Averaging Phase (M sweeps):                                                │
│     - Continues site exchanges.                                                │
│     - Samples macroscopic state into MCSampler at every decorrelation interval.│
│                                                                                │
│ OUTPUT: Raw time series for Energy, Hmix, and CFs in MCSampler.                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Statistical Analysis: MCSStatisticsProcessor (Calculation Layer)               │
│ - Accepts raw time series from MCSampler.                                      │
│ - Invokes Sokal Windowing determining Integrated Autocorrelation Time (τ_int). │
│ - Applies Block Averaging computing Standard Error of the Mean (SEM).          │
│ - Executes Jackknife resampling to compute unbiased Heat Capacity (Cv) margins.│
│                                                                                │
│ OUTPUT: MCResult fortified with rigorous averages and error limits.            │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Finalization: Validation and Packing                                           │
│ - ThermodynamicWorkflow maps MCResult to ThermodynamicResult.                  │
│ - CalculationService dispatches result (paired with error bars) to the UI.     │
└────────────────────────────────────────────────────────────────────────────────┘
```

---

## 8. Testing Implications

The architecture enables strict testing isolation according to the layer semantics:

- **Model Layer Tests (Pure Unit Tests)**: No mocks needed. Evaluators (`CVMGibbsModel`, `CVMSolver`, `MCSampler`) are tested via direct mathematical assertions (`assert result.converged`, `assert sampler.getSampleCount() == N`). Run deterministically in CI.
- **Calculation Layer Tests (Integration Tests)**: Will test orchestrator lifecycles. Can mock `ModelSession` and `ProgressEvent` sinks. `MCSStatisticsProcessor` is verified purely via generated mock time-series data ensuring correct algorithmic averaging / Jackknifing.
- **UI Layer Tests (System Tests)**: Testing component visibility, proper metadata extraction, and workflow launch behaviors.
