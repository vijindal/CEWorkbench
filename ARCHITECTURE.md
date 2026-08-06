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


---

## 9. Appendix: CVM Calculation Trace (Complete Method-Level Walkthrough)

A detailed, method-by-method trace of a single CVM calculation from CLI entry to formatted output — the concrete worked example behind section 6's summary dataflow. Useful when debugging a specific step or onboarding to the CVM code path in depth.

### Phase 1: CLI Entry & Parameter Collection
**Entry Point:** `Main.main(String[] args)` → `Main.runCalcMin()` or `Main.runCalcFixed()` or **TRACE STARTS HERE for temperature scan:**

#### Phase 1 Classes & Methods
```
org.ce.ui.cli.Main
├─ main(String[] args)
│  └─ Creates ModelSpecifications(elements, structure, model, engineConfig)
│  └─ Creates CalculationSpecifications(property, mode)
│  └─ service.executeScan(modelSpecs, calcSpecs, sink, eventSink)
│  └─ ResultFormatter.table(results)
│
└─ Helper Methods:
   ├─ runCalcMin(elements, structure, model, temp, composition[])
   │  └─ service.execute(modelSpecs, calcSpecs, sink, eventSink) → single point
   │  └─ ResultFormatter.fullBlock(result)
   │
   ├─ runCalcFixed(elements, structure, model, temp, composition[], cfs[])
   │  └─ service.execute(...) with FIXED_CORRELATIONS parameter
   │
   └─ viewHamiltonian(elements, structure, model)
      └─ appCtx.getCecWorkflow().loadAndValidateCEC()
```

**Data Objects Created:**
- `ModelSpecifications`: {elements, structure, model, engineConfig}
- `CalculationSpecifications`: {property, mode, parameters (T, x, MCS params)}

---

### Phase 2: Unified Service Dispatch (Calculation Layer)
**Entry:** `CalculationService.execute()` or `CalculationService.executeScan()`

#### Phase 2a: Service - Session Management
```
org.ce.calculation.workflow.CalculationService
├─ execute(modelSpecs, calcSpecs, textSink, eventSink)
│  └─ getOrBuildSession(modelSpecs, textSink)  ← SESSION BUILD
│  └─ runAnalysis(modelSpecs, calcSpecs, textSink, eventSink)
│
├─ executeScan(modelSpecs, calcSpecs, textSink, eventSink)
│  └─ getOrBuildSession(modelSpecs, textSink)  ← SESSION BUILD
│  └─ lineScan.scan1D(session, vars, ...)
│
├─ getOrBuildSession(specs, sink) : ModelSession
│  ├─ Check cache: cachedSession
│  └─ sessionBuilder.build(systemId, engineConfig, sink)
│     ├─ ModelSession.Builder.build()  ← DETAILED IN PHASE 3
│     └─ Caches result in cachedSession
│
└─ runSinglePoint(session, request) : ThermodynamicResult
   └─ thermoWorkflow.runCalculation(session, request)
```

#### Phase 2b: Analysis Workflow Dispatch
```
CalculationService.runAnalysis()
├─ Parse varying parameters (T, X compositions)
├─ Determine scan type:
│  ├─ If 0D: runSinglePoint(session, request)
│  ├─ If 1D: lineScan.scan1D(session, ...)
│  └─ If 2D: gridScan.scan2D(session, ...)
│
└─ All scans delegate to ThermodynamicWorkflow
```

**Data Objects:**
- `ThermodynamicRequest`: {temperature, composition[], progressSink, eventSink, mcsParams}

---

### Phase 3: Model Session Construction (Critical Path)
**Entry:** `ModelSession.Builder.build(systemId, engineConfig, progressSink)`

This is the **most expensive phase** — performs cluster identification, C-matrix generation, and basis resolution.

#### Phase 3a: Cluster Identification (4-stage pipeline)

```
org.ce.model.ModelSession.Builder.build()
│
├─ [STAGE 0] Load Inputs
│  ├─ InputLoader.parseClusterFile(disorderedClusterFile) → List<Cluster>
│  ├─ InputLoader.parseSpaceGroup(disorderedSymFile) → SpaceGroup
│  ├─ InputLoader.parseClusterFile(orderedClusterFile) → List<Cluster>
│  └─ InputLoader.parseSpaceGroup(orderedSymFile) → SpaceGroup
│
├─ [STAGE 1 & 2] Cluster & CF Basis Identification
│  ├─ ClusterCFIdentificationPipeline.run(
│  │  ├─ disorderedClusters, disorderedOperations
│  │  ├─ orderedClusters, orderedOperations
│  │  ├─ transformationMatrix, translationVector
│  │  ├─ numComponents
│  │  └─ progressSink
│  │  )
│  ├─ Implements: Berechnung der Basis-Funktionen (BF) für den geordneten Zustand
│  └─ Returns: PipelineResult (holds CF basis identification)
│
├─ [STAGE 3] Orthogonal C-Matrix Generation
│  ├─ CMatrixPipeline.run(
│  │  ├─ clusterIdentResult
│  │  ├─ cfIdentResult
│  │  ├─ disorderedClusters
│  │  ├─ numComponents
│  │  └─ progressSink
│  │  )
│  ├─ Generates: Orthogonal C-matrix (inverse of correlation function basis)
│  └─ Returns: CMatrixData {cmat, lcv, wcv, cfBasisIndices}
│
├─ [STAGE 4] CVCF Basis Generation & Resolution
│  ├─ CvCfBasis.generate(
│  │  ├─ structure
│  │  ├─ pipelineResult
│  │  ├─ cmatrixData
│  │  ├─ model
│  │  └─ progressSink
│  │  )
│  ├─ Applies: Cluster Variation Method (CVM) basis transformation
│  │           (CVCF = Cluster Variation Correlation Functions)
│  └─ Returns: CvCfBasis {numNonPointCfs, totalCfs, cvcfCMatrixData, ...}
│
└─ [STAGE 5] Hamiltonian Loading & Resolution
   ├─ CECManagementWorkflow.loadAndValidateCEC(
   │  ├─ hamiltonianId (or hamiltonianId_CVCF for CVM)
   │  └─ progressSink
   │  )
   ├─ Loads ECI parameters from hamiltonian.json
   └─ Returns: CECEntry {elements, structure, model, cecTerms[]}
```

#### Phase 3b: Session Assembly

```
ModelSession constructor (package-private)
├─ systemId: SystemId
├─ clusterData: PipelineResult (from Stage 1-3)
├─ cecEntry: CECEntry (loaded Hamiltonian)
├─ resolvedHamiltonianId: String
├─ cvcfBasis: CvCfBasis (CVCF basis)
├─ engineConfig: EngineConfig
├─ numComponents: int
└─ Methods:
   ├─ label() : String
   ├─ systemId() : SystemId
   ├─ numComponents() : int
   ├─ hamiltonianId() : String
```

**Key Paths:**
- File paths derived via: `SystemId.structure`, `SystemId.model`, `SystemId.elements`
- Cluster files: `inputs/clus/{structure}-{model}.txt`
- Symmetry files: `inputs/sym/{structure}-SG.txt`
- Hamiltonian: `hamiltonians/{hamiltonianId}/hamiltonian.json`

---

### Phase 4: Thermodynamic Workflow Dispatch
**Entry:** `ThermodynamicWorkflow.runCalculation(session, request)`

```
org.ce.calculation.workflow.thermo.ThermodynamicWorkflow
│
├─ runCalculation(session, request) : ThermodynamicResult
│  ├─ Validate inputs (T > 0, composition sums to 1.0)
│  ├─ Switch on engineType:
│  │  ├─ "CVM" → runCvm(session, request)
│  │  └─ "MCS" → runMcs(session, request)
│  └─ Return: ThermodynamicResult
│
└─ CVM Path Only (detailed below)
```

---

### Phase 5: CVM Physics Evaluation & Minimization

#### Phase 5a: CVM Setup (First Call Only)

```
ThermodynamicWorkflow.runCvm(session, request)
│
├─ [First initialization] Initialize CVMGibbsModel:
│  ├─ CVMGibbsModel cvmModel = new CVMGibbsModel()
│  ├─ cvmModel.initialize(
│  │  ├─ elements
│  │  ├─ structure
│  │  ├─ model
│  │  ├─ cecEntry (from session)
│  │  └─ progressSink
│  │  )
│  │
│  └─ CVMGibbsModel.initialize() orchestrates:
│     ├─ Create ClusterIdentificationRequest from system identity
│     ├─ Load cluster/symmetry files (STAGE 0 — local copy in CVMGibbsModel)
│     ├─ Run ClusterCFIdentificationPipeline (STAGE 1-2)
│     ├─ Run CMatrixPipeline (STAGE 3)
│     ├─ Run CvCfBasis.generate() (STAGE 4)
│     ├─ Extract CVCF C-Matrix from basis (mdCvcf = basisRef.cvcfCMatrixData)
│     └─ Populate internal state:
│        ├─ tcdis, mhdis, kb, mh, lc (cluster geometry)
│        ├─ cmat, lcv, wcv (C-matrix and weights)
│        ├─ ncf (num non-point CFs)
│        ├─ tcf (total CFs = ncf + numComponents)
│        └─ basis (CvCfBasis reference)
│
└─ Print CVM header & input parameters
```

#### Phase 5b: CVM Solver Invocation

```
ThermodynamicWorkflow.runCvm() continuation
│
├─ CVMGibbsModel.getEquilibriumState(temperature, composition, tolerance, ...)
│  │
│  └─ Delegates to internal solver (cached)
│     ├─ Creates or reuses: CVMSolver solver = new CVMSolver()
│     │
│     └─ solver.minimize(
│        ├─ cvmModel (physics evaluator)
│        ├─ moleFractions (composition, length K)
│        ├─ temperature
│        ├─ tolerance (typically 1.0e-5)
│        ├─ progressSink (iteration progress)
│        └─ eventSink (progress events)
│        )
│        Returns: CVMSolver.EquilibriumResult
```

**Note:** `CVMGibbsModel` also caches the last result to avoid re-initialization on subsequent calls at same (T, x).

---

### Phase 5c: Newton-Raphson Minimization Loop

```
org.ce.model.cvm.CVMSolver.minimize()
│
├─ Initialize:
│  ├─ n = model.getNcf()  (number of non-point CVCF variables)
│  ├─ u[] = model.computeRandomCFs(moleFractions)  (initial guess)
│  ├─ trace = List<IterationSnapshot>  (debug output)
│  └─ MAX_ITER = 20
│
├─ For iter = 0 to MAX_ITER:
│  │
│  ├─ [Physics Evaluation]
│  │  ├─ ModelResult current = model.evaluate(u, moleFractions, temperature)
│  │  │
│  │  └─ CVMGibbsModel.evaluate(u, composition, T) : ModelResult
│  │     │
│  │     └─ evaluateInternal(
│  │        ├─ u[] (non-point CVCF variables, length ncf)
│  │        ├─ moleFractions (composition, length K)
│  │        ├─ temperature
│  │        ├─ cecEntry (Hamiltonian from session)
│  │        ├─ basis (CvCfBasis from session)
│  │        ├─ tcdis, mhdis, kb, mh, lc (cluster geometry)
│  │        ├─ cmat, lcv, wcv (C-matrix)
│  │        └─ ncf (num non-point CFs)
│  │        )
│  │        Returns: ModelResult {G, H, S, Gu[], Guu[][], Hu[], Su[], Suu[][], cfs[]}
│  │
│  └─ evaluateInternal() computation:
│     │
│     ├─ [1] Evaluate temperature-dependent ECI parameters
│     │  ├─ eci[] = CECEvaluator.evaluate(
│     │  │  ├─ cecEntry (loaded Hamiltonian)
│     │  │  ├─ temperature
│     │  │  ├─ basis (CVCF basis)
│     │  │  └─ "CVM" (engine type)
│     │  │  )
│     │  └─ eci[l] = a[l] - T * b[l]  (temperature-dependent effective cluster interactions)
│     │
│     ├─ [2] Build full correlation function vector
│     │  ├─ uFull[] = CMatrixPipeline.buildFullCVCFVector(
│     │  │  ├─ u[] (non-point CFs, length ncf)
│     │  │  ├─ moleFractions (point CFs, length K)
│     │  │  └─ ncf
│     │  │  )
│     │  └─ uFull = [u[0..ncf-1], moleFractions[0..K-1]]  (length tcf = ncf + K)
│     │
│     ├─ [3] Evaluate cluster variables (CV)
│     │  ├─ cv[][][] = CMatrixPipeline.evaluateCVs(
│     │  │  ├─ uFull[] (full CF vector)
│     │  │  ├─ cmat (C-matrix structure)
│     │  │  ├─ lcv (CV indices per cluster)
│     │  │  ├─ tcdis (num cluster types)
│     │  │  └─ lc (num clusters per type)
│     │  │  )
│     │  └─ cv[t][j][incv] = cluster variable value
│     │
│     ├─ [4] Compute Enthalpy (linear term)
│     │  ├─ H = Σ_l eci[l] * u[l]  (sum over non-point CFs only)
│     │  └─ Hu[l] = eci[l]  (gradient)
│     │
│     ├─ [5] Compute Entropy (non-linear, involves logarithms)
│     │  ├─ For each cluster variable: S -= R_GAS * weight * (cv * log(cv))
│     │  ├─ Su[l] = ∂S/∂u[l]  (computed via C-matrix chain rule)
│     │  └─ Suu[l1][l2] = ∂²S/∂u[l1]∂u[l2]  (Hessian, symmetric)
│     │
│     ├─ [6] Compute Gibbs free energy
│     │  ├─ G = H - T*S
│     │  ├─ Gu[l] = eci[l] - T * Su[l]  (gradient)
│     │  └─ Guu[l1][l2] = -T * Suu[l1][l2]  (Hessian)
│     │
│     └─ Return: ModelResult with all derivatives
│
│  ├─ [Convergence Check #1: Gradient]
│  │  ├─ errf = Σ|Gu[l]|  (L1 norm of gradient)
│  │  └─ If errf ≤ tolerance:
│  │     └─ CONVERGED → Return EquilibriumResult(state=CVMEquilibriumState, converged=true, ...)
│  │
│  ├─ [Newton Step]
│  │  ├─ Solve: Guu * p = -Gu  (linear system for step direction)
│  │  │  ├─ p[] = LinearAlgebra.solve(Guu, negGu)
│  │  │  └─ Uses LU decomposition via JAMA or similar
│  │  │
│  │  ├─ Calculate step size limit (physical bounds check)
│  │  │  ├─ alpha = model.calculateStepLimit(u, p, moleFractions)
│  │  │  └─ Ensures correlation functions remain in valid range
│  │  │
│  │  └─ Update trial point
│  │     └─ u_new = u + alpha * p
│  │
│  └─ [Convergence Check #2: Step Size]
│     ├─ errx = Σ|alpha * p[l]|  (L1 norm of step)
│     └─ If errx ≤ TOLX (1.0e-12):
│        └─ CONVERGED (x-convergence) → Return EquilibriumResult
│
├─ End loop
│
└─ Return EquilibriumResult:
   ├─ state: CVMEquilibriumState (u, G, H, S, etc.)
   ├─ converged: boolean
   ├─ iterations: int
   ├─ finalGradientNorm: double
   └─ trace: List<IterationSnapshot>  (debug snapshots)
```

---

### Phase 6: Equilibrium State Assembly

```
org.ce.model.cvm.CVMEquilibriumState (immutable DTO)
├─ Constructor: CVMEquilibriumState(u[], modelResult, temperature, R_GAS)
├─ Fields:
│  ├─ correlationFunctions[] : equilibrium CFs
│  ├─ gibbsEnergy : G at equilibrium
│  ├─ enthalpy : H at equilibrium
│  ├─ entropy : S at equilibrium
│  └─ (derived from modelResult)
└─ Methods:
   ├─ getGibbsEnergy() : double
   ├─ getEnthalpy() : double
   ├─ getEntropy() : double
   ├─ getCFs() : double[]
   └─ getMolarHeat() : double
```

---

### Phase 7: Result Assembly & Formatting

```
ThermodynamicWorkflow.runCvm() completion
│
├─ Extract equilibrium state
│  ├─ CVMEquilibriumState eqState = solverResult.state
│  ├─ Validate convergence: if (!solverResult.converged) throw
│  └─ Print results
│
├─ Create ThermodynamicResult (immutable DTO)
│  ├─ temperature: double
│  ├─ composition: double[]
│  ├─ gibbsEnergy: double
│  ├─ enthalpy: double
│  ├─ entropy: double
│  ├─ stdEnthalpy: double (NaN for CVM)
│  ├─ heatCapacity: double (NaN for CVM)
│  ├─ optimizedCFs: double[] (equilibrium non-point CFs)
│  ├─ avgCFs: null (MCS only)
│  └─ stdCFs: null (MCS only)
│
└─ Return ThermodynamicResult
```

---

### Phase 8: Scan Loop (if applicable)

For temperature scans (`executeScan()`), results from Phase 7 are collected in a loop:

```
org.ce.calculation.workflow.thermo.ScanWorkflows.LineScan
│
├─ scan1D(session, varyingParam, ...)
│  ├─ For each scan point (T or X composition):
│  │  ├─ Create ThermodynamicRequest with current parameters
│  │  ├─ Call ThermodynamicWorkflow.runCalculation(session, request)
│  │  └─ Collect ThermodynamicResult in list
│  │
│  └─ Return: List<ThermodynamicResult>
│
└─ ScanWorkflows.GridScan (for 2D scans)
   ├─ scan2D(session, varyingT, varyingX, ...)
   └─ Return: List<List<ThermodynamicResult>>
```

---

### Phase 9: Result Formatting & Output

```
org.ce.calculation.ResultFormatter
├─ table(List<ThermodynamicResult>) : String
│  └─ Formats as ASCII table with columns:
│     ├─ T (K)
│     ├─ x_1, x_2, ... (compositions)
│     ├─ G (J/mol)
│     ├─ H (J/mol)
│     └─ S (J/(mol·K))
│
└─ fullBlock(ThermodynamicResult) : String
   └─ Detailed output block for single point:
      ├─ System ID
      ├─ Composition
      ├─ G, H, S values
      ├─ CFs (if CVM)
      └─ Statistics (if MCS)
```

---

### Summary: Complete Class & Method Hierarchy

#### Entry Points
```
org.ce.ui.cli.Main
  main(String[])
  runCalcMin(...)
  runCalcFixed(...)
  viewHamiltonian(...)

org.ce.ui.gui.WorkbenchContext
  setSystem(elements, structure, model)
  calculateCVM(calcSpecs)
```

#### Calculation Service Layer
```
org.ce.calculation.workflow.CalculationService
  execute(modelSpecs, calcSpecs, ...)
  executeScan(modelSpecs, calcSpecs, ...)
  executeGridScan(modelSpecs, calcSpecs, ...)
  getOrBuildSession(specs, sink)
  runSinglePoint(session, request)
  runAnalysis(modelSpecs, calcSpecs, ...)

org.ce.calculation.workflow.thermo.ThermodynamicWorkflow
  runCalculation(session, request)
  runCvm(session, request)
  runMcs(session, request)
  printCvmHeader(sink)
  printInputParameters(...)
  printEquilibriumResults(...)
  validateInputs(T, composition)
  validateConvergence(result, sink)

org.ce.calculation.workflow.thermo.ScanWorkflows
  LineScan.scan1D(session, var, ...)
  GridScan.scan2D(session, varT, varX, ...)
  FiniteSizeScan.run(session, ...)
```

#### Model Session Construction
```
org.ce.model.ModelSession
  label() : String
  numComponents() : int
  Builder.build(systemId, engineConfig, sink)

org.ce.model.storage.Workspace
  hamiltonianFile(id) : String
  clustersFile(structure, model) : String

org.ce.model.cluster.InputLoader
  parseClusterFile(path) : List<Cluster>
  parseSpaceGroup(path) : SpaceGroup

org.ce.model.cluster.ClusterCFIdentificationPipeline
  run(...) : PipelineResult
  
org.ce.model.cluster.CMatrixPipeline
  run(...) : CMatrixData
  buildFullCVCFVector(u[], composition, ncf) : double[]
  evaluateCVs(uFull[], cmat, lcv, ...) : double[][][]

org.ce.model.cvm.CvCfBasis
  generate(structure, prResult, cmatData, model, sink) : CvCfBasis
  numNonPointCfs : int
  totalCfs() : int
```

#### CVM Physics Evaluator
```
org.ce.model.cvm.CVMGibbsModel
  initialize(elements, structure, model, cecEntry, sink)
  evaluate(u[], composition, T) : ModelResult
  evaluateInternal(u[], composition, T, ...) : ModelResult
  getEquilibriumState(T, composition, tol, sink, eventSink)
  computeRandomCFs(composition) : double[]
  calculateStepLimit(u, p, composition) : double
  getNcf() : int

org.ce.model.cvm.CVMSolver
  minimize(model, composition, T, tol, sink, eventSink) : EquilibriumResult
  EquilibriumResult:
    state : CVMEquilibriumState
    converged : boolean
    iterations : int
    finalGradientNorm : double
    trace : List<IterationSnapshot>

org.ce.model.cvm.CVMEquilibriumState
  getGibbsEnergy() : double
  getEnthalpy() : double
  getEntropy() : double
  getCFs() : double[]
  getMolarHeat() : double
```

#### Hamiltonian Management
```
org.ce.model.hamiltonian.CECEntry
  elements : String
  structure : String
  model : String
  cecTerms[] : CECTerm[]

org.ce.model.hamiltonian.CECEvaluator
  evaluate(cecEntry, T, basis, engineType) : double[]
  (returns eci[l] = a[l] - T * b[l] for each term)

org.ce.calculation.workflow.CECManagementWorkflow
  loadAndValidateCEC(sink, hamiltonianId) : CECEntry
  scaffoldFromClusterData(hamiltonianId, elements, structure, model)
```

#### Results & Data Transfer
```
org.ce.model.ThermodynamicResult
  temperature : double
  composition : double[]
  gibbsEnergy : double
  enthalpy : double
  entropy : double
  stdEnthalpy : double
  heatCapacity : double
  optimizedCFs : double[]
  avgCFs : double[]
  stdCFs : double[]

org.ce.calculation.CalculationSpecifications
  set(parameter, value)
  getOrDefault(parameter) : value

org.ce.calculation.ResultFormatter
  table(List<ThermodynamicResult>) : String
  fullBlock(ThermodynamicResult) : String
```

#### Progress & Events
```
org.ce.model.ProgressEvent
  EngineStart(engineType, totalSteps)
  CvmIteration(iter, G, gradNorm, H, S, u[])
  McSweep(step, totalSteps, E, acceptRate, isEquil, cfs[])

java.util.function.Consumer<String> progressSink
java.util.function.Consumer<ProgressEvent> eventSink
```

---

### Data Flow Diagram (Simplified)

```
┌──────────────────┐
│   CLI / GUI      │
│  Entry (Main)    │
└────────┬─────────┘
         │ ModelSpecifications
         │ CalculationSpecifications
         ▼
┌──────────────────────────────┐
│  CalculationService          │
│  (Dispatch & Session Mgmt)   │
└────────┬─────────────────────┘
         │ getOrBuildSession()
         ▼
┌──────────────────────────────┐
│  ModelSession.Builder.build()│
│  (Expensive: ~Stages 0-4)    │
├──────────────────────────────┤
│ ├─ InputLoader               │
│ ├─ ClusterCFIdentification   │
│ ├─ CMatrixPipeline           │
│ ├─ CvCfBasis.generate()      │
│ └─ CECManagementWorkflow     │
└────────┬─────────────────────┘
         │ ModelSession (cached)
         ▼
┌────────────────────────────────┐
│ ThermodynamicWorkflow          │
│ (Calculation Dispatcher)       │
└────────┬─────────────────────┘
         │ if engine=="CVM"
         ▼
┌────────────────────────────────┐
│ CVMGibbsModel.initialize()     │
│ CVMGibbsModel.getEquilibrium() │
└────────┬─────────────────────┘
         │
         ▼
┌────────────────────────────────┐
│ CVMSolver.minimize()           │
│ (Newton-Raphson Loop)          │
├────────────────────────────────┤
│ For each iteration:            │
│ ├─ evaluate()                  │
│ │  ├─ CECEvaluator.evaluate()  │
│ │  └─ evaluateInternal()       │
│ │     ├─ buildFullCVCFVector() │
│ │     ├─ evaluateCVs()         │
│ │     └─ compute G, H, S       │
│ ├─ Check convergence           │
│ └─ Newton step                 │
└────────┬──────────────────────┘
         │ EquilibriumResult
         ▼
┌────────────────────────────────┐
│ CVMEquilibriumState            │
│ (Immutable equilibrium state)  │
└────────┬──────────────────────┘
         │
         ▼
┌────────────────────────────────┐
│ ThermodynamicResult            │
│ (Single point result)          │
└────────┬──────────────────────┘
         │ (if scan: loop back)
         ▼
┌────────────────────────────────┐
│ ResultFormatter                │
│ (table() or fullBlock())       │
└────────┬──────────────────────┘
         │
         ▼
      Output
```

---

### Key Design Patterns

#### 1. **Session Caching**
- `ModelSession` is built once, cached, reused for all calculations in same system
- Eliminates redundant cluster identification & I/O

#### 2. **Physics Evaluator + Optimizer Separation**
- `CVMGibbsModel` = physics (evaluates G, H, S, gradients)
- `CVMSolver` = algorithm (Newton-Raphson loop)

#### 3. **Temperature-Dependent ECI**
- ECI parameters computed on-the-fly: `eci[l] = a[l] - T * b[l]`
- Enables efficient temperature scans without re-optimization

#### 4. **Correlation Function Representation**
- Full vector: `uFull[] = [non-point CFs, point CFs (composition)]`
- Minimization over non-point CFs only; point CFs fixed by composition constraint

#### 5. **Entropy via Cluster Variables**
- Entropy computed from cluster variables `cv[t][j][incv]`
- C-matrix chain rule connects CF gradients to CV gradients

#### 6. **CVCF Basis Transformation**
- CVM always uses CVCF basis (not orthogonal)
- Basis pre-computed in `ModelSession`, passed to evaluator

---

### Critical Performance Paths

#### Expensive (One-time per system)
1. **Cluster Identification** (~1-5s for BCC)
   - Geometry matching, symmetry filtering, CF basis construction

2. **C-Matrix Generation** (~0.5-2s)
   - Orthogonal basis generation for clusters

3. **CVCF Basis Resolution** (~0.1-0.5s)
   - Basis lookup/transformation from registry

#### Cheap (Per calculation point)
1. **Newton-Raphson Minimization** (10-20 iterations, ~100ms-1s)
   - Typical convergence: 5-10 iterations for well-behaved systems

2. **Physics Evaluation** (per iteration, ~10-50ms)
   - ECI lookup (linear, O(1))
   - CV evaluation (C-matrix × CF vector, O(n²))
   - Entropy integral (O(clusters))

---

### Testing & Verification Points

1. **CVMGibbsModel verification suite** (`CvmPhysicsVerifier.java`)
   - Validates entropy calculations against analytical solutions
   - Checks gradient/Hessian finite differences

2. **Convergence diagnostics** (printed to progress sink)
   - Iteration count, gradient norm, step size, property values

3. **Session integrity checks** (`ModelSession.Builder`)
   - Validates cluster data loaded correctly
   - Confirms C-matrix dimensions match basis

4. **Hamiltonian validation** (`CECManagementWorkflow.loadAndValidateCEC()`)
   - Checks all ECI terms present
   - Validates term names match cluster indexing

