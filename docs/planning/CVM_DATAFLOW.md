# CVM Free Energy Calculation - Complete Dataflow

## True Entry Point (from CLI)
**File:** [Main.java:38](app/src/main/java/org/ce/ui/cli/Main.java#L38)

---

## Full Dataflow Diagram (CLI → Minimization → Energy Evaluation)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STAGE 0: CLI ENTRY POINT                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│ Main.java:38 main(args)                                                      │
│ MODE: "type2" (thermodynamic calculation)                                    │
│ SYSTEM: Nb-Ti / BCC_A2 / T  (or A-B / BCC_B2 / T default)                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STAGE 1: WORKFLOW SETUP                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│ Main.java:151                                                                │
│   CalculationService.runLineScanTemperature(                                │
│       CLUSTER_ID, HAMILTONIAN_ID,                                           │
│       composition [0.5, 0.5],                                               │
│       tStart=300K, tEnd=2000K, tStep=100K,                                  │
│       engineType="CVM"                                                      │
│   )                                                                          │
│                                                                              │
│ OUTPUT: List<ThermodynamicResult>                                           │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                STAGE 2: TEMPERATURE SCAN (LineScanWorkflow)                 │
├─────────────────────────────────────────────────────────────────────────────┤
│ LineScanWorkflow.scanTemperature()                                          │
│                                                                              │
│ For each temperature T in [300K, 400K, ..., 2000K]:                         │
│   │                                                                          │
│   └─► ThermodynamicWorkflow.runCalculation()  [calls for each T]           │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│            STAGE 3: LOAD DATA & DISPATCH TO ENGINE                          │
├─────────────────────────────────────────────────────────────────────────────┤
│ ThermodynamicWorkflow.runCalculation()                                      │
│                                                                              │
│ 1. Load cluster data from ClusterDataStore                                  │
│    └─ AllClusterData: Stage 1-3 topology                                    │
│                                                                              │
│ 2. Load Hamiltonian (ECI) from CECManagementWorkflow                        │
│    └─ CECEntry: a[l], b[l] coefficients for eci[l] = a[l] + b[l]·T        │
│                                                                              │
│ 3. Evaluate ECI at temperature:  eci[l] = a[l] + b[l]·T                    │
│                                                                              │
│ 4. Create ThermodynamicInput bundle & dispatch to CVMEngine                │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│              STAGE 4: ENGINE SETUP (CVMEngine.compute)                      │
├─────────────────────────────────────────────────────────────────────────────┤
│ CVMEngine:97                                                                │
│                                                                              │
│ 1. Validate input (cluster data, composition, temperature)                 │
│                                                                              │
│ 2. Create CVMInput from AllClusterData:                                     │
│    ├─ Stage 1: ClusterIdentificationResult (tcdis, mhdis, kb, mh, lc)     │
│    ├─ Stage 2: CFIdentificationResult (tcf, ncf, lcf)                      │
│    └─ Stage 3: CMatrixResult (cfBasisIndices, cmat, lcv, wcv)             │
│                                                                              │
│ 3. Create CVMPhaseModel.create()                                            │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
╔═════════════════════════════════════════════════════════════════════════════╗
║          STAGE 5: MINIMIZATION (CVMPhaseModel & N-R Solver)                ║
║   THIS IS WHERE THE LOOP HAPPENS — CVMFreeEnergy.evaluate() is called      ║
║   ONCE PER ITERATION OF THE NEWTON-RAPHSON SOLVER                          ║
╠═════════════════════════════════════════════════════════════════════════════╣
║ CVMPhaseModel.create() [line 198]                                          ║
║   → Sets ECI, temperature, composition                                      ║
║   → Calls ensureMinimized() [line 212]                                     ║
║                                                                              ║
║   ┌──────────────────────────────────────────────────────────┐              ║
║   │ CVMPhaseModel.minimize() [line 440]                      │              ║
║   │                                                          │              ║
║   │  Calls: NewtonRaphsonSolverSimple.solve()               │              ║
║   │          [NewtonRaphsonSolverSimple.java:194]           │              ║
║   │                                                          │              ║
║   └──────────────────────────────────────────────────────────┘              ║
║                                    │                                        ║
║                                    ▼                                        ║
║   ┌──────────────────────────────────────────────────────────┐              ║
║   │  NewtonRaphsonSolverSimple.minimize() [line 259]        │              ║
║   │                                                          │              ║
║   │  1. Initialize:                                          │              ║
║   │     u = getURand(data)  ← random-state initial guess    │              ║
║   │     cv = updateCV(data, u)                              │              ║
║   │                                                          │              ║
║   │  2. MAIN N-R LOOP: for iter = 1 to maxIter             │              ║
║   │     ┌─────────────────────────────────────────────┐     │              ║
║   │     │ ITERATION iter                              │     │              ║
║   │     │ ════════════════════════════════════════    │     │              ║
║   │     │                                             │     │              ║
║   │     │ A) EVALUATE ENERGY & DERIVATIVES           │     │              ║
║   │     │    vals = usrfun(data, u, Gu, Guu)       │     │              ║
║   │     │           [line 284, calls CVMFreeEnergy  │     │              ║
║   │     │            .evaluate() at line 376]       │     │              ║
║   │     │                                             │     │              ║
║   │     │    INPUT: u[ncf] (current CFs)            │     │              ║
║   │     │    OUTPUT: G, H, S, Gu[ncf], Guu[ncf×ncf]│     │              ║
║   │     │                                             │     │              ║
║   │     │ B) CHECK CONVERGENCE (gradient norm)      │     │              ║
║   │     │    gradNorm = ||Gu||  [line 290-294]     │     │              ║
║   │     │    if (gradNorm < tolerance)              │     │              ║
║   │     │       → CONVERGED! Return result          │     │              ║
║   │     │                                             │     │              ║
║   │     │ C) SOLVE LINEAR SYSTEM: Guu·du = -Gu     │     │              ║
║   │     │    du = LinearAlgebra.solve(Guu, -Gu)    │     │              ║
║   │     │    [line 315]                             │     │              ║
║   │     │                                             │     │              ║
║   │     │ D) STEP LIMITING (keep CVs positive)      │     │              ║
║   │     │    stpmax = stpmx(data, u, du, cv)       │     │              ║
║   │     │    [line 322]                             │     │              ║
║   │     │                                             │     │              ║
║   │     │ E) UPDATE CFs:                             │     │              ║
║   │     │    u[i] += stpmax · du[i]                 │     │              ║
║   │     │    [line 325-327]                         │     │              ║
║   │     │                                             │     │              ║
║   │     │ F) UPDATE CVs (next iteration):           │     │              ║
║   │     │    cv = updateCV(data, u)                 │     │              ║
║   │     │    [line 330]                             │     │              ║
║   │     │                                             │     │              ║
║   │     │ G) CHECK STEP SIZE (convergence test 2)   │     │              ║
║   │     │    if (stepNorm < TOLX)                   │     │              ║
║   │     │       → Check if truly converged...       │     │              ║
║   │     │                                             │     │              ║
║   │     │ H) CONTINUE to next iteration...          │     │              ║
║   │     └─────────────────────────────────────────────┘     │              ║
║   │                                                          │              ║
║   │  3. Return CVMSolverResult with:                        │              ║
║   │     - equilibriumCFs (converged u)                      │              ║
║   │     - G, H, S (at convergence)                         │              ║
║   │     - iterations, gradientNorm                         │              ║
║   │     - iterationTrace (diagnostics for all iters)       │              ║
║   │                                                          │              ║
║   └──────────────────────────────────────────────────────────┘              ║
║                                    │                                        ║
║                                    ▼                                        ║
║   ┌──────────────────────────────────────────────────────────┐              ║
║   │ FINAL EVALUATION at equilibrium                          │              ║
║   │ [CVMPhaseModel.minimize():473]                           │              ║
║   │                                                          │              ║
║   │ equilibrium = CVMFreeEnergy.evaluate(                   │              ║
║   │     equilibriumCFs,  moleFractions, T, ECI, ... )      │              ║
║   │                                                          │              ║
║   │ (This is ONE MORE CALL after convergence for record)    │              ║
║   │                                                          │              ║
║   └──────────────────────────────────────────────────────────┘              ║
║                                                                              ║
║ Result: Cache equilibrium state in CVMPhaseModel                           ║
║         isMinimized = true                                                 ║
║                                                                              ║
╚═════════════════════════════════════════════════════════════════════════════╝
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│              STAGE 6: RETURN TO THERMODYNAMIC WORKFLOW                      │
├─────────────────────────────────────────────────────────────────────────────┤
│ CVMEngine.compute() [line 150]                                              │
│   Returns: EquilibriumState {T, composition, H, G}                          │
│                                                                              │
│ ThermodynamicWorkflow.runCalculation()                                      │
│   Returns: ThermodynamicResult                                              │
│                                                                              │
│ LineScanWorkflow.scanTemperature()                                          │
│   Collects all results for each temperature                                 │
│                                                                              │
│ CalculationService.runLineScanTemperature()                                 │
│   Returns: List<ThermodynamicResult>  [one per T]                           │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     STAGE 7: DISPLAY RESULTS                                │
├─────────────────────────────────────────────────────────────────────────────┤
│ Main.java:157-163                                                           │
│                                                                              │
│ Print table:                                                                │
│   T (K)     G (J/mol)         H (J/mol)                                     │
│   ─────────────────────────────────────────                                │
│   300.0     -1234.5678        -567.8901                                     │
│   400.0     -1345.6789        -678.9012                                     │
│   ...                                                                        │
│   2000.0    -9876.5432        -5432.1098                                    │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## DETAILED STAGE 5: Inside CVMFreeEnergy.evaluate() (Called Each N-R Iteration)

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                         CVMFreeEnergy.evaluate()                                 ║
║                    Called by usrfun() at each N-R iteration                      ║
╚══════════════════════════════════════════════════════════════════════════════════╝
                                    │
                    ┌───────────────┴───────────────┐
                    │                               │
                    ▼                               ▼
        ┌─────────────────────────┐    ┌────────────────────────────────────┐
        │ buildFullCFVector()     │    │ evaluate()                         │
        │ ClusterVariableEvaluator│    │ ClusterVariableEvaluator          │
        │ [Line 113]              │    │ [Line 167]                         │
        └─────────────────────────┘    └────────────────────────────────────┘
                    │                               │
                    ▼                               ▼
        ┌─────────────────────────┐    ┌────────────────────────────────────┐
        │ RMatrixCalculator       │    │ Matrix-Vector Multiply             │
        │ .buildBasis()           │    │ (C-matrix × uFull)                 │
        │ [Line 120]              │    │ No further calls                   │
        │                         │    │                                    │
        │ INPUT:                  │    │ INPUT:                             │
        │ - numElements (K)       │    │ - uFull[tcf]                       │
        │                         │    │ - cmat[t][j][v][tcf+1]            │
        │ OUTPUT:                 │    │ - lcv[t][j]                        │
        │ - basis[K]             │    │ - tcdis, lc[t]                    │
        │   (symmetric integers)  │    │                                    │
        │   K=2: [-1, 1]         │    │ OUTPUT:                            │
        │   K=3: [-1, 0, 1]      │    │ - cv[tcdis][lc[t]][lcv[t][j]]     │
        └─────────────────────────┘    │   (cluster variables)             │
                    │                  │                                    │
                    └──────────────────┴────────────┘                     │
                                       │                                   │
                                       ▼                                   │
                    ┌──────────────────────────────┐                      │
                    │ Build Point CFs              │                      │
                    │ (indices ncf...tcf-1)        │                      │
                    │                              │                      │
                    │ pointCF[k] = Σᵢ xᵢ·tᵢᵏ⁺¹    │                      │
                    │                              │                      │
                    │ INPUT:                       │                      │
                    │ - moleFractions[K]          │                      │
                    │ - basis[K]                  │                      │
                    │ - cfBasisIndices[col][0]    │                      │
                    │                              │                      │
                    │ OUTPUT:                      │                      │
                    │ - uFull[tcf]                │                      │
                    │   (complete CF vector)      │                      │
                    └──────────────────────────────┘                      │
                                       │                                   │
                                       └───────────────┬───────────────────┘
                                                       │
                                    ┌──────────────────┴──────────────────┐
                                    │                                     │
                                    ▼                                     ▼
                    ╔═══════════════════════════════╗      ╔═══════════════════════════════╗
                    ║   ENTHALPY CALCULATION        ║      ║   ENTROPY CALCULATION         ║
                    ║   (Inline in evaluate)        ║      ║   (Inline in evaluate)        ║
                    ╠═══════════════════════════════╣      ╠═══════════════════════════════╣
                    ║ H = Σₜ mhdis[t] · Σₗ ECI[l]  ║      ║ Loop: for each CV value       ║
                    ║        · u[l]                  ║      ║ - Check if cv > EPS (1e-6)    ║
                    ║                               ║      ║ - If yes: cv·ln(cv)           ║
                    ║ Hcu[l] = mhdis[t] · ECI[l]  ║      ║ - If no: smooth C² extension ║
                    ║ (first derivative)            ║      ║                               ║
                    ║                               ║      ║ S = -R Σₜ kb[t]·mhdis[t]·Σⱼ ║
                    ║ Hcuu = 0                      ║      ║        mh[t][j]·Σᵥ wcv·cv·  ║
                    ║ (second derivative/Hessian)   ║      ║        ln(cv)                 ║
                    ║                               ║      ║                               ║
                    ║ INPUT:                        ║      ║ Scu[l] = -R Σₜ kb[t]·ms[t]·  ║
                    ║ - u[ncf] (optimisation vars) ║      ║          Σⱼ mh[t][j]·Σᵥ      ║
                    ║ - ECI[tcf]                    ║      ║          wcv·cmat[v][l]·    ║
                    ║ - mhdis[tcdis]               ║      ║          ln(cv)              ║
                    ║ - lcf[t][j] (CF counts)      ║      ║                               ║
                    ║                               ║      ║ Scuu[l1][l2] = -R Σₜ kb[t]· ║
                    ║ OUTPUT:                       ║      ║               ms[t]·Σⱼ mh·  ║
                    ║ - Hval (scalar)              ║      ║               Σᵥ wcv·cmat[l₁]║
                    ║ - Hcu[ncf]                   ║      ║               ·cmat[l₂]/cv  ║
                    ║                               ║      ║                               ║
                    ║                               ║      ║ INPUT:                        ║
                    ║                               ║      ║ - cv[tcdis][lc[t]][lcv[t][j]]║
                    ║                               ║      ║ - cmat[t][j][v][tcf]         ║
                    ║                               ║      ║ - kb[tcdis]                  ║
                    ║                               ║      ║ - mhdis[tcdis]               ║
                    ║                               ║      ║ - mh[t][j]                   ║
                    ║                               ║      ║ - wcv[t][j][v]               ║
                    ║                               ║      ║ - temperature (K)             ║
                    ║                               ║      ║ - R_GAS = 8.3144598 J/(mol·K)║
                    ║                               ║      ║ - ncf                         ║
                    ║                               ║      ║                               ║
                    ║                               ║      ║ OUTPUT:                       ║
                    ║                               ║      ║ - Sval (scalar)              ║
                    ║                               ║      ║ - Scu[ncf]                   ║
                    ║                               ║      ║ - Scuu[ncf][ncf]             ║
                    ╚═══════════════════════════════╝      ╚═══════════════════════════════╝
                                    │                                     │
                                    └──────────────────┬──────────────────┘
                                                       │
                                                       ▼
                    ╔═══════════════════════════════════════════════════════╗
                    ║       GIBBS ENERGY COMBINATION                       ║
                    ║       (Inline in evaluate)                           ║
                    ╠═══════════════════════════════════════════════════════╣
                    ║ G = H - T·S                                           ║
                    ║                                                       ║
                    ║ Gcu[l] = Hcu[l] - T·Scu[l]                           ║
                    ║                                                       ║
                    ║ Gcuu[l1][l2] = -T·Scuu[l1][l2]                       ║
                    ║                (since Hcuu = 0)                      ║
                    ║                                                       ║
                    ║ INPUT:                                                ║
                    ║ - Hval, Hcu[ncf]                                     ║
                    ║ - Sval, Scu[ncf], Scuu[ncf][ncf]                     ║
                    ║ - temperature (K)                                     ║
                    ║                                                       ║
                    ║ OUTPUT:                                               ║
                    ║ - Gval (scalar)                                      ║
                    ║ - Gcu[ncf]  (gradient/first derivative)              ║
                    ║ - Gcuu[ncf][ncf]  (Hessian/second derivative)        ║
                    ╚═══════════════════════════════════════════════════════╝
                                                       │
                                                       ▼
                    ╔═══════════════════════════════════════════════════════╗
                    ║              FINAL RESULT                            ║
                    ║              EvalResult                              ║
                    ╠═══════════════════════════════════════════════════════╣
                    ║ OUTPUTS:                                              ║
                    ║ - G       (Gibbs energy of mixing)                   ║
                    ║ - H       (Enthalpy of mixing)                       ║
                    ║ - S       (Entropy of mixing)                        ║
                    ║ - Gcu[ncf]    (gradient: ∂G/∂u)                      ║
                    ║ - Gcuu[ncf²]  (Hessian: ∂²G/∂u²)                     ║
                    ╚═══════════════════════════════════════════════════════╝
```

---

## Key Classes and Responsibilities

| Class | File | Role |
|-------|------|------|
| **Main** | `ui/cli/Main.java` | CLI entry point; dispatches to workflows |
| **CalculationService** | `workflow/CalculationService.java` | Service facade; orchestrates line/grid/point scans |
| **LineScanWorkflow** | `workflow/thermo/LineScanWorkflow.java` | Temperature/composition scanning |
| **ThermodynamicWorkflow** | `workflow/thermo/ThermodynamicWorkflow.java` | Loads data; dispatches to engines (CVM or MCS) |
| **CVMEngine** | `domain/engine/cvm/CVMEngine.java` | Engine interface impl; validates, creates CVMPhaseModel |
| **CVMPhaseModel** | `domain/engine/cvm/CVMPhaseModel.java` | Stateful thermodynamic model; triggers minimization via N-R solver |
| **NewtonRaphsonSolverSimple** | `domain/engine/cvm/NewtonRaphsonSolverSimple.java` | **The N-R minimization loop** — calls `usrfun()` repeatedly |
| **CVMFreeEnergy** | `domain/engine/cvm/CVMFreeEnergy.java` | **Free energy evaluator** — called once per N-R iteration |
| **ClusterVariableEvaluator** | `domain/cluster/ClusterVariableEvaluator.java` | Builds full CF vector; evaluates CVs from CFs |
| **RMatrixCalculator** | `domain/cluster/RMatrixCalculator.java` | Constructs basis vectors for multi-component systems |
| **LinearAlgebra** | `domain/cluster/LinearAlgebra.java` | Gaussian elimination (solves `Guu · du = -Gu`) |

---

## Execution Summary

### Temperature Scan Scenario (Main Example)
```
Main.java:151
  runLineScanTemperature(clusterId, hamiltonianId,
                        composition=[0.5,0.5],
                        T: 300K→2000K, step=100K,
                        engineType="CVM")
  │
  └─ LineScanWorkflow.scanTemperature()
      │
      └─ For each T in [300, 400, 500, ..., 2000]:  (18 temperatures)
          │
          └─ ThermodynamicWorkflow.runCalculation(T)
              │
              ├─ Load cluster data (once per T)
              ├─ Load Hamiltonian (once per T)
              ├─ Evaluate ECI(T) = a + b·T (once per T)
              │
              └─ CVMEngine.compute()
                  │
                  └─ CVMPhaseModel.create(cvmInput, eci, T, composition)
                      │
                      └─ ensure
Minimized()
                          │
                          └─ NewtonRaphsonSolverSimple.solve()
                              │
                              ├─ getURand()        ← initial guess (random state)
                              │
                              └─ N-R LOOP: iter = 1 to max_iter (typically 20-80 iters)
                                  │
                                  ├─ usrfun() [line 284]
                                  │   │
                                  │   └─ CVMFreeEnergy.evaluate() ◄─ CALLED ONCE PER ITERATION
                                  │       │
                                  │       ├─ buildFullCFVector()
                                  │       ├─ evaluate() CVs
                                  │       ├─ Compute H
                                  │       ├─ Compute S (nonlinear, with smooth extension)
                                  │       └─ Compute G = H - T·S
                                  │
                                  ├─ Check convergence: ||∇G|| < tolerance
                                  │   [line 304]
                                  │
                                  ├─ Solve: ∇²G · du = -∇G
                                  │   [line 315, LinearAlgebra.solve()]
                                  │
                                  ├─ Step limiting: compute stpmax  [line 322]
                                  │
                                  ├─ Update CFs: u += stpmax · du  [line 325]
                                  │
                                  └─ If converged → BREAK

                              └─ Final evaluation (line 473 in CVMPhaseModel):
                                  CVMFreeEnergy.evaluate()  ◄─ CALLED ONCE MORE


TOTAL CALLS TO CVMFreeEnergy.evaluate() FOR THIS SCENARIO:
  = (# temperatures) × (avg # N-R iterations)
  = 18 × 40 (typical)
  ≈ 720 evaluations for a full temperature scan
```

### Single-Point Calculation (Simpler)
```
CalculationService.runSinglePoint(clusterId, hamiltonianId,
                                  T=1000K, composition=[0.5,0.5],
                                  engineType="CVM")
  │
  └─ ThermodynamicWorkflow.runCalculation(1000K)
      │
      └─ CVMEngine.compute()
          └─ CVMPhaseModel.create()
              └─ NewtonRaphsonSolverSimple.solve()
                  │
                  ├─ N-R LOOP (typically 20-80 iterations)
                  │   └─ CVMFreeEnergy.evaluate()  ◄─ called per iteration
                  │
                  └─ Final evaluation
                      └─ CVMFreeEnergy.evaluate()  ◄─ called once more

TOTAL CALLS TO CVMFreeEnergy.evaluate():
  = avg # N-R iterations (typically 30-60)
```

---

## Detailed Call Chain

### Level 1: Entry Point
```
CVMFreeEnergy.evaluate(...)
├─ Input Parameters:
│  ├─ double[] u                    (non-point CF values, length ncf)
│  ├─ double[] moleFractions        (composition, length K)
│  ├─ int numElements               (K: number of components)
│  ├─ double temperature            (Kelvin)
│  ├─ double[] eci                  (effective cluster interactions)
│  ├─ List<Double> mhdis            (HSP cluster multiplicities)
│  ├─ double[] kb                   (Kikuchi-Baker entropy coefficients)
│  ├─ double[][] mh                 (normalized multiplicities)
│  ├─ int[] lc                      (ordered clusters per HSP type)
│  ├─ List<List<double[][]>> cmat   (C-matrix)
│  ├─ int[][] lcv                   (CV counts)
│  ├─ List<List<int[]>> wcv         (CV weights)
│  ├─ int tcdis                     (number of HSP cluster types)
│  ├─ int tcf                       (total number of CFs)
│  ├─ int ncf                       (number of non-point CFs)
│  ├─ int[][] lcf                   (CF count per type/group)
│  └─ int[][] cfBasisIndices        (basis-index decorations)
│
└─ Returns: EvalResult
```

### Level 2: Sub-calls

#### 2a. buildFullCFVector()
```
ClusterVariableEvaluator.buildFullCFVector(u, moleFractions, numElements,
                                           cfBasisIndices, ncf, tcf)
├─ Calls: RMatrixCalculator.buildBasis(numElements)
│  └─ Returns: double[] basis
│     For K=2: [-1, 1]
│     For K=3: [-1, 0, 1]
│     For K=4: [-2, -1, 1, 2]
│
├─ Computes: Point CFs from mole fractions and basis
│  ├─ pointCF[k] = Σᵢ moleFractions[i] · basis[i]^(k+1)
│  └─ Result: pointCFValues[nxcf]  where nxcf = tcf - ncf
│
└─ Returns: double[] uFull[tcf]
   ├─ uFull[0..ncf-1] = u (copied)
   └─ uFull[ncf..tcf-1] = placed point CF values
```

#### 2b. evaluate() - Cluster Variable Evaluation
```
ClusterVariableEvaluator.evaluate(uFull, cmat, lcv, tcdis, lc)
├─ Matrix-Vector Multiplication (innermost loop):
│  │
│  └─ for each HSP type t, group j, CV v:
│     cv[t][j][v] = Σₖ cmat[t][j][v][k] · uFull[k] + cmat[t][j][v][tcf]
│                   └─────────────────────────────────────┬──────────────────┘
│                                                   Linear combination
│
└─ Returns: double[][][] cv[tcdis][lc[t]][lcv[t][j]]
```

### Level 3: Inline Calculations (Terminal Points - no further calls)

#### 3a. Enthalpy Calculation [Lines 130-152]
```
Loop over HSP types and non-point CFs:
├─ cfOffset tracks position in CF list
├─ for each type t (except last/point type):
│  ├─ nCFsForType = Σⱼ lcf[t][j]
│  └─ for each CF l in type:
│     ├─ Hcu[l] = mhdis[t] · ECI[l]
│     └─ Hval += Hcu[l] · u[l]
│
└─ Output:
   ├─ Hval (double)
   └─ Hcu[ncf] (gradient)
```

#### 3b. Entropy Calculation [Lines 154-235]
```
Loop over HSP types, groups, and cluster variables:
├─ for each type t:
│  ├─ coeff_t = kb[t] · mhdis[t]
│  ├─ for each group j < lc[t]:
│  │  ├─ mh_tj = mh[t][j]
│  │  ├─ cm = cmat[t][j] (C-matrix for this group)
│  │  ├─ w = wcv[t][j] (weights)
│  │  ├─ nv = lcv[t][j] (number of CVs)
│  │  │
│  │  └─ for each CV v < nv:
│  │     ├─ cvVal = cv[t][j][v]
│  │     ├─ wv = w[v]
│  │     │
│  │     ├─ SMOOTH EXTENSION for small CVs:
│  │     │  ├─ if cvVal > 1e-6:
│  │     │  │  ├─ sContrib = cvVal · ln(cvVal)
│  │     │  │  ├─ logEff = ln(cvVal)
│  │     │  │  └─ invEff = 1/cvVal
│  │     │  │
│  │     │  └─ else (numerical stability):
│  │     │     ├─ sContrib = EPS·ln(EPS) + (1+ln(EPS))·d + 0.5·d²/EPS
│  │     │     ├─ logEff = ln(EPS) + d/EPS
│  │     │     └─ invEff = 1/EPS
│  │     │        where d = cvVal - EPS
│  │     │
│  │     ├─ prefix = coeff_t · mh_tj · wv
│  │     │
│  │     ├─ ENTROPY VALUE:
│  │     │  └─ Sval -= R · prefix · sContrib
│  │     │
│  │     ├─ GRADIENT (dS/du):
│  │     │  └─ for l < ncf:
│  │     │     ├─ cml = cm[v][l]
│  │     │     └─ Scu[l] -= R · prefix · cml · logEff
│  │     │
│  │     └─ HESSIAN (d²S/du²):
│  │        └─ for l1,l2 < ncf:
│  │           ├─ cml1 = cm[v][l1]
│  │           ├─ cml2 = cm[v][l2]
│  │           ├─ val = -R · prefix · cml1 · cml2 · invEff
│  │           ├─ Scuu[l1][l2] += val
│  │           └─ Scuu[l2][l1] += val (symmetric)
│  │
│
└─ Output:
   ├─ Sval (double)
   ├─ Scu[ncf] (gradient)
   └─ Scuu[ncf][ncf] (Hessian)
```

#### 3c. Gibbs Energy Combination [Lines 237-251]
```
Combine enthalpy and entropy:
├─ Gval = Hval - temperature · Sval
├─ for each l < ncf:
│  └─ Gcu[l] = Hcu[l] - temperature · Scu[l]
├─ for each l1,l2 < ncf:
│  └─ Gcuu[l1][l2] = -temperature · Scuu[l1][l2]
│                    (Hcuu = 0, so only T·Scuu term remains)
│
└─ Output:
   ├─ Gval (double)
   ├─ Gcu[ncf]
   └─ Gcuu[ncf][ncf]
```

---

## Data Structure Summary

### Input Dimensions
| Parameter | Type | Size | Description |
|-----------|------|------|-------------|
| `u` | `double[]` | ncf | Non-point CF values (optimization variables) |
| `moleFractions` | `double[]` | K | Mole fractions (Σ = 1) |
| `eci` | `double[]` | tcf | Effective cluster interactions |
| `mhdis` | `List<Double>` | tcdis | HSP multiplicities |
| `kb` | `double[]` | tcdis | Kikuchi-Baker coefficients |
| `mh` | `double[][]` | tcdis × max(lc[t]) | Normalized multiplicities |
| `cmat` | `List<List<double[][]>>` | tcdis × lc[t] × lcv[t][j] × (tcf+1) | C-matrix |
| `wcv` | `List<List<int[]>>` | tcdis × lc[t] × lcv[t][j] | CV weights |
| `lcv` | `int[][]` | tcdis × lc[t] | CV counts |

### Intermediate Data Structures
| Variable | Type | Size | Created in | Used in |
|----------|------|------|-----------|---------|
| `basis` | `double[]` | K | buildBasis() | buildFullCFVector |
| `pointCFValues` | `double[]` | nxcf | buildFullCFVector() | buildFullCFVector |
| `uFull` | `double[]` | tcf | buildFullCFVector() | evaluate() |
| `cv` | `double[][][]` | tcdis × lc[t] × lcv[t][j] | evaluate() | entropy calc |
| `Hval` | `double` | 1 | enthalpy calc | final result |
| `Hcu` | `double[]` | ncf | enthalpy calc | Gibbs energy |
| `Sval` | `double` | 1 | entropy calc | final result |
| `Scu` | `double[]` | ncf | entropy calc | Gibbs energy |
| `Scuu` | `double[][]` | ncf × ncf | entropy calc | Gibbs energy |

### Output Structure
```java
EvalResult {
  double G          // Gibbs free energy
  double H          // Enthalpy
  double S          // Entropy
  double[] Gcu      // Gradient (length ncf)
  double[][] Gcuu   // Hessian (ncf × ncf)
}
```

---

## Key Mathematical Operations

### 1. Basis Vector Construction
```
K=2:  [-1, 1]
K=3:  [-1, 0, 1]
K=4:  [-2, -1, 1, 2]
```

### 2. Point Correlation Functions
```
pointCF[k] = Σᵢ xᵢ · tᵢ^(k+1)
where tᵢ = basis[i], xᵢ = moleFractions[i]
```

### 3. Cluster Variables (Linear Combination)
```
cv[t][j][v] = Σₖ cmat[t][j][v][k] · uFull[k] + const
```

### 4. Enthalpy (Linear in CFs)
```
H = Σₜ mhdis[t] · Σₗ ECI[l] · u[l]
Hcu[l] = mhdis[t] · ECI[l]
Hcuu = 0 (linear → no curvature)
```

### 5. Entropy (Nonlinear via CV)
```
S = -R · Σₜ kb[t]·ms[t]·Σⱼ mh[t][j]·Σᵥ wcv·cv·ln(cv)

With smooth extension for cv ≤ ε = 1e-6:
  cv·ln(cv) ≈ ε·ln(ε) + (1+ln(ε))·d + 0.5·d²/ε
  where d = cv - ε
```

### 6. Gibbs Energy
```
G = H - T·S
Gcu = Hcu - T·Scu
Gcuu = -T·Scuu
```

---

## Computational Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| buildBasis | O(K) | K = numElements |
| buildFullCFVector | O(K·nxcf) | K powers + point CF placement |
| evaluate() | O(tcdis·lc·lcv·tcf) | Matrix-vector products |
| Enthalpy | O(tcdis·ncf) | Linear, single pass |
| Entropy | O(tcdis·lc·lcv·ncf²) | Dominated by Hessian assembly |
| Gibbs | O(ncf²) | Vector and matrix combinations |
| **Total** | **O(tcdis·lcv·ncf²)** | Entropy Hessian dominates |

---

## Numerical Considerations

1. **Smooth Entropy Extension** (Line 178-200)
   - Prevents log(0) when CV values are near zero
   - Uses C² continuous extension with ε = 1e-6
   - Critical for K ≥ 3 with all-zero initial guess

2. **Basis Vector** (RMatrixCalculator)
   - Symmetric integer sequence
   - Used for polynomial basis in R-matrix calculations
   - Consistent across all composition calculations

3. **Gas Constant** (Line 47)
   - R_GAS = 8.3144598 J/(mol·K)
   - ECI must be in J/mol for dimensional consistency

4. **Point CF Ordering**
   - Not necessarily ascending power order
   - Determined by CF identification pipeline
   - cfBasisIndices specifies power for each point CF

---

## WHERE IS MINIMIZATION? (Answer Summary)

### The N-R Minimization Loop Location

**File:** [NewtonRaphsonSolverSimple.java:259](app/src/main/java/org/ce/domain/engine/cvm/NewtonRaphsonSolverSimple.java#L259)

**Method:** `minimize(CVMData data, int maxIter, double tolerance)`

**Call Chain Leading to It:**
```
Main:151 → CalculationService → LineScanWorkflow
→ ThermodynamicWorkflow → CVMEngine
→ CVMPhaseModel.create():212
→ CVMPhaseModel.ensureMinimized():413
→ CVMPhaseModel.minimize():440
→ NewtonRaphsonSolverSimple.solve():228
→ NewtonRaphsonSolverSimple.minimize():259  ◄─ THE LOOP IS HERE
```

**What Happens in minimize():**

1. **Initialize** (lines 262-267):
   - u = getURand() — compute random-state CFs as initial guess
   - cv = updateCV() — evaluate cluster variables from initial u

2. **Main Loop** (lines 280-358):
   ```
   for iter = 1 to maxIter:
     1. usrfun(data, u, Gu, Guu)  ← calls CVMFreeEnergy.evaluate()
     2. Check convergence: if ||Gu|| < tolerance → exit (converged)
     3. Solve: Guu · du = -Gu  (LinearAlgebra.solve)
     4. Limit step size: stpmax = stpmx()
     5. Update CFs: u += stpmax · du
     6. Evaluate CVs: cv = updateCV()
     7. Check if step too small → may exit (converged or stalled)
   ```

3. **Return** CVMSolverResult with equilibrium CFs and trace

**How Many Times is CVMFreeEnergy.evaluate() Called?**

- **Per Single Point:** 30-100 iterations (depending on tolerance and initial guess)
- **Temperature Scan (300-2000 K, 100 K step):**
  - 18 temperatures × ~40 iterations = **~720 evaluations**
- **Each call cost:** O(tcdis · lcv · ncf²) — dominated by entropy Hessian

---

## Quick Reference: Entry Points for Different Scenarios

### Scenario 1: Temperature Scan (Most Common)
```
Main.java:151
  → CalculationService.runLineScanTemperature()
    → [For each T]
      → ThermodynamicWorkflow.runCalculation()
        → CVMEngine.compute()
          → CVMPhaseModel.minimize()  ◄─ Minimization happens here
            → NewtonRaphsonSolverSimple.minimize()  ◄─ THE LOOP
              → CVMFreeEnergy.evaluate() [many times]
```

### Scenario 2: Single Temperature Point
```
Main.java (or API)
  → CalculationService.runSinglePoint()
    → ThermodynamicWorkflow.runCalculation()
      → CVMEngine.compute()
        → CVMPhaseModel.minimize()
          → NewtonRaphsonSolverSimple.minimize()
            → CVMFreeEnergy.evaluate() [per iteration]
```

### Scenario 3: Direct Model Usage (Python/External API)
```
CVMPhaseModel.create(input, eci, T, composition)
  → ensureMinimized()
    → NewtonRaphsonSolverSimple.minimize()
      → CVMFreeEnergy.evaluate() [per iteration]
```

---

## Files to Study (Priority Order)

| Priority | File | Key Method | Purpose |
|----------|------|-----------|---------|
| **1** | [Main.java](app/src/main/java/org/ce/ui/cli/Main.java) | `main():38` | CLI entry point — start here |
| **2** | [CVMEngine.java](app/src/main/java/org/ce/domain/engine/cvm/CVMEngine.java) | `compute():32` | Engine that sets up CVM calculation |
| **3** | [CVMPhaseModel.java](app/src/main/java/org/ce/domain/engine/cvm/CVMPhaseModel.java) | `create():198` | Factory for model; triggers minimization |
| **4** | [NewtonRaphsonSolverSimple.java](app/src/main/java/org/ce/domain/engine/cvm/NewtonRaphsonSolverSimple.java) | `minimize():259` | **THE N-R LOOP (MINIMIZATION)** |
| **5** | [CVMFreeEnergy.java](app/src/main/java/org/ce/domain/engine/cvm/CVMFreeEnergy.java) | `evaluate():106` | Free energy evaluation (called per iteration) |
| **6** | [ClusterVariableEvaluator.java](app/src/main/java/org/ce/domain/cluster/ClusterVariableEvaluator.java) | `evaluate():167` | CV from CF evaluation |
| **7** | [LinearAlgebra.java](app/src/main/java/org/ce/domain/cluster/LinearAlgebra.java) | `solve()` | Solves linear system in N-R step |
