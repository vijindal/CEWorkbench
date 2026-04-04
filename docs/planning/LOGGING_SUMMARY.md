# CVM Calculation Logging Summary

## Quick Start: Enable FINE-Level Logging

The logging configuration file is already created at project root:

```bash
cd d:/codes/CEWorkbench

# Option 1: Via gradle with logging enabled
./gradlew run -x test -q

# Option 2: Direct Java with logging properties
java -Djava.util.logging.config.file=logging.properties \
  -jar app/build/libs/CEWorkbench-*.jar
```

This enables `FINE` level logging which shows all entry/exit traces with parameters.

---

## Complete Dataflow with Logging

### STAGE 2: Temperature Scan Orchestration
```
LineScanWorkflow.scanTemperature()
  ├─ LOG.info — ENTER with scan parameters
  ├─ For each temperature point T:
  │  └─ LOG.info — Point N/M: T=...K
  └─ LOG.info — EXIT: completed N points
```

**Example output:**
```
[INFO] LineScanWorkflow.scanTemperature — ENTER
[INFO]   clusterId: BCC_B2_T_bin
[INFO]   hamiltonianId: A-B_BCC_B2_T
[INFO]   composition: [0.500, 0.500]
[INFO]   T range: 1000.0 to 1000.0 K, step 100.0 K
[INFO]   engineType: CVM
[INFO]   Point 1/1: T=1000.0 K
```

---

### STAGE 3: Load Data & Dispatch to Engine
```
ThermodynamicWorkflow.runCalculation()
  ├─ LOG.info — ENTER with request details (T, composition, engineType)
  │
  ├─ loadThermodynamicData()
  │  ├─ LOG.fine — STAGE 3a: Load cluster data (Stages 1-3 topology)
  │  ├─ ClusterDataStore.load(clusterId)
  │  ├─ LOG.fine — ✓ Loaded <clusterId>
  │  │
  │  ├─ LOG.fine — STAGE 3b: Load Hamiltonian (ECI coefficients)
  │  ├─ CECManagementWorkflow.loadAndValidateCEC(clusterId, hamiltonianId)
  │  ├─ LOG.fine — ✓ Loaded <hamiltonianId>
  │  ├─ LOG.fine — elements: <elements>
  │  ├─ LOG.fine — structurePhase: <phase>
  │  ├─ LOG.fine — ECI terms: N (format: a + b*T)
  │  ├─ LOG.finer —   [0] <termName>: a=<val>, b=<val>
  │  ├─ LOG.finer —   [1] <termName>: a=<val>, b=<val>
  │  ├─ LOG.finer —   [2] <termName>: a=<val>, b=<val>
  │  ├─ LOG.finer —   ... (N-3 more terms)
  │  │
  │  └─ LOG.fine — STAGE 3c: Create ThermodynamicData bundle
  │     └─ LOG.fine — systemName: <systemName>
  │
  ├─ Create ThermodynamicInput bundle
  ├─ Select engine (CVM or MCS)
  │
  ├─ CVMEngine.compute(input)  ◄─ See STAGE 4 below
  │
  └─ LOG.info — EXIT: G=<value> J/mol
```

**Example output:**
```
[INFO] ThermodynamicWorkflow.runCalculation — ENTER
[INFO]   clusterId: BCC_B2_T_bin
[INFO]   hamiltonianId: A-B_BCC_B2_T
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.5, 0.5]
[INFO]   engineType: CVM
[FINE]   STAGE 3a: Load cluster data (Stages 1-3 topology)...
[FINE]     ✓ Loaded BCC_B2_T_bin
[FINE]   STAGE 3b: Load Hamiltonian (ECI coefficients)...
[FINE]     ✓ Loaded A-B_BCC_B2_T
[FINE]     elements: A-B
[FINE]     structurePhase: BCC_B2
[FINE]     ECI terms: 4 (format: a + b*T)
[FINER]       [0] e4AB: a=0.0, b=0.0
[FINER]       [1] e3AB: a=0.0, b=0.0
[FINER]       [2] e22AB: a=-260.0, b=0.0
[FINER]       ... (1 more terms)
[FINE]   STAGE 3c: Create ThermodynamicData bundle
[FINE]     systemName: A-B_BCC_B2
```

---

### STAGE 4: CVM Engine Thermodynamic Calculation
```
CVMEngine.compute(input)
  ├─ LOG.info — === CVMEngine.compute() START ===
  ├─ LOG.info — Input parameters:
  │  ├─   systemId: <id>
  │  ├─   systemName: <name>
  │  ├─   temperature: <T> K
  │  ├─   composition: [<x1>, <x2>]
  │  └─   numComponents: <K>
  │
  ├─ Validate C-Matrix
  ├─ LOG.fine — ✓ C-Matrix found
  │
  ├─ Validate composition
  ├─ LOG.fine — ✓ Composition valid (sum=1.000000000)
  │
  ├─ Validate temperature
  ├─ LOG.fine — ✓ Temperature valid
  │
  ├─ Load C-Matrix from AllClusterData
  ├─ LOG.fine — ✓ C-Matrix loaded
  │
  ├─ Create CVMInput with topology
  ├─ LOG.fine — ✓ CVMInput created
  │
  ├─ Evaluate ECI at temperature
  │  ├─ LOG.fine — Evaluating ECI at T=<T> K...
  │  ├─ LOG.fine —   Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
  │  ├─ LOG.fine —   CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
  │  ├─ LOG.fine —   Processing N CEC terms...
  │  ├─ LOG.finer —     eci[idx] (cfName) = a + b*T = value
  │  ├─ LOG.fine —   ECI array: [values]
  │  └─ LOG.fine — ✓ ECI evaluated (4 non-point terms)
  │
  ├─ Run CVM N-R minimization
  │
  └─ CVMPhaseModel.create() ◄─ See STAGE 5 below
     ├─ CVMPhaseModel.ensureMinimized()
     │  ├─ LOG.info — CVMPhaseModel.ensureMinimized — STARTING minimization
     │  ├─ LOG.info —   systemId: <id>
     │  ├─ LOG.info —   temperature: <T> K
     │  ├─ LOG.info —   composition: [<x1>, <x2>]
     │  ├─ LOG.info —   ncf: 4, tcf: 6
     │  ├─ LOG.fine —   tolerance: 1.0e-6
     │  ├─ LOG.fine —   ECI: [<values>]
     │  │
     │  └─ minimize()
     │     ├─ NewtonRaphsonSolverSimple.solve() ◄─ See STAGE 6 below
     │     │  ├─ LOG.info — NewtonRaphsonSolverSimple.solve — ENTER
     │     │  ├─ LOG.fine —   ncf=4, tcf=6, tcdis=5, T=1000.0
     │     │  ├─ LOG.fine —   moleFractions=[0.5, 0.5]
     │     │  ├─ LOG.fine —   maxIter=400, tolerance=1.0e-6
     │     │  │
     │     │  └─ minimize(data, maxIter, tolerance)
     │     │     ├─ LOG.fine — NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1000.0, xB=0.5, tolerance=1.0e-6
     │     │     ├─ LOG.fine —   Initial CF (random state): [<values>]
     │     │     │
     │     │     ├─ For each iteration:
     │     │     │  ├─ CVMFreeEnergy.evaluate()
     │     │     │  │  ├─ LOG.finer — CVMFreeEnergy.evaluate — ENTER: T=1000.0, ncf=4, tcf=6
     │     │     │  │  └─ LOG.finer — CVMFreeEnergy.evaluate — EXIT: G=<val>, H=<val>, S=<val>
     │     │     │  │
     │     │     │  └─ LOG.finest — iter   1: G=<G> H=<H> S=<S> ||Gu||=<norm>
     │     │     │  └─ LOG.finest — iter  20: G=<G> H=<H> S=<S> ||Gu||=<norm>
     │     │     │
     │     │     ├─ On convergence:
     │     │     │  └─ LOG.fine — NewtonRaphsonSolverSimple — CONVERGED: iter=1, ||Gu||=1.87e-10
     │     │     │
     │     │     └─ Return CVMSolverResult
     │     │
     │     │  └─ LOG.info — NewtonRaphsonSolverSimple.solve — EXIT: converged=true, iterations=1, ||Gu||=1.8745e-10
     │     │
     │     └─ Compute equilibrium thermodynamic values
     │        └─ LOG.info — CVMPhaseModel.ensureMinimized — COMPLETE
     │
     └─ LOG.info — ✓ CVM minimization converged in 1 iterations (||Gu||=1.8745e-10)
        ├─   G(eq) = 2.88145994e+03 J/mol
        ├─   H(eq) = 0.00000000e+00 J/mol
        ├─   S(eq) = -2.88145994e+00 J/(mol·K)
        │
        └─ LOG.info — === CVMEngine.compute() SUCCESS ===
```

---

## Complete Logging Chain Example

Full trace for one temperature point at **T = 1000 K, x_B = 0.5** with **FINE level** enabled:

```
[INFO] LineScanWorkflow.scanTemperature — ENTER
[INFO]   clusterId: BCC_B2_T_bin
[INFO]   hamiltonianId: A-B_BCC_B2_T
[INFO]   composition: [0.500, 0.500]
[INFO]   T range: 1000.0 to 1000.0 K, step 100.0 K
[INFO]   engineType: CVM
[INFO]   Point 1/1: T=1000.0 K
[INFO] ThermodynamicWorkflow.runCalculation — ENTER
[INFO]   clusterId: BCC_B2_T_bin
[INFO]   hamiltonianId: A-B_BCC_B2_T
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.5, 0.5]
[INFO]   engineType: CVM
[FINE]   STAGE 3a: Load cluster data (Stages 1-3 topology)...
[FINE]     ✓ Loaded BCC_B2_T_bin
[FINE]   STAGE 3b: Load Hamiltonian (ECI coefficients)...
[FINE]     ✓ Loaded A-B_BCC_B2_T
[FINE]     elements: A-B
[FINE]     structurePhase: BCC_B2
[FINE]     ECI terms: 4 (format: a + b*T)
[FINER]       [0] e4AB: a=0.0, b=0.0
[FINER]       [1] e3AB: a=0.0, b=0.0
[FINER]       [2] e22AB: a=-260.0, b=0.0
[FINER]       ... (1 more terms)
[FINE]   STAGE 3c: Create ThermodynamicData bundle
[FINE]     systemName: A-B_BCC_B2
[INFO] === CVMEngine.compute() START ===
[INFO] Input parameters:
[INFO]   systemId: BCC_B2_T_bin
[INFO]   systemName: A-B_A2
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.500000, 0.500000]
[INFO]   numComponents: 2
[FINE] Validating C-Matrix...
[FINE] ✓ C-Matrix found
[FINE] Validating composition...
[FINE] ✓ Composition valid (sum=1.000000000)
[FINE] Validating temperature...
[FINE] ✓ Temperature valid
[FINE] Loading C-Matrix from AllClusterData...
[FINE] ✓ C-Matrix loaded
[FINE] Creating CVMInput with topology...
[FINE] ✓ CVMInput created
[FINE] Evaluating ECI at T=1000.0 K...
[FINE]   Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
[FINE]   CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
[FINE]   Processing 4 CEC terms...
[FINER]     eci[2] (v22AB) = -260.0 + 0.0*1000.0 = -2.6000e+02
[FINER]     eci[3] (v21AB) = -390.0 + 0.0*1000.0 = -3.9000e+02
[FINE]   ECI array: [[2]=-2.60e+02, [3]=-3.90e+02]
[FINE] ✓ ECI evaluated (4 non-point terms)
[INFO] Running CVM N-R minimization...
[INFO] CVMPhaseModel.ensureMinimized — STARTING minimization
[INFO]   systemId: BCC_B2_T_bin
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.5, 0.5]
[INFO]   ncf: 4, tcf: 6
[FINE]   tolerance: 1.0e-6
[FINE]   ECI: [[2]=-2.60e+02, [3]=-3.90e+02]
[FINE]   Calling NewtonRaphsonSolverSimple.solve()...
[INFO] NewtonRaphsonSolverSimple.solve — ENTER
[FINE]   ncf=4, tcf=6, tcdis=5, T=1000.0
[FINE]   moleFractions=[0.5, 0.5]
[FINE]   maxIter=400, tolerance=1.0e-6
[FINE] NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1000.0, xB=0.5, tolerance=1.0e-6
[FINE]   Initial CF (random state): [0.0625, 0.0, 0.25, 0.25]
[FINEST] NewtonRaphsonSolverSimple — iter   1: G=2.88145994e+03 H=0.00000000e+00 S=-2.88145994e+00 ||Gu||=1.87e-10
[FINE] NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=1, ||Gu||=1.87e-10
[FINE]   Newton-Raphson solver returned
[FINE]   Cached 1 iteration snapshots
[FINE]   Solver converged!
[FINE]   Equilibrium CFs (non-point): [0.0625, 0.0, 0.25, 0.25]
[FINE]   Evaluating CVMFreeEnergy at equilibrium CFs...
[FINER] CVMFreeEnergy.evaluate — ENTER: T=1000.0, ncf=4, tcf=6
[FINER] CVMFreeEnergy.evaluate — EXIT: G=2.88145994e+03, H=0.00000000e+00, S=-2.88145994e+00
[FINE]   Thermodynamic values: G=2.88145994e+03, H=0.00000000e+00, S=-2.88145994e+00
[INFO] CVMPhaseModel.minimize — SUCCESS: T=1000.0 K, x=[0.5, 0.5], iterations=1, ||dG||=1.87e-10, G=2881.4599 J/mol, elapsed=10 ms
[INFO] CVMPhaseModel.ensureMinimized — COMPLETE
[INFO] ✓ CVM minimization converged in 1 iterations (||Gu||=1.87e-10)
[INFO]   G(eq) = 2.88145994e+03 J/mol
[INFO]   H(eq) = 0.00000000e+00 J/mol
[INFO]   S(eq) = -2.88145994e+00 J/(mol·K)
[INFO] === CVMEngine.compute() SUCCESS ===
[INFO] ThermodynamicWorkflow.runCalculation — EXIT: G=2.8815e+03 J/mol
[INFO] LineScanWorkflow.scanTemperature — EXIT: completed 1 points
```

---

## Log Levels Reference

| Level | Shows | Usage |
|-------|-------|-------|
| `SEVERE` | Errors and exceptions | Framework errors |
| `WARNING` | Non-convergence, singular Hessian | Solver warnings |
| `INFO` | Major milestones, entry/exit | Workflow and engine milestones |
| `FINE` | Parameter summaries, step completion | Data loading, validation steps |
| `FINER` | ECI term details, intermediate values | Free energy evaluation details |
| `FINEST` | Every N-R iteration snapshot | Very detailed iteration tracing |

**Default (no config):** Shows `INFO` and above

**With logging.properties:** Shows `FINE` and above (recommended for full tracing)

---

## Configuration Files

**Location:** `d:/codes/CEWorkbench/logging.properties`

```properties
# Global level
.level = FINE

# CVM-specific loggers
org.ce.domain.engine.cvm.CVMEngine.level = FINE
org.ce.domain.engine.cvm.CVMPhaseModel.level = FINE
org.ce.domain.engine.cvm.NewtonRaphsonSolverSimple.level = FINE
org.ce.domain.engine.cvm.CVMFreeEnergy.level = FINER
org.ce.domain.cluster.ClusterVariableEvaluator.level = FINER

# Workflow loggers
org.ce.workflow.thermo.ThermodynamicWorkflow.level = INFO
org.ce.workflow.thermo.LineScanWorkflow.level = INFO

# Console handler
handlers = java.util.logging.ConsoleHandler
java.util.logging.ConsoleHandler.level = FINEST
java.util.logging.ConsoleHandler.formatter = java.util.logging.SimpleFormatter
java.util.logging.SimpleFormatter.format = [%4$-7s] %5$s%n
```

---

## Summary

- **All major methods** have entry/exit logging with input/output parameters
- **FINE level** shows data loading, validation, parameter summaries
- **FINER level** shows free energy evaluation details
- **FINEST level** shows every N-R iteration
- Configuration is **automatic** via `logging.properties` file
- Full trace captures **complete dataflow** from scan request through CVM convergence
