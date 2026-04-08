# CVM Calculation Pipeline: Data Flow & Class Trace

This document describes the complete calculation pipeline for the Cluster Variation Method (CVM) in the CEWorkbench project. It traces the flow of data from the initial user request to the final thermodynamic results.

## 1. High-Level Workflow Overview

The CVM calculation is managed by the `ThermodynamicWorkflow` and executed by the `CVMEngine`. The process is divided into three main phases:
1.  **Preparation**: Loading structural data (HSP cluster types, multiplicities) and the Hamiltonian (ECI).
2.  **Model Setup**: Creating the `CVMInput` topology and evaluating the Hamiltonian at the target temperature.
3.  **Minimization**: Running the Newton-Raphson solver to find the equilibrium state by minimizing the Gibbs free energy.

---

## 2. Detailed Data Flow & Class Trace

### Phase A: Preparation & Orchestration
**Entry Point:** `ThermodynamicWorkflow.runCalculation(request)`

| Class | Method | Responsibility |
| :--- | :--- | :--- |
| `ThermodynamicWorkflow` | `runCalculation` | Orchestrates the entire run: resolves Hamiltonian ID, loads data, and calls the engine. |
| `ThermodynamicWorkflow` | `resolveHamiltonianIdForEngine` | Checks for `_CVCF` versions of the Hamiltonian (required for modern CVM). |
| `ThermodynamicWorkflow` | `loadThermodynamicData` | Loads `AllClusterData` from store and validates the `CECEntry` (Hamiltonian). |
| `CVMEngine` | `compute` | The core engine loop. Triggers "Always-Fresh" identification (Stage 1-3). |
| `CVMEngine` | `resolveClusterData` | Forces a fresh run of the `ClusterIdentificationWorkflow` to ensure diagnostic consistency. |
| `ClusterIdentificationWorkflow` | `identify` | Runs the full structural identification pipeline (Stage 1: Clusters, Stage 2: CFs, Stage 3: C-Matrix). |

### Phase B: Mapping & Initialization
**Location:** `CVMEngine.compute`

| Class | Method | Responsibility |
| :--- | :--- | :--- |
| `CvCfBasis.Registry` | `get` | Retrieves the static CVCF basis transformation (T-matrix) for the given structure (e.g., BCC_A2). |
| `CECEvaluator` | `evaluate` | Evaluates ECI (a + bT) and maps them to the basis indices by matching names (e.g., `v4AB`, `v21AB`). |
| `CVMEngine` | `validateCmatEciConsistency` | Ensures the C-matrix columns in `AllClusterData` match the labels required by the `CvCfBasis`. |
| `CVMGibbsModel` | `getInitialGuess` | Computes the random-state (disordered) correlation functions to seed the solver. |
| `CvCfBasis` | `computeRandomState` | Transforms the random-state orthogonal CFs to the CVCF basis using $v = T^{-1} \cdot u$. |

### Phase C: Physical Evaluation & Minimization
**Location:** `CVMSolver.minimize`

| Class | Method | Responsibility |
| :--- | :--- | :--- |
| `CVMSolver` | `minimize` | The Newton-Raphson loop. Iteratively updates the state $u$ until the gradient $\nabla G$ is below tolerance. |
| `CVMGibbsModel` | `evaluate` | Computes the physical observables: Gibbs energy ($G$), Enthalpy ($H$), Entropy ($S$), and their derivatives ($\nabla G$, $\mathbf{H}_{G}$). |
| `ClusterVariableEvaluator` | `buildFullCVCFVector` | Combines independent optimization variables ($v$) with fixed composition ($x$) into a single vector. |
| `ClusterVariableEvaluator` | `evaluate` | Multiplies the C-matrix by the CF vector to get cluster probabilities ($\rho_{tjv} = \mathbf{C} \cdot \mathbf{v}$). |
| `CVMGibbsModel` | `evaluateInternal` | Performs the Shannon entropy summation over all cluster varieties ($\sum \rho \ln \rho$) using the $k_b$ and $m_{tj}$ coefficients. |
| `LinearAlgebra` | `solve` | Solves the linear system $\mathbf{H}_{G} \cdot \mathbf{p} = -\nabla G$ to find the Newton step direction. |
| `CVMGibbsModel` | `calculateStepLimit` | Calculates the maximum step size $\alpha$ to keep cluster probabilities in the physical range $[0, 1]$. |

---

## 3. Printing & Logging Points

The pipeline provides granular feedback to the user and logs.

### Console/UI Output (via `progressSink`)
- **Workflow Load**:
    - `ThermodynamicWorkflow`: Prints Stage 3a/b/c loading status.
    - `ThermodynamicWorkflow.logFullCecTable`: Lists all ECI terms (a, b).
- **Engine Setup**:
    - `CVMEngine`: Prints "CVM THERMODYNAMIC CALCULATION" banner.
    - `AllClusterData.printSummary`: Large block showing all identified clusters and multiplicities.
    - `CVMEngine`: Prints input parameters (Temperature, Composition).
- **Minimization Trace**:
    - `CVMSolver.minimize`: Prints one line per iteration: 
      `iter   0  |∇G| = 1.234e-01  G = -12345.6789  H = -11000.0000  S =  5.123456`
- **Final Results**:
    - `CVMEngine`: Prints "EQUILIBRIUM RESULTS" (G, H, S).
    - `ThermodynamicWorkflow`: Prints finalized J/mol results and equilibrium Correlation Functions.

### System Logs (via `java.util.logging.Logger`)
- **Mapping Details**:
    - `CECEvaluator`: Logs `[CVM-EXTRACT]` showing exactly which Hamiltonian term mapped to which basis index.
    - `CECEvaluator`: Logs `[CVM-MAPPING] Unmapped CFs` warnings if the Hamiltonian is incomplete.
- **Structural Validation**:
    - `CVMEngine`: Logs warning if falling back to orthogonal Hamiltonian instead of CVCF.
- **Trace Exit**:
    - `ThermodynamicWorkflow`: Logs `EXIT: G=...` for performance/timing monitoring.
