# CVM Calculation Pipeline: Detailed Data Flow & Execution Trace

This document provides a step‑by‑step trace of the Cluster Variation Method (CVM) calculation as it runs in **CEWorkbench**.  It follows the data from the initial user request through cluster identification, random‑state initialization, full‑CF construction, evaluation of cluster variables, and the Newton‑Raphson minimisation loop.

---

## 1. Entry Point

- **Class:** `ThermodynamicWorkflow`
- **Method:** `runCalculation(ClusterIdentificationRequest request)`
- **What happens:**
  1. Parses the request (temperature, composition, Hamiltonian ID).
  2. Calls `resolveHamiltonianIdForEngine` to ensure a **CVCF** version of the Hamiltonian is available.
  3. Loads `AllClusterData` via `loadThermodynamicData`.
  4. Instantiates `CVMEngine` and invokes `engine.compute(request)`.

---

## 2. Engine Core (`CVMEngine`)

| Step | Method | Responsibility |
|------|--------|----------------|
| **2.1** | `compute` | Orchestrates the three‑stage CVM run (identification → initialization → minimisation). |
| **2.2** | `resolveClusterData` | Forces a fresh run of `ClusterIdentificationWorkflow.identify` to guarantee that the cluster list, CF basis, and C‑matrix are up‑to‑date. |
| **2.3** | `validateCmatEciConsistency` | Checks that the columns of the C‑matrix match the Hamiltonian’s basis indices. |

---

## 3. Cluster Identification (`ClusterIdentificationWorkflow`)

1. **Stage 1 – Cluster Generation** – discovers all symmetry‑unique clusters for the crystal structure.
2. **Stage 2 – CF Basis Construction** – builds `cfBasisIndices` (decoration patterns) for each non‑point CF.
3. **Stage 3 – C‑matrix Construction** – creates the `CMatrix.Result` containing the coefficient matrix used later.

The result is stored in `AllClusterData`, which now holds:
- `int ncf` – number of independent (non‑point) CFs.
- `int tcf` – total CF count (non‑point + K‑1 point CFs).
- `int[][] cfBasisIndices` – basis‑index decorations for each CF column.
- The **C‑matrix** itself.

---

## 4. Random‑State Initialization (`CVMGibbsModel.getInitialGuess`)

```java
double[] uRand = ClusterVariableEvaluator.computeRandomCFs(
        moleFractions, numElements, cfBasisIndices, ncf, tcf);
```
- **Inputs**: composition (`moleFractions`), number of components (`numElements`), basis decorations, `ncf`, `tcf`.
- **Output**: `uRand` – random‑state values for all **non‑point** CFs.
- This corresponds to the Mathematica rule `uRandRules`.

---

## 5. Full CF Vector Construction (`CVMGibbsModel.buildFullCFVector`)

```java
double[] uFull = ClusterVariableEvaluator.buildFullCFVector(
        uRand, moleFractions, numElements, cfBasisIndices, ncf, tcf);
```
- Copies the random non‑point CFs.
- Computes the K‑1 point CFs (`σ^k`) from composition using `ClusterMath.buildBasis`.
- Places point CFs in the correct columns as dictated by `cfBasisIndices`.
- **Result**: `uFull` – length `tcf`, ready for C‑matrix multiplication.

---

## 6. Evaluation of Cluster Variables (`ClusterVariableEvaluator.evaluate`)

```java
double[][][] cv = ClusterVariableEvaluator.evaluate(
        uFull, cmat, lcv, tcdis, lc);
```
- Multiplies each row of the C‑matrix by `uFull` (adds a constant term if present).
- Produces `cv[t][j][v]` – the probability of each cluster variety.

---

## 7. Gibbs‑Model Calculations (`CVMGibbsModel.evaluate`)

1. **Entropy** – Shannon sum over all `cv` values.
2. **Enthalpy** – Linear combination of ECI coefficients with the CF vector.
3. **Gibbs Energy** – `G = H - T * S`.
4. **Gradient (`∇G`)** and **Hessian (`H_G`)** – required for Newton‑Raphson.

---

## 8. Newton‑Raphson Minimisation (`CVMSolver.minimize`)

```
while (norm(grad) > tolerance) {
    // Solve H_G * p = -grad
    double[] p = LinearAlgebra.solve(H_G, negate(grad));
    // Determine max step α to keep 0 ≤ ρ ≤ 1
    double α = calculateStepLimit(cv, p);
    // Update CF vector
    for (int i = 0; i < uFull.length; i++) uFull[i] += α * p[i];
    // Re‑evaluate entropy, enthalpy, gradient, Hessian
    evaluate(uFull);
    logIteration(iter, gradNorm, G, H, S);
}
```
- **Logging**: each iteration prints `iter |∇G| = … G = … H = … S = …` (see lines 66‑68 of the original trace).
- Convergence ends when the gradient norm falls below the configured tolerance.

---

## 9. Final Output

After convergence the engine prints:
- **Equilibrium Gibbs free energy (G)**
- **Enthalpy (H)**
- **Entropy (S)**
- **Equilibrium CFs** – extracted from the final `uFull` vector.
- These values are returned to `ThermodynamicWorkflow`, which forwards them to the UI / CLI.

---

## 10. Logging & Diagnostics

| Component | Log Message Example | Purpose |
|-----------|--------------------|---------|
| `ThermodynamicWorkflow` | `ThermodynamicWorkflow: loading data…` | High‑level progress. |
| `AllClusterData.printSummary` | Lists identified clusters, multiplicities, and C‑matrix dimensions. |
| `CECEvaluator` | `[CVM-EXTRACT] term v4AB → basis index 3` | Shows Hamiltonian‑to‑basis mapping. |
| `CVMEngine` | `CVM THERMODYNAMIC CALCULATION` banner. |
| `CVMSolver` | `iter 0 |∇G| = 1.23e‑01 G = -12345.6789 …` | Iteration‑level trace. |
| `ThermodynamicWorkflow` | `EXIT: G=…` | Final summary. |

---

## 11. Mermaid Diagram – High‑Level Flow

```mermaid
flowchart TD
    A[User Request] --> B[ThermodynamicWorkflow.runCalculation]
    B --> C[ClusterIdentificationWorkflow.identify]
    C --> D[AllClusterData (clusters, cfBasis, C‑matrix)]
    D --> E[ClusterVariableEvaluator.computeRandomCFs]
    E --> F[ClusterVariableEvaluator.buildFullCFVector]
    F --> G[ClusterVariableEvaluator.evaluate]
    G --> H[CVMGibbsModel.evaluate (G, H, S, ∇G, H_G)]
    H --> I[CVMSolver.minimize (Newton‑Raphson)]
    I --> J[Converged uFull]
    J --> K[ThermodynamicWorkflow returns results]
    K --> L[UI / CLI display]
```

---

*This trace reflects the current implementation after the latest repository pull (2026‑04‑08). Any future refactoring should keep the method signatures and data‑flow relationships documented here.*
