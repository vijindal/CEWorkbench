# TYPE-1a: Cluster & Correlation Function Identification - Complete Dataflow

## Overview

TYPE-1a is the **Structural Identification Pipeline** (Stages 1-4) of the CVM calculation.
It transforms crystal geometry (atom coordinates) into the mathematical data structures (C-Matrix, CF orbits, Multiplicities) required for thermodynamic minimization.

In the current architecture, this pipeline is **Always-Fresh**, meaning it is triggered automatically at the start of every CVM calculation to ensure consistency between the Hamiltonian basis and the structural topology.

---

## Entry Points

### 1. CLI Entry (Manual Trigger)
**File:** `src/main/java/org/ce/ui/cli/Main.java`
**Method:** `main(args)` with `mode="type1a"`

### 2. Engine Entry (Automated Trigger)
**File:** `src/main/java/org/ce/domain/engine/cvm/CVMEngine.java`
**Method:** `resolveClusterData(input)`
Triggers the full pipeline on-the-fly for every CVM computation.

---

## Complete Execution Flow

```
╔════════════════════════════════════════════════════════════════════════════════╗
║                           ORCHESTRATION LAYER                                  ║
╚════════════════════════════════════════════════════════════════════════════════╝
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Main.java / CVMEngine.java                                                     │
│                                                                                │
│ 1. Build ClusterIdentificationRequest                                          │
│    - disorderedClusterFile: "clus/BCC_A2-T.txt" (Maximal clusters/HSP)         │
│    - orderedClusterFile:    "clus/BCC_B2-T.txt" (Maximal clusters/Phase)       │
│    - symmetryGroup:         "BCC_A2-SG"         (Symmetry operations)          │
│                                                                                │
│ 2. Call ClusterIdentificationWorkflow.identify(request)                        │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
╔════════════════════════════════════════════════════════════════════════════════╗
║                   CLUSTER IDENTIFICATION WORKFLOW                              ║
║             (ClusterIdentificationWorkflow.identify)                           ║
╚════════════════════════════════════════════════════════════════════════════════╝
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ [STAGE 1]: Load Resources (InputLoader)                                        │
│  - Parse .txt files into List<Cluster>                                         │
│  - Parse symmetry files into SpaceGroup operations                             │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ [STAGE 1c]: Cluster Identification (ClusterIdentifier)                         │
│  - Analyze HSP and Phase symmetry orbits.                                      │
│  - OUTPUT: ClusterIdentificationResult                                         │
│    - tcdis (types), tc (total), mhdis (multiplicities), kb (entropy coeffs)    │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ [STAGE 2]: CF Identification (CFIdentifier)                                    │
│  - Enumerate Correlation Function (CF) orbits with K-component basis.          │
│  - Map phase-specific CFs back to HSP symmetry types.                          │
│  - OUTPUT: CFIdentificationResult (ncf, tcf, lcf[t][j])                        │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ [STAGE 3]: Orthogonal C-Matrix (CMatrix)                                       │
│  - Build foundational C-matrix using site-occupancy polynomials.               │
│  - OUTPUT: CMatrix.Result (Orthogonal basis)                                   │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ [STAGE 4]: CVCF Transformation (CvCfBasis & Transformer)                       │
│  - 1. Retrieve T-matrix from CvCfBasis.Registry (e.g., BCC_A2, Model T).       │
│  - 2. Perform matrix multiplication: C_cvcf = C_ortho * T_matrix.              │
│  - 3. Compute Random State CFs in CVCF basis.                                  │
│  - OUTPUT: Final AllClusterData (Ready for Thermodynamic Solver)               │
└────────────────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Stage Breakdown

### STAGE 1c: Cluster Symmetry Analysis
**Class:** `ClusterIdentifier`
- **Multiplicities ($m_h$):** Calculated based on the size of symmetry orbits in the disordered phase.
- **Kikuchi-Baker ($k_b$):** Calculated via matrix inversion of the cluster containment table ($N_{ij}$). These coefficients adjust for the overcounting of entropy in the CVM sum.

### STAGE 2: CF Orbit Identification
**Class:** `CFIdentifier`
- **Decorations:** Every site in a cluster is decorated with basis functions ($s_1, s_2, \dots$).
- **Grouping:** CFs are grouped by their HSP type (e.g., all "v4AB" variants are grouped under the same type).

### STAGE 3: C-Matrix Building (Orthogonal)
**Class:** `CMatrix.buildOrthogonal`
- **Substitution Rules:** Maps every possible atomic configuration of a cluster to a linear combination of Correlation Functions.
- **Weights ($w_{cv}$):** Counts how many atomic configurations map to the same set of CF values (ensuring symmetry-equivalent configurations are grouped).

### STAGE 4: CVCF Basis Transformation
**Class:** `CvCfBasisTransformer`
- **Thermodynamic Basis:** Transforms the orthogonal C-matrix into a chemically intuitive basis (CVCF) where CFs represent physical interactions like chemical order or composition.
- **Consistency Check:** Validates that the transformed column labels match the Hamiltonian's required indices.

---

## Data Structures (The "Pipe")

### AllClusterData (Container)
| Component | Source | Purpose |
| :--- | :--- | :--- |
| `disClusterResult` | Stage 1 | Multiplicities and Entropy coefficients ($k_b$). |
| `disCFResult` | Stage 2 | CF orbit naming and indexing. |
| `cMatrixResult` | Stage 3/4 | The transformed C-Matrix used to calculate $\rho_{tjv}$. |
| `vRandEqui` | Stage 4 | The "initial guess" CF vector for the solver. |

---

## Key Files to Study

| Priority | File | Role |
| :--- | :--- | :--- |
| 1 | `ClusterIdentificationWorkflow.java` | The master pipeline orchestrator (Stages 1-4). |
| 2 | `ClusterIdentificationRequest.java` | Parameter container and resource auto-locator. |
| 3 | `CMatrix.java` | Logic for orthogonal matrix building and polynomial substitution. |
| 4 | `CvCfBasisTransformer.java` | The final transformation step from Orthogonal to CVCF. |
| 5 | `ClusterIdentifier.java` | Geometric symmetry and Kikuchi-Baker calculation. |

---

## Summary: What TYPE-1a Produces
The pipeline produces a single **`AllClusterData`** object. This object acts as the "Static Model" for CVM. It contains no temperature or composition data—it is the pure geometric and algebraic framework onto which physical parameters (ECI, T, x) are later projected.
