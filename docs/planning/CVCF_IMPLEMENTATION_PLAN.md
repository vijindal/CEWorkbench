# Plan: Implement CVCF Basis for BCC_A2

## Context

The current CVM implementation uses an orthogonal (Inden) basis where CECs are not invariant when extending from binary to ternary/quaternary systems. The Jindal-Lele (2025) CALPHAD paper introduces the **CVCF basis**: CFs chosen from Cluster Variables that are invariant across all system orders — subsystem CECs can be used directly in higher-order systems.

**Two changes are needed:**

1. **Entropy:** The cmatrix (which maps CFs to CV probabilities for entropy computation) must be transformed from the old orthogonal basis to the new CVCF basis using user-provided transformation rules.

2. **Enthalpy:** Instead of using the old orthogonal `H = Σ ECI_l · u_l`, the enthalpy is computed **directly** via equation 30 of the paper:

```
Hmc = Σ_{P<Q}      (e21PQ·v21PQ + e22PQ·v22PQ + e3PQ·v3PQ + e4PQ·v4PQ)
    + Σ_{P<Q<R}    (e3PQR1·v3PQR1 + e3PQR2·v3PQR2 + e3PQR3·v3PQR3
                  + e4PQR1·v4PQR1 + e4PQR2·v4PQR2 + e4PQR3·v4PQR3)
    + Σ_{P<Q<R<S}  (e4PQRS1·v4PQRS1 + e4PQRS2·v4PQRS2 + e4PQRS3·v4PQRS3)
```

Multiplicities are **absorbed into the ECIs** — users provide effective ECIs (e21PQ, e22PQ, etc.) that already include multiplicity.

**Key invariant:** The orthogonal basis remains an inner implementation detail only. Users see and provide CECs exclusively in the CVCF basis with names like `v21NbTi`, `v22NbTi`, `v3NbTi`, `v4NbTi`, `v3NbTiV1`, etc.

CF names are **hardcoded per (structure, numComponents)** — no generalization.

---

## Critical Files

| Role | Path |
|------|------|
| C-matrix builder (inner layer, keep) | `app/src/main/java/org/ce/domain/cluster/CMatrixBuilder.java` |
| C-matrix result | `app/src/main/java/org/ce/domain/cluster/CMatrixResult.java` |
| Free energy (entropy + enthalpy) | `app/src/main/java/org/ce/domain/cvm/CVMFreeEnergy.java` |
| CV evaluator | `app/src/main/java/org/ce/domain/cluster/ClusterVariableEvaluator.java` |
| Phase model | `app/src/main/java/org/ce/domain/cvm/CVMPhaseModel.java` |
| ECI storage | `app/src/main/java/org/ce/domain/hamiltonian/CECEntry.java` |

---

## New Files Created

```
app/src/main/java/org/ce/domain/cluster/cvcf/
  CvCfBasis.java                        ← CF name list + T matrix + const vector
  BccA2CvCfTransformations.java         ← Hardcoded T matrices (binary/ternary/quaternary)
  CvCfBasisTransformer.java             ← Applies T to old cmatrix → new cmatrix
  CvCfHamiltonianEvaluator.java         ← Computes Hmc directly via eq. 30
  CvCfCFRegistry.java                   ← Maps CF name → column index for (structure, ncomp)
```

---

## Implementation Steps Completed

### Step 1 ✅ — `CvCfBasis` (data class)

Holds all data for one (structure, numComponents) combination:

```java
public class CvCfBasis {
    public final String structurePhase;   // "BCC_A2"
    public final int numComponents;       // 2, 3, or 4
    public final List<String> cfNames;    // ordered list of all CVCF CF names
    public final int numNonPointCfs;      // number of non-point CFs (optimization variables)
    public final double[][] T;            // [numOldCfs][numNewCfs] transformation matrix
}
```

**CF name ordering for BCC_A2 (all components in alphabetical order):**

| ncomp | Non-point CFs | Point CFs | Total | Order |
|-------|---------------|-----------|-------|-------|
| 2 (A,B) | v4AB, v3AB, v22AB, v21AB | xA, xB | 6 | Tetra → Tri → II → I → Point |
| 3 (A,B,C) | v4AB, v4AC, v4BC, v4ABC1/2/3, v3AB, v3AC, v3BC, v3ABC1/2/3, v22AB, v22AC, v22BC, v21AB, v21AC, v21BC | xA, xB, xC | 21 | Tetra → Tri → II → I → Point |
| 4 (A,B,C,D) | (all tetra + tri + pair subsets in same pattern) | xA, xB, xC, xD | 55 | Tetra → Tri → II → I → Point |

Point CF count = K (all K components included explicitly).

---

### Step 2 ✅ — `BccA2CvCfTransformations` (hardcoded T matrices)

Each T matrix row corresponds to one old orthogonal CF (u[type][group][1] in user notation), each column to one new CVCF CF.

**Binary T (6 rows × 6 columns):**

```
Columns: v4AB, v3AB, v22AB, v21AB, xA, xB
```

**Ternary T (21 rows × 21 columns):**

```
Columns: v4AB, v4AC, v4BC, v4ABC1/2/3, v3AB, v3AC, v3BC, v3ABC1/2/3, v22AB, v22AC, v22BC, v21AB, v21AC, v21BC, xA, xB, xC
```

**Quaternary T (55 rows × 55 columns):**

```
Columns: v4AB, v4AC, v4AD, v4BC, v4BD, v4CD, v4ABC1/2/3, ... (all 55 CFs in tetra→tri→II→I→point order)
```

All matrices transcribed verbatim from user-provided Mathematica transformation rules.

---

### Step 3 ✅ — `CvCfBasisTransformer`

Transforms the old orthogonal cmatrix to the new CVCF cmatrix.

**Mathematical relation** per cluster group (t, j), per CV row v:

```
CV[v] = Σ_{k_old} cmat_old[v][k_old] · u_old[k_old] + cmat_old[v][tcf]
      = Σ_{k_old} cmat_old[v][k_old] · (Σ_{k_new} T[k_old][k_new]·v_new[k_new])
                + cmat_old[v][tcf]
      = Σ_{k_new} cmat_new[v][k_new] · v_new[k_new] + const_new[v]
```

Where:
- `cmat_new[v][k_new] = Σ_{k_old} cmat_old[v][k_old] · T[k_old][k_new]`
- `const_new[v] = cmat_old[v][tcf]`

```java
public CMatrixResult transform(CMatrixResult oldResult, CvCfBasis basis) {
    // For each (t, j): newCmat[t][j] = oldCmat[t][j][:, :tcf] × T
    //                  newConst[t][j] = oldCmat[t][j][:, :tcf] × constVec + oldConst[t][j]
    // Returns CMatrixResult with new cmat and cfNames=basis.cfNames
}
```

---

### Step 4 ✅ — `CvCfHamiltonianEvaluator`

Computes Hmc directly from eq. 30 (no reference to old orthogonal basis).

The new CVCF CFs (v21AB, v22AB, etc.) are read from the **same** cmatrix-evaluated CV vector used by entropy — they are rows of the new cmatrix applied to the CF optimization vector.

```java
public double computeH(Map<String, Double> cvCfValues, Map<String, Double> ecis,
                        List<String[]> componentPairs,  // [A,B], [B,C], [A,C]
                        List<String[]> componentTriples,
                        List<String[]> componentQuads) {
    double H = 0;
    // Binary terms: Σ_{P<Q} (e21PQ·v21PQ + e22PQ·v22PQ + e3PQ·v3PQ + e4PQ·v4PQ)
    for (String[] pq : componentPairs) {
        H += ecis.getOrDefault("e21"+pq[0]+pq[1], 0.0) * cvCfValues.getOrDefault(..., 0.0);
        // ... similar for e22, e3, e4
    }
    // Ternary and quaternary terms analogously
    return H;
}
```

Missing ECIs default to 0 (direct inheritance: no entry needed for zero interactions).

---

### Step 5 ✅ — `CvCfCFRegistry`

Maps CF names to column indices in CVCF basis for consistent lookup.

```java
public int indexOf(String cfName) {
    Integer idx = nameToIndex.get(cfName);
    if (idx == null) {
        throw new IllegalArgumentException("CF name not found in " + structurePhase);
    }
    return idx;
}
```

Factory methods:
- `forBccA2Binary()`
- `forBccA2Ternary()`
- `forBccA2Quaternary()`

---

## Remaining Integration Steps

### Dimensional facts (binary BCC_A2)

| Variable | Orthogonal | CVCF |
|---|---|---|
| `tcf` | 6 | 6 (same) |
| `ncf` | 4 | 4 (same) |
| `nxcf` | 2 | 2 (K components) |
| Point CFs | σ₀=1, σ₁=-xA+xB | vFull[ncf+0]=xA, vFull[ncf+1]=xB |
| `cfBasisIndices` | required | null (unused) |

T matrix is **square** (6×6 for binary). After `CvCfBasisTransformer.transform()`, the returned `CMatrixResult` has `cfBasisIndices = null`.

---

### Step 6 ✅ — `CMatrixBuilder.java`: Apply CVCF transformation at Stage-3

**File:** `app/src/main/java/org/ce/domain/cluster/CMatrixBuilder.java`

Added `CvCfBasis basis` parameter to `build()`. Builds orthogonal result then returns transformed:
```java
CMatrixResult orthResult = new CMatrixResult(cmat, lcv, wcv, cfBasisIndices);
return CvCfBasisTransformer.transform(orthResult, basis);
```

Updated callers to pass `BccA2CvCfTransformations.binaryBasis()`:
- `ClusterIdentificationWorkflow.java`
- `Main.java`
- `DataPreparationPanel.java`

**Status:** ✅ done

---

### Step 7 ✅ — `ClusterVariableEvaluator.java`: Add CVCF vector methods

**File:** `app/src/main/java/org/ce/domain/cluster/ClusterVariableEvaluator.java`

Added `buildFullCVCFVector(double[] v, double[] x, int ncf, int tcf)` — point variables = literal mole fractions.

Added `computeRandomCVCFs(double[] x, CvCfBasis basis)` — applies Tinv to orthogonal random state.

Added `Tinv` field to `CvCfBasis.java` (nullable) and `BINARY_T_INV` (6×6) to `BccA2CvCfTransformations.java`. `binaryBasis()` now passes Tinv.

**Status:** ✅ done

---

### Step 8 ✅ — `CVMFreeEnergy.java`: Drop cfBasisIndices, call buildFullCVCFVector

**File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMFreeEnergy.java`

Removed `numElements`, `lcf`, `cfBasisIndices` from `evaluate()`. Renamed `u` → `v`.
Calls `buildFullCVCFVector(v, moleFractions, ncf, tcf)`. Enthalpy simplified to `H = Σ eci[l]*v[l]`.
Binary convenience overload also updated.

**Status:** ✅ done

---

### Step 9 ✅ — `NewtonRaphsonSolverSimple.java`: Update CVMData and solver integration

**File:** `app/src/main/java/org/ce/domain/engine/cvm/NewtonRaphsonSolverSimple.java`

In `CVMData`:
- Removed `numElements`, `lcf`, `cfBasisIndices`, `cfRank`; added `CvCfBasis basis`
- Simplified constructor signature

Updates completed:
- `getURand()` → calls `ClusterVariableEvaluator.computeRandomCVCFs(data.moleFractions, data.basis)` ✅
- `updateCV()` → calls `ClusterVariableEvaluator.buildFullCVCFVector(u, data.moleFractions, data.ncf, data.tcf)` ✅
- `usrfun()` → updated `CVMFreeEnergy.evaluate()` call signature (removed `numElements`, `lcf`, `cfBasisIndices`) ✅
- `solve()` overloads → replaced params with `CvCfBasis basis` ✅

Minimization loop, step computation, convergence criteria, line search: **unchanged** ✅

**Status:** ✅ done

---

### Step 10 — `CVMPhaseModel.java`: Add CvCfBasis to CVMInput

**File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMPhaseModel.java`

- Add `CvCfBasis basis` to `CVMInput` inner class
- Remove `cfBasisIndices` field (no longer needed)
- Pass `basis` when constructing `CVMData`
- Update `buildFullCFVector()` helper → delegate to `buildFullCVCFVector()`

**Status:** ⏳ pending

---

### Step 11 — `CVMEngine.java`: Inject CvCfBasis

**File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMEngine.java`

Resolve and inject basis into `CVMInput`:
```java
CvCfBasis basis = BccA2CvCfTransformations.forBinary();  // hardcoded for now
CVMInput cvmInput = new CVMInput(..., basis);
```

TODO (future): determine from `systemId` + `numComponents` via a registry.

**Status:** ⏳ pending

---

### Step 12 — ECI Naming and Loading

New `hamiltonian.json` format for binary BCC_A2 Nb-Ti:
```json
{
  "cecTerms": [
    {"name": "e4AB", "a": 0.0, "b": 0.0},
    {"name": "e3AB", "a": 0.0, "b": 0.0},
    {"name": "e22AB", "a": -260.0, "b": 0.0},
    {"name": "e21AB", "a": -390.0, "b": 0.0}
  ]
}
```

`CvCfCFRegistry.indexOf(cfName)` maps each ECI name to its CF column index. Missing terms → ECI = 0.

**Legacy support:** Existing `"CF_N"` style entries can be converted using `ECI_new = T^T · ECI_old`.

**Status:** ⏳ pending

---

## Verification Checklist

- [ ] T matrix consistency (binary): At equimolar composition (xA=xB=0.5), evaluate old orthogonal CFs via `CMatrixResult.evaluateRandomCFs()`. Then compute the same values using T matrix applied to CVCF random values. Must match.

- [ ] CV value consistency (binary Nb-Ti): Run CV evaluation with old cmat and new cmat at same CF values (related by T). CV probabilities must be identical.

- [ ] Full CVM binary run: Nb-Ti BCC_A2 with CVCF ECIs. Compare phase boundary with existing orthogonal-basis result — physics must be identical.

- [ ] CEC inheritance (ternary): Nb-Ti-V ternary using only binary CECs (ternary ECIs = 0). CVM should run without error, producing a continuous bcc solid solution (no ternary-specific phase separation).

---

## Status

**Completed:**
- ✅ CvCfBasis data class
- ✅ BccA2CvCfTransformations with all three T matrices
- ✅ CvCfBasisTransformer implementation
- ✅ CvCfHamiltonianEvaluator implementation
- ✅ CvCfCFRegistry implementation
- ✅ CF name and ECI name lists for all system orders

**Completed integration:**
- ✅ Step 6 — CMatrixBuilder + 3 callers
- ✅ Step 7 — ClusterVariableEvaluator CVCF methods + CvCfBasis Tinv + BccA2 Tinv
- ✅ Step 8 — CVMFreeEnergy

**Pending integration:**
- ⏳ Step 9 — NewtonRaphsonSolverSimple: CVMData, getURand, updateCV, usrfun
- ⏳ Step 10 — CVMPhaseModel: CVMInput, cfBasisIndices → basis
- ⏳ Step 11 — CVMEngine: inject binaryBasis()
- ⏳ Step 12 — ECI loading from hamiltonian.json with CVCF names
- ⏳ Verification tests (matrix consistency, CV values, phase behavior)

---

## References

- **Paper:** Jindal & Lele (2025), CALPHAD — CVCF basis invariance across system orders
- **Implementation:** User-provided Mathematica transformation rules for BCC_A2 (binary, ternary, quaternary)
- **Key Invariant:** CECs are provided in CVCF basis; orthogonal basis is internal only
