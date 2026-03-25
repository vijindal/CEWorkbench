# CVCF Implementation: Complete Status and Data Flow

**Date:** March 25, 2026
**Overall Status:** Phase 1 ✅ Complete | Phase 2-5 ⏳ Pending
**Progress:** 6/10 steps complete (60%)

---

## 🏗️ Architecture Decision: Dual-Pipeline with Option A

**Decision:** Maintain both orthogonal basis and CVCF basis pipelines, with CVCF as CVM-internal only.

| Component | Old Basis | CVCF Basis | Notes |
|-----------|-----------|-----------|-------|
| **hamiltonian.json** | Old basis format (CF_N indexing) | ✓ Same format | No changes needed to I/O |
| **CVMPhaseModel** | Uses cmat_old | Optionally uses cmat_new | CvCfIntegration param controls |
| **CVMFreeEnergy** | Standard computation | CVCF-aware computation | Routes based on integration param |
| **MCS Engine** | Uses old-basis ECIs | ✓ Uses old-basis ECIs | Unchanged, no basis knowledge |
| **Results** | G, H, S (from old basis) | G, H, S (basis-invariant) | Same thermodynamic properties |

**Key Points:**
1. **ECIs stay in old basis format** — hamiltonian.json unchanged, no format versioning
2. **MCS unaffected** — continues to use old-basis ECIs for Metropolis sampling
3. **CVCF is opt-in for CVM** — pass `CvCfIntegration` to CVMPhaseModel to activate
4. **Basis-invariant results** — both paths produce identical thermodynamic properties (H, S, G)
5. **Lower risk** — no changes to I/O, storage, or external engines

**When to use each:**
- **Old basis:** Default, current behavior, maximum stability
- **CVCF basis:** For systems where cluster variable representation is cleaner (future optimization)

---

## Complete CVM Calculation Data Flow

This shows a typical CVM thermodynamic calculation from user input to results, with CVCF integration points marked.

### 🔵 OLD BASIS (Current Implementation)

```
┌─────────────────────────────────────────────────────────────────┐
│                    USER / WORKFLOW LAYER                         │
│  CalculationService.run(systemId, T, composition, eccFile)      │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│              CLUSTER DATA LOADING (Stages 1-3)                   │
│  AllClusterData.load(systemId)                                   │
│  ├─ Stage 1: ClusterIdentificationResult (clusters, HSP mults)   │
│  ├─ Stage 2: CFIdentificationResult (CF count, lcf)              │
│  └─ Stage 3: CMatrixResult (cmat_old, lcv, wcv, cfBasisIndices) │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│               ECI LOADING (Old Orthogonal Basis)                 │
│  HamiltonianStore.load(system, "hamiltonian.json")               │
│  └─ Parse CECEntry[] with CF_N naming (old basis)                │
│     └─ eci[0..ncf-1] indexed by CF type (not by name)            │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│              CVMPhaseModel CREATION & MINIMIZATION                │
│  CVMPhaseModel.create(input, eci, T, moleFractions)              │
│  ├─ Store: cmat_old, eci (flat array), mhdis, kb, mh, etc      │
│  └─ ensureMinimized()                                            │
│      └─ NewtonRaphsonSolverSimple.solve(...)                     │
│          ├─ Initialize u (random or zero)                        │
│          └─ Iterate:                                             │
│              ├─ CVMFreeEnergy.evaluate(u, ..., eci, cmat_old)    │
│              │  ├─ Build uFull (add point CFs from composition)  │
│              │  ├─ Compute CV via cmat: CV = cmat·uFull + const  │
│              │  ├─ Entropy: S = -R·Σ (kb·ms·mh·w·cv·ln(cv))      │
│              │  ├─ Enthalpy: H = Σ_t (mhdis[t]·Σ_l eci[l]·u[l])  │
│              │  ├─ Gibbs: G = H - T·S                            │
│              │  ├─ Gradient: dG/du via cmat chain rule           │
│              │  └─ Hessian: d²G/du² via cmat chain rule          │
│              └─ Update u via Newton step                         │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│              RESULTS EXTRACTION & RETURN                          │
│  EquilibriumState {G, H, S, x, T, phase, ...}                    │
│  └─ WorkflowLayer returns to UI/storage                          │
└─────────────────────────────────────────────────────────────────┘
```

### 🟢 CVCF BASIS (After Full Integration)

**Design Decision (Phase 1, Option A):** CVCF basis is **CVM-internal only**. MCS and all I/O remain on old orthogonal basis.

- ECIs in hamiltonian.json stay in old basis format (no change)
- MCS engine uses old basis ECIs unchanged (no change)
- CVCF transformation happens only within CVMPhaseModel when optional
- Simpler, lower-risk, maintains backward compatibility

```
┌─────────────────────────────────────────────────────────────────┐
│                    USER / WORKFLOW LAYER                         │
│  CalculationService.run(systemId, T, composition, eccFile)      │
└──┬───────────────────────────────────────────────────────────────┘
   │
   ▼
┌──────────────────────────────────────────────────────────────────┐
│  CLUSTER DATA + ECI LOADING (Unified for both paths)             │
│  AllClusterData.load() → cmat_old, clusters, CFs                 │
│  HamiltonianStore.load() → eci[] (indexed by cluster type)       │
│                           (OLD BASIS FORMAT - unchanged)         │
└──┬──────────────────────────────────────────────────────────────┘
   │
   ├─────────────────────────┬─────────────────────────┐
   │ (OLD PATH)              │ (CVCF PATH)             │
   ▼                         ▼                         │
┌─────────────────────┐  ┌────────────────────────┐   │
│ CVMPhaseModel.run   │  │ 1. CvCfIntegration ✅  │   │
│  (no CVCF)          │  │    Transform cmat_old  │   │
│ ├─ cmat_old         │  │    → cmat_new (CVCF)   │   │
│ ├─ eci[]            │  │                        │   │
│ └─ Solve normally   │  │ 2. CVMPhaseModel.run   │   │
│                     │  │    (with CvCfIntegr.)  │   │
│                     │  │ ├─ cmat_new (CVCF)     │   │
│                     │  │ ├─ eci[] (old basis)   │   │
│                     │  │ │  [NOT transformed]   │   │
│                     │  │ ├─ Compute CV via cmN  │   │
│                     │  │ ├─ Entropy: same algo  │   │
│                     │  │ ├─ Enthalpy: same algo│   │
│                     │  │ │ (basis-invariant)    │   │
│                     │  │ └─ Minimize normally   │   │
└─────────────────────┘  └────────────────────────┘   │
   │                         │                       │
   │                         │ (Both produce same    │
   │                         │  thermodynamic       │
   │                         │  properties)         │
   └─────────────────┬───────┘                       │
                     ▼                               │
         ┌──────────────────────────┐                │
         │ EquilibriumState         │◄───────────────┘
         │ {G, H, S, x, T, phase}   │
         └──────────────────────────┘
                                │     └─ Map: "e4NbTi" → 1.5e-3     │
                                │                                   │
                                │  5. CVMPhaseModel w/ CVCF ⏳       │
                                │     ├─ Store: cmat_new, eci_map   │
                                │     └─ ensureMinimized()          │
                                │         └─ NewtonRaphsonSolver    │
                                │             ├─ Initialize u       │
                                │             └─ Iterate:           │
                                │                 ├─ CVMFreeEnergy- │
                                │                 │  CVCF.evaluate()│
                                │                 │  ⏳ [NEW PATH]   │
                                │                 │  ├─ Build uFull │
                                │                 │  ├─ Compute CV  │
                                │                 │  │  (unchanged) │
                                │                 │  ├─ Entropy:    │
                                │                 │  │  (unchanged) │
                                │                 │  ├─ Enthalpy:   │
                                │                 │  │  H = Σ e·v   │
                                │                 │  │  via         │
                                │                 │  │  CvCfHamil.. │
                                │                 │  │  Evaluator ⏳ │
                                │                 │  ├─ Gibbs:      │
                                │                 │  │  (unchanged) │
                                │                 │  ├─ Gradient ⏳ │
                                │                 │  │  dH/du via   │
                                │                 │  │  chain rule  │
                                │                 │  └─ Hessian ⏳  │
                                │                 │     (unchanged) │
                                │                 └─ Newton step    │
                                │                                   │
                                │  6. Results extraction ✅          │
                                │     EquilibriumState              │
                                │                                   │
└────────────────────────────────────────────────────────────────┘
```

---

## CVCF Integration: Step-by-Step Status

### ✅ **COMPLETED (Phase 1)**

#### **Step 1 ✅ — CvCfBasis Data Class**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/CvCfBasis.java`
- **Status:** Complete and tested
- **What it does:** Holds T matrix, CF names, non-point CF count
- **Used by:** CvCfIntegration, CvCfBasisTransformer

#### **Step 2 ✅ — BccA2CvCfTransformations (T Matrices)**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/BccA2CvCfTransformations.java`
- **Status:** Complete with all three matrices
- **Contents:**
  - Binary T (6×6): BINARY_T, BINARY_CF_NAMES, BINARY_ECI_NAMES
  - Ternary T (21×21): TERNARY_T, TERNARY_CF_NAMES, TERNARY_ECI_NAMES
  - Quaternary T (55×55): QUATERNARY_T, QUATERNARY_CF_NAMES, QUATERNARY_ECI_NAMES
- **Factory methods:** `binaryBasis()`, `ternaryBasis()`, `quaternaryBasis()`

#### **Step 3 ✅ — CvCfBasisTransformer**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/CvCfBasisTransformer.java`
- **Status:** Complete and tested
- **Formula:** `cmat_new = cmat_old · T`
- **Used by:** CvCfIntegration

#### **Step 4 ✅ — CvCfHamiltonianEvaluator**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/CvCfHamiltonianEvaluator.java`
- **Status:** Complete
- **Computes:** Hmc via equation 30 (Jindal & Lele 2025)
- **Input:** CF name maps, ECI maps, component pairs/triples/quads
- **Output:** Scalar H value (direct enthalpy)
- **Note:** Missing ECIs default to 0 (CEC inheritance)

#### **Step 5 ✅ — CvCfCFRegistry**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/CvCfCFRegistry.java`
- **Status:** Complete
- **Maps:** CF name → column index
- **Factory:** `forBccA2Binary()`, `forBccA2Ternary()`, `forBccA2Quaternary()`

#### **Step 6 ✅ — CvCfIntegration Adapter [NEW - Just Completed]**
- **File:** `app/src/main/java/org/ce/domain/cluster/cvcf/CvCfIntegration.java`
- **Status:** Complete with full test coverage
- **Key methods:**
  - `forBccA2Binary(cmat, "Nb", "Ti")` → Transform + organize
  - `flatEciToMap(eci[])` → Array to CF-name map
  - `mapEciToFlat(map)` → CF-name map to array
  - `getComponentPairs()`, `getComponentTriples()`, `getComponentQuads()`
- **Test coverage:** 11 unit tests, all passing

#### **Step 7 ✅ — Comprehensive Test Suite [NEW - Just Completed]**
- **File:** `app/src/test/java/org/ce/domain/cluster/cvcf/CvCfBasisTransformationTest.java`
- **Status:** 11/11 tests passing
- **Covers:**
  - Basis structure (dimensions, CF names)
  - T matrix validity
  - ECI array ↔ map conversions
  - CEC inheritance (missing ECIs → 0)
  - Component alphabetical ordering

---

### ⏳ **REMAINING WORK (Phases 2-5)**

#### **Step 8 ⏳ — CVMPhaseModel Integration (Phase 2)**
- **File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMPhaseModel.java`
- **What needs to change:**
  1. Add optional `CvCfIntegration` parameter to constructor
  2. Store both old and new cmat (for backward compatibility)
  3. Store ECIs as both flat array and Map<String, Double>
  4. Add method to detect/select CVCF vs old basis
  5. Pass correct cmat to NewtonRaphsonSolverSimple

- **Code location:** Lines ~200-275 (CVMPhaseModel constructor & field initialization)
- **Impact:** Enables CVCF to be used throughout solver without breaking old code
- **Status:** Blocked on this step

#### **Step 9 ⏳ — Enthalpy Integration (Phase 2)**
- **Files affected:**
  - `CVMFreeEnergy.java` (lines ~130-152: enthalpy computation)
  - `NewtonRaphsonSolverSimple.java` (uses CVMFreeEnergy)

- **What needs to change:**
  1. Detect whether to use old or CVCF enthalpy
  2. If CVCF: Extract CV values by name from cmat evaluation
  3. Call `CvCfHamiltonianEvaluator.computeH(cvMap, eciMap, pairs, ...)`
  4. If old: Keep existing H computation

- **Code location:** CVMFreeEnergy.evaluate() line ~130-152
- **Status:** Blocked on Step 8

#### **Step 10 ⏳ — Enthalpy Gradient (Phase 3)**
- **File:** `CVMFreeEnergy.java` (lines ~206-217: gradient computation)
- **What needs to change:**
  1. For CVCF: Compute dH/du via chain rule
  2. Formula: dH/du_k = Σ_{i,j} e_ij · d(v_ij)/du_k
  3. Since d(cv)/du comes from cmat (as chains), reuse existing gradient infrastructure

- **Key insight:** Hessian is unchanged (H is linear in CFs)
- **Status:** Blocked on Steps 8-9

#### **Step 11 ⏳ — ECI Loading from File (Phase 4) [Simplified with Option A]**
- **Files affected:** None (no changes required!)

- **Design Decision (Option A):**
  - ECIs in hamiltonian.json stay in **old orthogonal basis format** (unchanged)
  - No format changes to CECEntry or HamiltonianStore
  - CVCF transformations are internal to CVMPhaseModel only
  - MCS engine continues to use old-basis ECIs unchanged

- **Rationale:**
  - Simpler implementation (no format versioning needed)
  - Backward compatible (existing hamiltonian.json files work as-is)
  - Lower risk (no changes to I/O layer)
  - ECIs are physically the same, just represented differently per basis

- **Current hamiltonian.json format (STAYS UNCHANGED):**
  ```json
  {
    "structurePhase": "BCC_A2",
    "elements": "Nb-Ti",
    "cecTerms": [
      {"name": "CF_0", "a": 0.0, "b": 0.0},
      {"name": "CF_1", "a": -390.0, "b": 0.0}
    ]
  }
  ```

- **Status:** Blocked on Steps 8-9 (but scope is now zero-effort once those complete)

#### **Step 12 ⏳ — Verification Tests (Phase 5)**
- **What needs to verify:**
  1. **T matrix consistency:** Equimolar evaluation matches
  2. **CV values:** Old basis CV evaluation == CVCF basis CV evaluation
  3. **Full CVM run:** Binary Nb-Ti phase boundary unchanged
  4. **CEC inheritance:** Ternary with binary ECIs only produces continuous solution
  5. **Gradient accuracy:** Numerical vs analytical gradients match

- **Test location:** New test class `CvCfPhaseModelIntegrationTest.java`
- **Status:** Blocked on Steps 8-10

---

## Integration Dependency Chain

```
     COMPLETED ✅              REMAINING ⏳

     ┌─────────────────┐
     │ CvCfBasis       │
     │ & T Matrices    │
     │ (Steps 1-2) ✅  │
     └────────┬────────┘
              │
              ▼
     ┌─────────────────┐
     │ Transformer &   │
     │ Evaluator       │
     │ (Steps 3-5) ✅  │
     └────────┬────────┘
              │
              ▼
     ┌─────────────────┐
     │ CvCfIntegration │
     │ Adapter         │
     │ (Steps 6-7) ✅  │◄─── YOU ARE HERE
     └────────┬────────┘
              │
              ▼
     ┌──────────────────────────┐
     │ CVMPhaseModel &          │
     │ Enthalpy Integration     │
     │ (Steps 8-9) ⏳            │
     └────────┬─────────────────┘
              │
              ▼
     ┌──────────────────────────┐
     │ Gradient & Enthalpy      │
     │ Derivatives              │
     │ (Step 10) ⏳              │
     └────────┬─────────────────┘
              │
              ▼
     ┌──────────────────────────┐
     │ ECI Loading & Legacy     │
     │ Support                  │
     │ (Step 11) ⏳              │
     └────────┬─────────────────┘
              │
              ▼
     ┌──────────────────────────┐
     │ Verification & Testing   │
     │ (Step 12) ⏳              │
     └──────────────────────────┘
```

---

## Files Overview

### Core CVCF Classes (All Complete ✅)

| File | Lines | Status | Role |
|------|-------|--------|------|
| CvCfBasis.java | 80 | ✅ | Data holder for basis definition |
| BccA2CvCfTransformations.java | 380+ | ✅ | Hardcoded T matrices & CF names |
| CvCfBasisTransformer.java | 82 | ✅ | Applies transformation T to cmatrix |
| CvCfHamiltonianEvaluator.java | 90 | ✅ | Computes enthalpy via eq. 30 |
| CvCfCFRegistry.java | 75 | ✅ | CF name → index mapping |
| **CvCfIntegration.java** | **220** | **✅** | **[NEW] Orchestration layer** |

### Integration Points (Remaining ⏳)

| File | Lines | Change | Impact |
|------|-------|--------|--------|
| CVMPhaseModel.java | 200-275 | Add CVCF option | Step 8 blocker |
| CVMFreeEnergy.java | 130-217 | Route enthalpy & gradient | Steps 9-10 |
| NewtonRaphsonSolverSimple.java | All | Transparent (uses CVMFreeEnergy) | Auto-handled |
| HamiltonianStore.java | 80+ | Support CVCF CF names | Step 11 |

### Test Files

| File | Tests | Status |
|------|-------|--------|
| CvCfBasisTransformationTest.java | 11 | ✅ All passing |
| CvCfPhaseModelIntegrationTest.java | — | ⏳ To be created |

---

## Quick Reference: What Each Phase Does

| Phase | Steps | Duration | Blocking | Key Output |
|-------|-------|----------|----------|------------|
| **Phase 1** | 1-7 | ✅ Complete | Nothing | CvCfIntegration ready for use |
| **Phase 2** | 8-9 | ~2-3 days | None | CVMPhaseModel accepts CVCF |
| **Phase 3** | 10 | ~1 day | Phase 2 | Analytical gradients working |
| **Phase 4** | 11 | ~1 day | Phase 2 | File I/O supporting CVCF |
| **Phase 5** | 12 | ~2-3 days | Phases 3-4 | All tests passing, verified |

---

## How to Proceed to Phase 2

**Next immediate task:**

1. Open `CVMPhaseModel.java`
2. Add this to the constructor:
   ```java
   private final CvCfIntegration cvCfIntegration;  // null if using old basis
   ```
3. Modify `minimize()` to:
   - Use `cvCfIntegration.getTransformedCmat()` if CVCF is active
   - Otherwise use `cmat` (old basis)
4. Pass correct cmat to `NewtonRaphsonSolverSimple.solve()`
5. Update ECI handling to use `integration.mapEciToFlat()` if needed

This unblocks Phase 3 (enthalpy integration) and Phase 4 (ECI loading).

---

## Summary Table

```
┌──────────────┬────────────────┬──────────────┬──────────────┐
│ Phase        │ Completion     │ Status       │ Next Action  │
├──────────────┼────────────────┼──────────────┼──────────────┤
│ Phase 1      │ ✅ 100% (7/7)  │ COMPLETE     │ [Done]       │
│ Phase 2      │ ⏳ 0% (0/2)    │ READY        │ Start now    │
│ Phase 3      │ ⏳ 0% (0/1)    │ BLOCKED      │ After Ph 2   │
│ Phase 4      │ ⏳ 0% (0/1)    │ BLOCKED      │ After Ph 2   │
│ Phase 5      │ ⏳ 0% (0/1)    │ BLOCKED      │ After Ph 3-4 │
├──────────────┼────────────────┼──────────────┼──────────────┤
│ **OVERALL**  │ **✅ 60% (7/10)** | **ON TRACK** | **PHASE 2** |
└──────────────┴────────────────┴──────────────┴──────────────┘
```

**Overall Progress:** 6 of 10 steps complete (60%)
**Estimated Time to Full Completion:** 5-7 days (assuming daily work)
**Critical Path:** Phase 2 → Phase 3 → Phase 5 (Phases 4 independent)

