# Dual-Pipeline Architecture: Old Basis vs CVCF Basis

**Decision Date:** March 25, 2026
**Architecture Pattern:** Optional feature with backward compatibility

---

## 1. Why Two Pipelines?

- **CVCF basis** offers mathematically invariant representation (Jindal & Lele 2025)
- **Old orthogonal basis** is current production-stable code
- **Both are equivalent** — same thermodynamic properties, different internal representation
- **Gradual adoption** — CVCF is opt-in, not forced

---

## 2. How They Coexist

### Shared Layer (Identical)
```
┌─────────────────────────────────────────────┐
│  Input / Output                              │
├─────────────────────────────────────────────┤
│ • hamiltonian.json (old basis ECIs)         │
│ • Cluster data (Stage 1-3, unchanged)       │
│ • Composition, Temperature inputs           │
├─────────────────────────────────────────────┤
│ Output: EquilibriumState                    │
│ {G, H, S, composition, T, phase}            │
│ (Basis-invariant — identical either way)    │
└─────────────────────────────────────────────┘
```

### Computation Branching (Inside CVMPhaseModel)
```
CVMPhaseModel constructor receives:
  • CMatrixResult cmat_old (always)
  • CvCfIntegration (optional, nullable)

if (cvcfIntegration != null) {
    // CVCF PATH
    cmat_new = cvcfIntegration.getTransformedCmat()
    use cmat_new for CV computation
    use CvCfHamiltonianEvaluator for H
} else {
    // OLD PATH (default)
    cmat = cmat_old
    use existing H computation
}

// Both paths then:
compute S via Kikuchi-Barker (unchanged)
compute G = H - T·S
minimize via Newton-Raphson
```

---

## 3. MCS Engine (Unchanged)

```
MCSEngine {
  input.cec.cecTerms  ← Always old basis ECIs from hamiltonian.json
  input.clusterData   ← Cluster coordination lists

  MCSRunner.run() {
    LocalEnergyCalc.totalEnergy(config, emb, eci, orbits)
    // eci is indexed by cluster type (integer)
    // No basis awareness needed
  }

  EquilibriumState output
}
```

**Key:** MCS gets ECI values from the same hamiltonian.json file as CVM. No format change, no conversion needed. MCS remains completely agnostic to basis.

---

## 4. File Locations and Scope

| File | Change? | Reason |
|------|---------|--------|
| hamiltonian.json | ❌ No | ECIs stay in old basis format |
| HamiltonianStore.java | ❌ No | Reads unchanged JSON format |
| CECEntry.java | ❌ No | No new fields needed |
| MCSEngine.java | ❌ No | Uses old-basis ECIs unchanged |
| LocalEnergyCalc.java | ❌ No | No basis awareness |
| MCSRunner.java | ❌ No | No changes |
| CVMPhaseModel.java | ✅ **Yes** (Phase 2) | Add optional CvCfIntegration param |
| CVMFreeEnergy.java | ✅ **Yes** (Phase 2-3) | Route H/dH computation based on basis |
| NewtonRaphsonSolverSimple.java | ⚠️ Minimal (Phase 2) | Just receives different cmat, logic unchanged |
| CvCfIntegration.java | ✅ **New** (Phase 1 ✅) | Already created |
| BccA2CvCfTransformations.java | ✅ **Complete** (Phase 1 ✅) | T matrices, CF names |
| CvCfBasisTransformer.java | ✅ **Complete** (Phase 1 ✅) | C-matrix transformation |

---

## 5. Data Flow Comparison

### OLD BASIS PATH
```
hamiltonian.json (old basis)
        ↓
   eci[0..n] (indexed by cluster type)
        ↓
  CVMPhaseModel (no CvCfIntegration)
        ├─ cmat_old
        ├─ eci[] unchanged
        └─ Standard H computation: H = Σ eci · u
        ↓
  EquilibriumState
```

### CVCF BASIS PATH
```
hamiltonian.json (old basis)
        ↓
   eci[0..n] (indexed by cluster type)
        ↓
  CvCfIntegration
        ├─ Transform: cmat_old → cmat_new
        └─ Keep: eci[] as-is
        ↓
  CVMPhaseModel (with CvCfIntegration)
        ├─ cmat_new (CVCF basis)
        ├─ eci[] unchanged (still old basis indices)
        └─ CVCF H computation: via CvCfHamiltonianEvaluator
           (converts CV values to CF names internally)
        ↓
  EquilibriumState
  (same properties as old basis path)
```

---

## 6. ECI Index Mapping

**Old Basis Indexing** (current, unchanged):
```
eci[0] → CF_0
eci[1] → CF_1
eci[2] → CF_2
...
eci[n-1] → CF_n-1
```

**CVCF Basis** (internal only, not stored):
```
e4NbTi  ← internally computed from eci[] via C-matrix evaluation
e3NbTi
e21NbTi
...
```

The conversion happens **inside CVMPhaseModel** via C-matrix evaluation, not via file I/O.

---

## 7. Phase Breakdown with Option A

| Phase | CVM | MCS | I/O | Status |
|-------|-----|-----|-----|--------|
| **Phase 0** (before) | Old basis only | Uses old ECIs | Old format | ✅ Current state |
| **Phase 1** | + CVCF infra ready | (unchanged) | (unchanged) | ✅ COMPLETE |
| **Phase 2** | ± CvCfIntegration optional | (unchanged) | (unchanged) | ⏳ In progress |
| **Phase 3** | ± Enthalpy gradients | (unchanged) | (unchanged) | ⏳ Pending |
| **Phase 4** | (no changes) | (unchanged) | (no changes) | ✅ ZERO EFFORT |
| **Phase 5** | ± Verification tests | (unchanged) | (unchanged) | ⏳ Pending |

---

## 8. Risk Analysis

### Low Risk
- ✅ MCS completely isolated from CVCF changes
- ✅ I/O format unchanged
- ✅ Old basis path untouched (can revert if needed)
- ✅ Both compute identical thermodynamic properties

### Moderate Risk
- ⚠️ CVMPhaseModel dual cmat handling (but minimal logic change)
- ⚠️ CVMFreeEnergy routing (but clean switch)

### Mitigation
- Run both CVM paths on same system, verify H/S/G match
- Phase 5 has explicit equivalence tests
- Rollback easy: just don't use CvCfIntegration param

---

## 9. Future Evolution (Optional)

If CVCF becomes preferred in future:
- Move hamiltonian.json to CVCF basis format
- Add `basis` field: `"basis": "CVCF"` or `"basis": "OLD"`
- Inverse-transform old→new at load time if needed
- **But this is not required for Phase 2-5 completion**

---

## 10. Summary: Option A Benefits

| Benefit | Why This Matters |
|---------|------------------|
| **No hamiltonian.json changes** | Existing data stays valid indefinitely |
| **MCS stays unchanged** | No new bugs, no revalidation needed |
| **CVCF is opt-in** | Old code paths guaranteed stable |
| **Basis-invariant output** | Both pipelines produce identical results |
| **Lower complexity** | Fewer moving parts = fewer bugs |
| **Faster development** | Fewer phases, fewer integration points |
| **Easier testing** | Comparison path is same code (old basis) |

