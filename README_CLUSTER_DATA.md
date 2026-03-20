# Cluster Data: Structure & Storage Guide

## Quick Overview

The CE Thermodynamics Workbench generates and stores **cluster data** in a four-layer crystallographic cascade. This document explains what the layers are, why they're needed, and how we optimize storage.

---

## Two Key Documents

### 1. [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md)
**Read this to understand WHY we generate four layers**

- **Layer 1:** HSP (Highest Symmetric Phase) Clusters
- **Layer 2:** Given Phase Clusters
- **Layer 3:** HSP Correlation Functions
- **Layer 4:** Given Phase Correlation Functions

Each layer builds on the previous to establish symmetry foundations and decorations correctly. All four are generated during the Type-1a workflow.

### 2. [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md)
**Read this to understand HOW we reduce file size by 96%**

- Layers 1–3 are workflow intermediates, not needed at runtime
- We serialize only Layer 4 (and KB coefficients from Layer 1)
- ~96% file size reduction: 345 KB → 13.6 KB
- Implementation: `@JsonIgnore` annotations + backward compatibility

---

## The Paradox (Resolved)

**Question:** Why store all four layers if only Layer 4 is used at runtime?

**Answer:**
- We *generate* all four layers for **correctness** (symmetry validation, KB coefficient derivation, CF grouping)
- We *store* only what's needed for **runtime** (Layer 4 CFs + KB coefficients)
- The workflow discards transient layers after validation

---

## File Structure

```
cluster_data.json
├─ disorderedClusterResult (Layer 1 + 2 combined)
│  ├─ disClusterData          ← Layer 1: HSP cluster data [runtime: multiplicity, KB coeff]
│  └─ ordClusterData          ← Layer 2: Given phase clusters [transient: validation only]
│
├─ disorderedCFResult (Layer 3 + 4 combined)
│  ├─ disCFData               ← Layer 3: HSP CFs [transient: reference only]
│  └─ lcf, tcf, ncf, ...      ← Layer 4: Given phase CF scalars [runtime: used by CVM/MCS]
│
└─ cmatrixResult              ← C-Matrix (expensive to recompute, always kept)
```

---

## Runtime Flow

```
user runs CVM/MCS calculation
  ↓
CalculationService loads cluster_data.json
  ↓
ThermodynamicEngine reads:
  • Cluster multiplicity (Layer 1)
  • KB coefficients (Layer 1)
  • CF scalars: lcf, tcf, ncf (Layer 4)
  • C-Matrix
  ↓
Computes equilibrium thermodynamics
  ↓
Returns EquilibriumState
```

Layers 2 and 3 are **never accessed**. Optimization removes them entirely from storage.

---

## Workflow Generation Flow

```
Type-1a: Cluster Identification Workflow
  ↓
Step 1: Identify HSP clusters
  └─ Output: Layer 1
  ↓
Step 2: Identify given-phase clusters
  └─ Output: Layer 2
  ├─ Validate against Layer 1
  └─ Keep during workflow for validation
  ↓
Step 3: Generate HSP correlation functions
  └─ Output: Layer 3
  ├─ Build reference basis
  └─ Keep during workflow for CF grouping
  ↓
Step 4: Generate given-phase correlation functions
  └─ Output: Layer 4
  ├─ Decorate by component
  ├─ Group by HSP clusters
  └─ **This is what gets saved**
  ↓
Before saving: @JsonIgnore prevents Layers 2–3 serialization
  ↓
Result: cluster_data.json (96% smaller)
```

---

## Key Design Principles

1. **Correctness over storage:** We generate all four layers to ensure proper symmetry handling
2. **Minimize persistence:** We store only what's needed at runtime
3. **Backward compatibility:** Old JSON files with all four layers still load correctly (`@JsonIgnoreProperties(ignoreUnknown = true)`)
4. **No computation changes:** All optimization is via annotations; no algorithms are modified

---

## Implementation Status

✅ **Completed:**
- Architecture documented
- Optimization plan designed (96% reduction target)
- `@JsonIgnore` annotations identified for specific getters
- Backward compatibility strategy in place

**Next steps:**
- Apply annotations to domain classes
- Run Type-1a workflow and verify new JSON size
- Verify CVM/MCS still work correctly with optimized data
- Test loading old cluster_data.json files

---

## For Developers

### When working on cluster identification (Type-1a)

Ensure you understand:
- Why both HSP and given-phase clusters are computed (validation, KB coefficients)
- Why both HSP and given-phase CFs are generated (grouping, reference structure)
- That Layers 2–3 are **transient** — validate during generation, discard before saving

### When working on thermodynamic engines (CVM/MCS)

Remember:
- You only need Layer 4 CF data at runtime
- Multiplicity and KB coefficients come from Layer 1
- All other layers are pre-computation artifacts

### When optimizing storage

See [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) for:
- Which getters to `@JsonIgnore`
- How to add backward compatibility
- Projected file size reductions per section

---

## References

- [SIMPLIFIED_ARCHITECTURE.md](SIMPLIFIED_ARCHITECTURE.md) — overall system design
- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — detailed layer explanation
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — technical optimization plan
