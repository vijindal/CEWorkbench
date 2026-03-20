# Cluster Data Structure: Four-Layer Architecture

## Overview

The cluster identification workflow (Type-1a) generates cluster data in a **four-layer crystallographic cascade**. Each layer builds upon the previous, establishing symmetry foundations before introducing decoration-specific details.

---

## The Four Layers

### Layer 1: HSP Clusters (Highest Symmetric Phase)

**Example:** A2 phase (for B2 ternary system)

- **Symmetry:** Highest possible symmetry
- **Decoration:** Undecorated — all sites use only `s[1]` (site occupation) operators
- **Purpose:** Foundation for all higher-order systems
- **Output:**
  - Multiplicity information for each unique cluster
  - KB (Korringa-Kohn-Rostoker) coefficients
  - Cluster topology and correlation basis

**Why it matters:** KB coefficients derived from the HSP are invariant across all phases in the system. Computing them once at the highest symmetry avoids redundant calculations.

---

### Layer 2: Given Phase Clusters

**Example:** B2 ternary (actual phase being studied)

- **Symmetry:** Lower than HSP, but cluster topology identical to HSP
- **Decoration:** Undecorated at this stage (equivalent to binary structure in terms of cluster geometry)
- **Purpose:** Establishes clusters for the specific phase
- **Relationship to Layer 1:** Same clusters, same multiplicity, different symmetry context

**Why it matters:** Validates that the target phase's cluster set is consistent with the HSP foundation.

---

### Layer 3: HSP Correlation Functions

**Example:** A2 ternary CFs (for B2 ternary system)

- **Basis:** Correlation functions derived from HSP clusters
- **Decoration:** Undecorated (simpler CF basis)
- **Purpose:** Reference basis functions from the highest-symmetry phase
- **Output:**
  - List of unique correlation functions
  - Scalar summaries (count, average distance, etc.)
  - Grouping structure

**Why it matters:** CFs are organized and indexed relative to the HSP phase. This provides a consistent reference frame across all phases in the system.

---

### Layer 4: Given Phase Correlation Functions

**Example:** B2 ternary CFs (actual phase)

- **Basis:** Correlation functions for the target phase
- **Decoration:** Fully decorated according to number of components
- **Purpose:** Feature space for the CE (Cluster Expansion) model
- **Output:**
  - Component-dependent correlation functions
  - CFs grouped by corresponding HSP clusters
  - All information needed for CE enthalpy calculations

**Why it matters:** These are the actual basis functions used in thermodynamic modeling. They must be grouped by HSP clusters to maintain structural and physical consistency.

---

## Workflow Sequence

```
Input: Crystal structure (target phase) + HSP definition
  ↓
Step 1: Identify HSP clusters → Layer 1
  ├─ Extract multiplicity
  └─ Compute KB coefficients
  ↓
Step 2: Identify given-phase clusters → Layer 2
  └─ Validate consistency with Layer 1
  ↓
Step 3: Generate HSP correlation functions → Layer 3
  └─ Build reference CF basis
  ↓
Step 4: Generate given-phase CFs → Layer 4
  ├─ Decorate CFs by component
  ├─ Group by HSP clusters
  └─ Ready for thermodynamic calculations
  ↓
Output: All four layers stored in cluster_data.json
```

---

## Why All Four Are Necessary During Workflow

1. **HSP establishes invariants** — KB coefficients computed once, used for all phases
2. **Phase clusters validate topology** — confirms target phase is compatible with HSP
3. **HSP CFs provide reference structure** — organizing principle for grouping lower-phase CFs
4. **Phase CFs are the actual feature space** — used by CVM and MCS engines

Without this sequence, you risk:
- Computing KB coefficients incorrectly for multi-component systems
- Losing the structural relationship between phases
- Having disorganized or redundant correlation functions

---

## Runtime Usage vs. Storage

**At runtime (CVM/MCS engines):**
- Only Layer 4 (given-phase CFs) is needed
- Layer 4 CF scalar summaries (tcf, ncf, lcf, etc.)
- Layer 1 KB coefficients
- Cluster multiplicity and site data

**Not read at runtime:**
- Layer 1 full cluster data (HSP clusters)
- Layer 2 given-phase clusters (topology already in Layer 4)
- Layer 3 correlation functions (HSP CFs)
- Layer 4 raw CF geometry

---

## Storage Optimization

Because Layers 1–3 are workflow intermediates and not needed at runtime, they are **excluded from serialization** using `@JsonIgnore` annotations. See [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) for details on the ~96% file size reduction.

The workflow still generates all four layers internally for correctness and validation. Only the runtime-essential data is persisted.

---

## Example: B2 Ternary System

| Layer | Phase | Clusters | Decorated? | KB Coeff? | Runtime Used? |
|-------|-------|----------|-----------|-----------|---------------|
| 1 | A2 (HSP) | High-sym parent | No | **Yes** | No (stored for ref) |
| 2 | B2 | Given phase | No | No | No (transient) |
| 3 | A2 (HSP) | High-sym parent | No | No | No (transient) |
| 4 | B2 | Given phase | **Yes** | No | **Yes** |

The actual CE enthalpy model uses **only Layer 4 CFs**, but Layer 4 was correctly constructed by validating against Layers 1–3.

---

## References

- [SIMPLIFIED_ARCHITECTURE.md](SIMPLIFIED_ARCHITECTURE.md) — overall system architecture
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — file size reduction strategy
