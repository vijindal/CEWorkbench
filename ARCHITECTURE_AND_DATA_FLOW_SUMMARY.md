# CE Thermodynamics Workbench: Complete Architecture & Data Flow Summary

## The Two Fundamental Tasks

The CE Thermodynamics Workbench performs **two classes of work**:

1. **Type-1a: Cluster Identification Workflow**
   - Generates reusable scientific data from crystal structure
   - Input: atomic positions, symmetry, composition
   - Output: `cluster_data.json` containing clusters, correlations, C-matrix

2. **Type-2: Thermodynamic Calculation**
   - Computes equilibrium phase diagrams using CE model
   - Input: `cluster_data.json` + Hamiltonian (ECI)
   - Output: phase diagram, thermodynamic properties at any (T, composition)

---

## Type-1a: Cluster Identification Workflow (Data Generation)

### Five Sequential Computation Stages

```
Input: Crystal structure, symmetry ops, site operators, components
  │
  ├─ STAGE 1: HSP (Highest Symmetric Phase) Clusters
  │   └─ genClusCoordList[disMaxClusCoord, disSymOpList, ...]
  │   └─ getNijTable, generateKikuchiBakerCoefficients
  │   Output: mhdis[itc], kbdis[itc], nijTable
  │
  ├─ STAGE 2: Given-Phase Clusters
  │   └─ genClusCoordList[maxClusCoord, symOpList, ...]
  │   └─ Coordinate transform, transClusCoordList
  │   Output: mh[itc][inc], lc[itc], (raw coords)
  │
  ├─ STAGE 3: HSP Correlation Functions
  │   └─ genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolList]
  │   Output: tcfdis, mcfdis, (CF coords)
  │
  ├─ STAGE 4: Given-Phase CFs (with grouping)
  │   └─ genClusCoordList, ordToDisordCoord, transClusCoordList
  │   └─ groupCFData[disClusData, disCFData, ordCFData, ...]  ← KEY STEP
  │   Output: cfCoordList, lcf[itc][inc], tcf, mcf, rcf
  │
  └─ C-MATRIX STAGE: Correlation Matrix
      └─ genSiteList, genCfSiteOpList, genSubstituteRules
      └─ genPRules (conventional basis)
      └─ genCV[ordClusData, cfData, ...]  ← Depends on Stages 1,2,4
      Output: cmat, lcv[itc][inc], wcv[itc][inc][incv], cv
```

### Storage Optimization

After generation, data is selectively serialized to `cluster_data.json`:

| Stage | Output | Serialized? | Reason |
|-------|--------|-------------|--------|
| 1 | `mhdis`, `kbdis`, `nijTable` | ✓ YES | Essential for entropy/enthalpy |
| 2 | `mh`, `lc`, (raw coords) | ✓/✗ Partial | Keep scalars, drop coordinates |
| 3 | `tcfdis`, `mcfdis`, (CF coords) | ✗ NO | Workflow-only, used for grouping |
| 4 | `cfCoordList`, `lcf`, `tcf`, `mcf`, `rcf` | ✓ YES | Defines feature space |
| C | `cmat`, `lcv`, `wcv`, `cv` | ✓ YES | Expensive, needed at runtime |

**Result:** ~96% file size reduction (345 KB → 13.6 KB)

---

## Type-2: Thermodynamic Calculation Workflow (Using Generated Data)

### CVM Engine Flow

```
Load cluster_data.json
  ├─ mhdis[i], kbdis[i]           ← Stage 1
  ├─ mh[itc][inc], lc[itc]        ← Stage 2
  ├─ lcf[itc][inc]                ← Stage 4
  ├─ cmat, lcv, wcv, cv           ← C-Matrix
  │
Load Hamiltonian (CE model)
  ├─ eList[i] (ECI parameters)
  │
For each (composition xB, temperature T):
  │
  ├─ Solve CVM equations via Lagrange multipliers
  │   └─ Find occupation variables {σᵢ, σᵢσⱼ, ...}
  │
  ├─ Compute correlations
  │   └─ uList[i] = cmat · occupations
  │   └─ uA[i], uB[i] = reference states
  │
  ├─ ENTROPY (Sc)
  │   └─ Sc = -R * Sum[kbdis[itc] * mhdis[itc] * ...]
  │       Uses: mhdis, kbdis, mh, wcv, uList, lcv
  │
  ├─ ENTHALPY (Hc)
  │   └─ Hc = Sum[mhdis[i] * eList[i] * (uList[i] - reference[i])]
  │       Uses: mhdis, eList, uList, uA, uB
  │
  ├─ GIBBS ENERGY (Gc)
  │   └─ Gc = Hc - T*Sc
  │       Uses: H, S, T
  │
  └─ Output: EquilibriumState {Gc, Hc, Sc, properties}
```

---

## System Architecture (Four Layers)

```
            ┌─────────────────────────────────────┐
            │           UI Layer                  │
            │  ├─ GUI (Swing/JavaFX panels)      │
            │  └─ CLI (command-line interface)   │
            └────────────────┬────────────────────┘
                             │
            ┌────────────────┴────────────────────┐
            │      Workflow Layer                 │
            │  ├─ CalculationService             │
            │  ├─ ClusterIdentificationWorkflow  │
            │  └─ ThermodynamicWorkflow          │
            └────────────────┬────────────────────┘
                             │
            ┌────────────────┴────────────────────┐
            │      Domain Layer (Physics)         │
            │  ├─ cluster package                 │
            │  │   ├─ Cluster, ClusterData       │
            │  │   ├─ CorrelationBasis, CMatrix │
            │  │   └─ Site, Sublattice          │
            │  ├─ hamiltonian package            │
            │  │   ├─ Hamiltonian               │
            │  │   └─ CEHamiltonian (ECI)       │
            │  ├─ engine package                 │
            │  │   ├─ ThermodynamicEngine       │
            │  │   ├─ CVMEngine                 │
            │  │   └─ MCSEngine                 │
            │  └─ result package                 │
            │      ├─ EquilibriumState          │
            │      └─ EngineMetrics             │
            └────────────────┬────────────────────┘
                             │
            ┌────────────────┴────────────────────┐
            │     Storage Layer (Disk I/O)       │
            │  ├─ ClusterDataStore               │
            │  ├─ HamiltonianStore               │
            │  ├─ Workspace                      │
            │  └─ JSON serialization             │
            └─────────────────────────────────────┘
```

**Golden Rule:** Domain ← everything else (nothing points backward)

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         USER INPUT                              │
│  Crystal structure, symmetry, components, number of sites      │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
            ┌──────────────────────────────┐
            │   TYPE-1a WORKFLOW           │
            │ (Cluster Identification)     │
            │                              │
            │ Five stages:                 │
            │ 1. HSP clusters              │
            │ 2. Given-phase clusters      │
            │ 3. HSP CFs                   │
            │ 4. Given-phase CFs + group   │
            │ 5. C-matrix computation      │
            └──────────────────┬───────────┘
                               │
                    ┌──────────▼──────────┐
                    │ Optimize & Serialize│
                    │ cluster_data.json   │
                    └──────────┬──────────┘
                               │
        ┌──────────────────────┴──────────────────────┐
        │                                             │
        ▼                                             ▼
┌────────────────────┐                    ┌─────────────────────┐
│  cluster_data.json │                    │ CE Hamiltonian      │
│  (~13.6 KB)        │                    │ (fitted ECI params) │
│                    │                    │                     │
│ ├─ mhdis, kbdis   │                    │ ├─ eList[i]        │
│ ├─ mh, lc         │                    │ └─ Symmetry info   │
│ ├─ cfCoordList    │                    └─────────────────────┘
│ ├─ lcf, tcf, mcf  │
│ └─ cmat, lcv, wcv │
└────────────────────┘
        │
        │  Load
        ▼
    ┌────────────────────────────────┐
    │  TYPE-2 WORKFLOW               │
    │ (Thermodynamic Calculation)    │
    │                                │
    │ For each (T, xB):             │
    │ 1. Solve CVM equations         │
    │ 2. Compute correlations (uList)│
    │ 3. Calculate entropy (Sc)      │
    │ 4. Calculate enthalpy (Hc)     │
    │ 5. Calculate Gibbs (Gc)        │
    └────────────────┬───────────────┘
                     │
                     ▼
            ┌─────────────────────┐
            │ Thermodynamic State │
            │ {Gc, Hc, Sc, u[i]} │
            │ @(T, xB)            │
            └────────┬────────────┘
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
   ┌─────────────┐        ┌────────────────┐
   │Phase Diagram│        │Thermodynamic   │
   │             │        │Properties Plot │
   │T vs xB      │        │ (H, S, Cp...)  │
   └─────────────┘        └────────────────┘
```

---

## Key Data Structures and Their Sources

### From Type-1a (Cluster Data)

| Variable | Stage | Type | Usage | Size |
|----------|-------|------|-------|------|
| `mhdis[itc]` | 1 | Vector | Entropy, enthalpy weight | small |
| `kbdis[itc]` | 1 | Vector | Entropy coefficient | small |
| `nijTable` | 1 | Matrix | KB recomputation (rare) | small |
| `mh[itc][inc]` | 2 | 2D array | Entropy scaling | small |
| `lc[itc]` | 2 | Vector | Loop bounds | small |
| `lcf[itc][inc]` | 4 | 2D array | CF organization | small |
| `tcf` | 4 | Scalar | Total CF count | 1 |
| `cmat` | C-Matrix | Matrix | CF computation | medium |
| `lcv[itc][inc]` | C-Matrix | 2D array | CV loop bounds | small |
| `wcv[itc][inc][incv]` | C-Matrix | 3D array | CF weights | small |

### From External Input (CE Model)

| Variable | Role | Source |
|----------|------|--------|
| `eList[i]` | ECI energy parameters | Hamiltonian (fitted) |
| `xB` | Composition | User input |
| `T` | Temperature | User input |

### Computed at Runtime

| Variable | Computed From | Used In |
|----------|---------------|---------|
| `occupations` | CVM Lagrange equations | CF computation |
| `uList[i]` | `cmat · occupations` | Entropy, enthalpy |
| `uA[i]`, `uB[i]` | Pure-phase occupations | Enthalpy reference |
| `Sc` | Entropy formula | Gibbs energy |
| `Hc` | Enthalpy formula | Gibbs energy |
| `Gc` | `Hc - T·Sc` | Phase stability |

---

## Critical Dependencies

### Why All Type-1a Stages Are Necessary

1. **Stage 1 → KB coefficients** — invariant across all phases
2. **Stage 2 → Multiplicity ratios** — given-phase structure mapping
3. **Stage 3 → CF reference basis** — defines feature space structure
4. **Stage 4 grouping** — requires Stages 1, 2, and 3 to organize CFs correctly
5. **C-Matrix** — requires Stages 1, 2, and 4 to project occupations to CFs

Each stage is a prerequisite for the next. Removing any would break the correctness of the final data.

### Why Stages 2 & 3 Can Be Removed From Storage

- Stage 2 raw coordinates → only scalars (`mh`, `lc`) used at runtime ✓ Keep scalars
- Stage 3 CF coordinates → only used for grouping in Stage 4 ✓ Done at generation

### Why Stages 1, 4, & C-Matrix Must Be Kept

- Stage 1 → `mhdis`, `kbdis` essential for entropy/enthalpy
- Stage 4 → `cfCoordList`, `lcf`, `mcf` define feature space
- C-Matrix → exponentially expensive to recompute, needed at every iteration

---

## Thermodynamic Model Equations

### Configurational Entropy (CVM)

```
Sc = -R * Σᵢₜc Σᵢₙc Σᵢₙᵥ (
  kbdis[itc] * mhdis[itc] * mh[itc][inc] *
  wcv[itc][inc][incv] * u[itc][inc][incv] * Log[u[itc][inc][incv]]
)
```

### Configurational Enthalpy (CE)

```
Hc = Σᵢ (
  mhdis[i] * eList[i] * (uList[i] - reference[i])
)

where:
reference[i] = (1 - xB) * uA[i] + xB * uB[i]
```

### Gibbs Free Energy

```
Gc = Hc - T * Sc
```

---

## File Organization

```
CEWorkbench/
├─ SIMPLIFIED_ARCHITECTURE.md           ← System overview
├─ CLUSTER_DATA_STRUCTURE.md            ← Why 4 layers
├─ CLUSTER_DATA_GENERATION_FLOW.md      ← How 4 stages (Mathematica)
├─ CMATRIX_CALCULATION_FLOW.md          ← C-matrix computation
├─ CVM_ENTROPY_CALCULATION.md           ← Entropy formula usage
├─ CVM_THERMODYNAMIC_PROPERTIES.md      ← Full H, S, Gc calculation
├─ COMPLETE_DATA_FLOW.md                ← Type-1a to Type-2 flow
├─ CLUSTER_DATA_STORAGE_OPTIMIZATION.md ← Implementation details
├─ README_CLUSTER_DATA.md               ← Quick reference
├─ ARCHITECTURE_AND_DATA_FLOW_SUMMARY.md (this file)
│
├─ src/main/java/org/ce/
│  ├─ domain/                    (Physics, algorithms)
│  │  ├─ cluster/
│  │  ├─ hamiltonian/
│  │  ├─ engine/
│  │  └─ result/
│  ├─ workflow/                  (Orchestration)
│  ├─ storage/                   (Disk I/O)
│  └─ ui/                        (GUI/CLI)
│
├─ cluster_data.json             (Generated: Stages 1-5 output)
└─ hamiltonian_data.json         (Input: CE model, ECI params)
```

---

## Reading Order

1. **Start here:** [SIMPLIFIED_ARCHITECTURE.md](SIMPLIFIED_ARCHITECTURE.md) — System design
2. **Understand data:** [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — Why 4 layers
3. **See implementation:** [CLUSTER_DATA_GENERATION_FLOW.md](CLUSTER_DATA_GENERATION_FLOW.md) — Mathematica code
4. **Learn usage:** [CVM_ENTROPY_CALCULATION.md](CVM_ENTROPY_CALCULATION.md) — Runtime entropy
5. **Complete picture:** [COMPLETE_DATA_FLOW.md](COMPLETE_DATA_FLOW.md) — Type-1a to Type-2

For implementation details:
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — Java changes needed
- [CVM_THERMODYNAMIC_PROPERTIES.md](CVM_THERMODYNAMIC_PROPERTIES.md) — Hc and Gc calculations

---

## Key Design Principles

1. **Domain drives architecture** — physics first, enterprise patterns second
2. **Minimize abstractions** — no repositories, factories, or adapters; just data stores
3. **Correctness before optimization** — generate all stages for validation, optimize storage
4. **Clean separation** — domain independent of workflow, storage, UI
5. **Scientific organization** — structure mirrors crystallographic symmetry

---

## Next Steps for Implementation

### For Java Developers

- [ ] Apply `@JsonIgnore` annotations to Stage 2 & 3 raw coordinates (per optimization plan)
- [ ] Add `@JsonIgnoreProperties(ignoreUnknown = true)` for backward compatibility
- [ ] Verify CVM/MCS still run correctly with optimized JSON
- [ ] Test loading old cluster_data.json files
- [ ] Measure new file size: target ~13.6 KB

### For Science Developers

- [ ] Verify Stage 4 grouping produces correct CF organization
- [ ] Check C-matrix dimensions and weights
- [ ] Validate entropy formula matches CVM theory
- [ ] Test enthalpy calculation against first-principles data
- [ ] Verify phase diagrams against experimental data

---

## References

All documentation files are referenced within for deep dives into specific topics.

