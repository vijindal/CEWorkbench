# CVM Entropy Calculation: Using Cluster Data at Runtime

## The Entropy Expression

```mathematica
Sc = -R * Sum[
  kbdis[[itc]] * msdis[[itc]] *
  Sum[
    mh[[itc]][[inc]] *
    Sum[
      wcv[[itc]][[inc]][[incv]] * CV[itc][inc][incv] * Log[CV[itc][inc][incv]],
      {incv, 1, lcv[[itc]][[inc]]}
    ],
    {inc, 1, lc[[itc]]}
  ],
  {itc, 1, tcdis}
];

Sc = Sc /. Table[msdis[[i]] -> mhdis[[i]], {i, tcdis}]
```

This is the **configurational entropy** computed by the CVM engine at each thermodynamic state point.

---

## Breakdown: What Data Comes From Where

### Loop Structure

The triple-nested sum iterates through a **hierarchy**:

```
For each HSP cluster type (itc = 1 to tcdis)
  For each cluster of that type (inc = 1 to lc[[itc]])
    For each correlation function in that cluster (incv = 1 to lcv[[itc]][[inc]])
      Compute: KB[itc] × m[itc] × mh[itc,inc] × w[itc,inc,incv] × CV[itc,inc,incv] × Log[CV]
```

This hierarchy reflects **exactly** how the cluster data was organized during Type-1a generation.

---

## Source of Each Variable

### **kbdis[[itc]]** — Kikuchi-Baker Coefficient

**From:** Stage 1 (HSP Clusters)
**Source class:** `ClusterIdentificationResult.disorderedClusterResult.kbCoefficients`
**Role:** Cluster-type weight in entropy summation

- One KB coefficient per HSP cluster type
- Derived from `nijTable` in Stage 1
- **Same across all phases** in the system
- Example: `kbdis = {1, -1/2, -1/2, 1/6, ...}`

---

### **msdis[[itc]]** — Disordered Phase Multiplicity

**From:** Stage 1 (HSP Clusters)
**Source class:** `ClusterIdentificationResult.disorderedClusterResult.disClusterData.multiplicities`
**Role:** Initial scaling factor (later replaced by mhdis)

- One multiplicity per HSP cluster type
- Count of equivalent clusters in the structure

**Note:** Immediately substituted with `mhdis` in the final step:
```mathematica
Sc = Sc /. Table[msdis[[i]] -> mhdis[[i]], {i, tcdis}]
```

---

### **mh[[itc]][[inc]]** — Multiplicity Ratio

**From:** Stage 2 (Given-Phase Clusters)
**Source class:** `ClusterIdentificationResult.orderedClusterResult.ordClusterData.multiplicities` (if stored)
**Computed as:** `mh = ordClusMList / mhdis`
**Role:** Cluster-specific scaling

- Two-level index: `[cluster type][cluster within type]`
- Ratio of given-phase multiplicity to HSP multiplicity
- Example: for B2 ternary vs. A2 parent: `mh[[1]] = {1, 1, 0.5, ...}`

---

### **lcv[[itc]][[inc]]** — CF Count

**From:** C-Matrix Calculation
**Source class:** `CMatrixResult.lcv` (length of correlation vectors)
**Role:** Loop limit for inner sum

- Two-level index: `[cluster type][cluster within type]`
- Number of distinct CFs for cluster `inc` of type `itc`
- Example: `lcv[[1]][[1]] = 5` means 5 CFs for point cluster in class 1

---

### **wcv[[itc]][[inc]][[incv]]** — CF Weight

**From:** C-Matrix Calculation
**Source class:** `CMatrixResult.wcv` (weights of correlation vectors)
**Role:** Normalization/scaling for each CF

- Three-level index: `[cluster type][cluster within type][CF index]`
- Positive weight for each CF in the C-matrix
- Often related to conventional basis normalization

---

### **CV[itc][inc][incv]** — Correlation Function Value

**From:** Runtime (CVM state variable)
**Computed:** From current occupation variables `{σᵢ}`
**Role:** Actual CF value at this thermodynamic state

- Three-level index: `[cluster type][cluster within type][CF index]`
- Example: `CV[1][1][1]` might be the point correlation (mean composition)
- Example: `CV[1][2][1]` might be a pair correlation
- Range: typically [-1, 1] in conventional basis

**Computation path:**
```
Occupation variables {σᵢ, σᵢσⱼ, ...}
  ↓
Site operator rules {pRules}
  ↓
Substitution rules {substituteRules}
  ↓
C-matrix projection {cmat · u}
  ↓
Correlation values CV[itc][inc][incv]
```

---

### **lc[[itc]]** — Cluster Count Per Class

**From:** Stage 2 (Given-Phase Clusters)
**Source class:** `ClusterIdentificationResult.orderedClusterResult.ordClusterData.orbitList` (length)
**Role:** Loop limit for middle sum

- One value per HSP cluster type
- Number of distinct clusters in that class
- Example: `lc = {1, 2, 3, ...}` means class 1 has 1 cluster, class 2 has 2, etc.

---

### **tcdis** — Total HSP Cluster Types

**From:** Stage 1 (HSP Clusters)
**Source class:** `ClusterIdentificationResult.disorderedClusterResult.disClusterData`
**Role:** Loop limit for outer sum

- Scalar: total distinct cluster types in HSP
- Example: `tcdis = 5` means summing over 5 different cluster types

---

## Data Integration Diagram

```
Type-1a generates and stores:
│
├─ Stage 1: HSP Clusters
│  └─ kbdis[itc]              (KB coefficients)
│  └─ mhdis[itc]              (HSP multiplicities)
│  └─ nijTable                (for recomputing KB if needed)
│
├─ Stage 2: Given-Phase Clusters
│  └─ mh[itc][inc]            (multiplicity ratios)
│  └─ lc[itc]                 (cluster count per class)
│
├─ Stage 4: Given-Phase CFs (grouped by HSP)
│  └─ [lcf organized by HSP cluster structure]
│
└─ C-Matrix Stage
   └─ lcv[itc][inc]           (CF count per cluster)
   └─ wcv[itc][inc][incv]     (CF weights)

Type-2 (CVM) uses at runtime:
│
├─ Load cluster_data.json
│
├─ Initialize:
│  ├─ kbdis, mhdis, mh, lc, lcv, wcv from disk
│  └─ tcdis = length(kbdis)
│
├─ Solve Lagrange multipliers to find occupation variables
│
├─ At each iteration:
│  ├─ Compute CV[itc][inc][incv] from occupation variables
│  └─ Evaluate entropy: Sc = -R * Sum[KB × m × mh × w × CV × Log[CV]]
│
└─ Extract thermodynamic properties
```

---

## Why This Organization Matters

### Cluster Type Hierarchy (itc)

All data is organized by **HSP cluster type** (`itc`). This means:

- Clusters of the **same symmetry class** contribute together to entropy
- KB coefficients are defined per symmetry class
- CF weights are consistent with cluster structure

### Cluster Within Class (inc)

Within each symmetry class, we iterate over **individual clusters** (`inc`). This allows:

- Different multiplicities for different cluster orientations (`mh[itc][inc]`)
- Different CF patterns for different cluster positions
- Correct averaging over orbit members

### CF Organization (incv)

CFs are indexed relative to their **parent cluster** (`incv`). This ensures:

- CFs are grouped by the cluster they describe
- Weights `wcv` match the CF ordering
- Loop limits `lcv` are well-defined

---

## The Final Substitution

```mathematica
Sc = Sc /. Table[msdis[[i]] -> mhdis[[i]], {i, tcdis}]
```

**Why replace `msdis` with `mhdis`?**

- `msdis` (disordered multiplicities) was used during derivation for clarity
- `mhdis` (HSP multiplicities) are the actual invariants used in CVM
- The KB coefficients were derived using HSP multiplicities
- Substitution ensures thermodynamic consistency

**Result:** Entropy correctly weighted by HSP multiplicity structure.

---

## Numerical Example (Hypothetical)

Given:
- `tcdis = 2` (two HSP cluster types)
- `lc = {1, 2}` (class 1: 1 cluster; class 2: 2 clusters)
- `lcv = {{3}, {2, 1}}` (class 1: 3 CFs; class 2: 2 CFs in first cluster, 1 in second)

Entropy calculation expands to:

```
Sc = -R * (
  kbdis[1] * mhdis[1] * (
    mh[1][1] * (wcv[1][1][1]*CV[1][1][1]*Log[...] + wcv[1][1][2]*CV[1][1][2]*Log[...] + wcv[1][1][3]*CV[1][1][3]*Log[...])
  )
  +
  kbdis[2] * mhdis[2] * (
    mh[2][1] * (wcv[2][1][1]*CV[2][1][1]*Log[...] + wcv[2][1][2]*CV[2][1][2]*Log[...])
    +
    mh[2][2] * (wcv[2][2][1]*CV[2][2][1]*Log[...])
  )
)
```

Each term scales by:
1. KB coefficient for the cluster type
2. HSP multiplicity
3. Relative multiplicity ratio
4. Weighted correlation values

---

## Verification: Data Consistency

When running CVM, ensure:

- [ ] `len(kbdis)` == `tcdis` — one KB per cluster type
- [ ] `len(mhdis)` == `tcdis` — one multiplicity per cluster type
- [ ] `len(mh)` == `tcdis` — all cluster types have multiplicity data
- [ ] `len(mh[itc])` == `lc[itc]` — multiplicity data matches cluster count
- [ ] `len(lcv[itc])` == `lc[itc]` — CF count data matches cluster count
- [ ] `len(wcv[itc][inc])` == `lcv[itc][inc]` — weight data matches CF count
- [ ] All `CV[itc][inc][incv]` values are in valid range (typically [-1, 1])

---

## References

- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — source of kbdis, mhdis, mh, lc, lcv, wcv
- [CLUSTER_DATA_GENERATION_FLOW.md](CLUSTER_DATA_GENERATION_FLOW.md) — how these values are computed
- [CMATRIX_CALCULATION_FLOW.md](CMATRIX_CALCULATION_FLOW.md) — source of lcv and wcv
- [SIMPLIFIED_ARCHITECTURE.md](SIMPLIFIED_ARCHITECTURE.md) — CVM engine as domain layer component
