# Complete Cluster Data Flow: Type-1a Generation → Type-2 Thermodynamics

## The Full Picture

This document shows **how cluster data flows** from generation (Type-1a) to usage (Type-2), and **why** every data element is structured the way it is.

---

## Type-1a: Cluster Identification Workflow

### Input
- Crystal structure (given phase)
- Symmetry operations
- Number of components
- Site occupation operators

### Output: `cluster_data.json`

---

## Stage 1: HSP Clusters

```mathematica
disClusData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolListBin];
nijTable = getNijTable[disClusCoordList, mhdis, disClusOrbitList];
kbdis = generateKikuchiBakerCoefficients[msdis, nijTable];
```

**Outputs stored:**
- `mhdis[itc]` — HSP multiplicities per cluster type
- `kbdis[itc]` — Kikuchi-Baker coefficients per cluster type
- `nijTable` — Korringa-Kohn-Rostoker table (rarely accessed at runtime)

**Why necessary:**
- KB coefficients are **invariants** derived once from highest-symmetry phase
- Valid for all phases in the system
- Required for CVM entropy calculation

**Runtime usage:** ✓ Used by CVM/MCS

---

## Stage 2: Given-Phase Clusters

```mathematica
clusData = genClusCoordList[maxClusCoord, symOpList, basisSymbolListBin];
clusCoordList = ordToDisordCoord[rotateMat, translateMat, clusCoordList];
ordClusData = transClusCoordList[disClusData, clusData, clusCoordList];
mh = ordClusMList / mhdis;
lc = Map[Length, ordClusCoordList];
```

**Outputs:**
- `mh[itc][inc]` — multiplicity ratios (given-phase / HSP)
- `lc[itc]` — cluster count per HSP class
- `ordClusCoordList`, `ordClusOrbitList` — detailed cluster structures

**Stored in serialized form?** ✗ No (optimization removes via `@JsonIgnore`)

**Why generated but not stored:**
- Needed to compute `mh` and `lc`
- Validates cluster-CF correspondence
- Not needed at runtime (CVM uses only `mh` and `lc` scalars)

**Runtime usage:** ✗ **Not used directly** (scalars extracted)

---

## Stage 3: HSP Correlation Functions

```mathematica
disCFData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolList];
tcfdis = disCFData[[5]] - 1;
mcfdis = disCFData[[2]];
```

**Outputs:**
- `tcfdis` — total CF types in HSP
- `mcfdis` — HSP CF multiplicities
- `disCFCoordList` — HSP CF coordinates

**Stored in serialized form?** ✗ No (optimization removes)

**Why generated but not stored:**
- Reference basis for grouping lower-phase CFs
- Defines feature space structure
- Not accessed at runtime (grouping is pre-computed)

**Runtime usage:** ✗ **Not used** (only as input to Stage 4)

---

## Stage 4: Given-Phase Correlation Functions (with grouping)

```mathematica
CFData = genClusCoordList[maxClusCoord, symOpList, basisSymbolList];
CFCoordList = ordToDisordCoord[rotateMat, translateMat, CFCoordList];
ordCFData = transClusCoordList[disCFData, CFData, CFCoordList];
cfData = groupCFData[disClusData, disCFData, ordCFData, ...];
lcf = readLength[cfCoordList];
tcf = Sum[...];
mcf = cfData[[2]];
```

**Outputs:**
- `cfCoordList` — CFs grouped by HSP cluster correspondence
- `lcf[itc][inc]` — CF count per cluster
- `tcf` — total CFs
- `mcf[itc][inc][incv]` — CF multiplicities

**Stored in serialized form?** ✓ Yes (critical runtime data)

**Why stored:**
- Defines the actual feature space for thermodynamic calculations
- Grouped by HSP clusters (organizing principle)
- Required for CVM/MCS at runtime

**Runtime usage:** ✓ **Used by CVM/MCS**

---

## C-Matrix Stage

```mathematica
maxClusSiteList = genSiteList[maxClusCoord];
cfSiteOpList = genCfSiteOpList[groupCfCoordList, maxClusSiteList];
substituteRules = genSubstituteRules[cfSiteOpList, cfSymbol];
pRules = genPRules[numSites, numComp, siteOcSymbol, siteOpSymbol];
cMatData = genCV[maxClusSiteList, ordClusData, cfData, ...];
cmat = cMatData[[1]];
lcv = cMatData[[2]];
wcv = cMatData[[3]];
cv = Table[Table[cmat[[i]][[j]].uList, ...], ...];
```

**Inputs (dependencies):**
- `ordClusData` from Stage 2
- `cfData` from Stage 4
- Cluster geometry and CF organization

**Outputs:**
- `cmat[itc][inc][incv]` — correlation matrix rows
- `lcv[itc][inc]` — CF count per cluster (redundant with Stage 4 but organized)
- `wcv[itc][inc][incv]` — CF weights
- `cv[itc][inc]` — organized correlation vectors

**Stored in serialized form?** ✓ Yes (fully retained)

**Why stored:**
- Exponentially expensive to recompute
- Required for CVM/MCS at runtime
- Storage cost negligible compared to generation cost

**Runtime usage:** ✓ **Used by CVM/MCS**

---

## Summary: What Gets Serialized

```
cluster_data.json (optimized: ~13.6 KB from 345 KB)

├─ disorderedClusterResult
│  ├─ mhdis[itc]                     ✓ (Stage 1)
│  ├─ nijTable                       ✓ (Stage 1)
│  ├─ kbCoefficients                 ✓ (Stage 1)
│  ├─ ordClusterData (raw coords)    ✗ (Stage 2 — removed)
│  ├─ ordClusMList → mh[itc][inc]    ✓ (Stage 2 — kept as scalar)
│  └─ orbitList → lc[itc]            ✓ (Stage 2 — kept as scalar)
│
├─ disorderedCFResult
│  ├─ tcfdis, mcfdis, rcfdis         ✗ (Stage 3 — removed)
│  ├─ disCFData (raw coords)         ✗ (Stage 3 — removed)
│  └─ groupedCFData                  ✗ (Stage 4 intermediate — removed)
│
├─ orderedCFResult (grouped)
│  ├─ cfCoordList                    ✓ (Stage 4 — kept)
│  ├─ lcf[itc][inc]                  ✓ (Stage 4 — kept)
│  ├─ tcf, ncf, mcf                  ✓ (Stage 4 — kept)
│  └─ rcf[itc][inc]                  ✓ (Stage 4 — kept)
│
└─ cmatrixResult
   ├─ cmat                           ✓ (C-Matrix stage)
   ├─ lcv[itc][inc]                  ✓ (C-Matrix stage)
   ├─ wcv[itc][inc][incv]            ✓ (C-Matrix stage)
   └─ cv[itc][inc]                   ✓ (C-Matrix stage)
```

---

## Type-2: CVM Thermodynamic Calculation

### Loading Data

```java
ClusterDataStore.load("cluster_data.json")
  → ClusterData object containing:
      - kbdis, mhdis, mh, lc
      - cfCoordList, lcf, tcf, mcf, rcf
      - cmat, lcv, wcv, cv
```

### CVM Entropy Calculation

At each thermodynamic state point:

```mathematica
Sc = -R * Sum[
  kbdis[[itc]] *           (* From Stage 1 *)
  mhdis[[itc]] *           (* From Stage 1 *)
  Sum[
    mh[[itc]][[inc]] *     (* From Stage 2 *)
    Sum[
      wcv[[itc]][[inc]][[incv]] * CV[itc][inc][incv] * Log[CV[itc][inc][incv]],
      (* wcv from C-Matrix, CV computed at runtime *)
      {incv, 1, lcv[[itc]][[inc]]}  (* lcv from C-Matrix *)
    ],
    {inc, 1, lc[[itc]]}  (* lc from Stage 2 *)
  ],
  {itc, 1, tcdis}  (* tcdis = length(kbdis) *)
];
```

**Data flow in entropy calculation:**

```
Occupation variables {σᵢ, σᵢσⱼ, ...}
  ↓ [pRules from Stage 4]
Site operators {s₁, s₂, ...}
  ↓ [substituteRules from C-Matrix]
CF contributions {CF[1,1,1], CF[1,1,2], ...}
  ↓ [cmat from C-Matrix]
Correlation values CV[itc][inc][incv]
  ↓ [multiply by wcv[[itc]][[inc]][[incv]]]
Weighted contributions
  ↓ [sum over incv, then inc, then itc]
  ↓ [scale by kbdis, mhdis, mh]
Entropy Sc
```

---

## Why The Four-Layer Architecture

Now we can see **exactly why** all four stages are necessary:

| Stage | Generates | Used By | Purpose |
|-------|-----------|---------|---------|
| 1 | `kbdis`, `nijTable` | Stage 2,4,C; CVM entropy | Thermodynamic foundation |
| 2 | `mh`, `lc` | Stage 4,C; CVM entropy | Cluster structure mapping |
| 3 | Intermediate CF basis | Stage 4 grouping | Reference organization |
| 4 | `cfCoordList`, `lcf`, `mcf` | C-Matrix,CVM entropy | Feature space definition |
| C | `cmat`, `lcv`, `wcv`, `cv` | CVM entropy at runtime | Correlation projection |

**Each stage is a prerequisite for the next:**
- Stage 1 → Invariant KB coefficients
- Stages 1,2 → Cluster mapping
- Stages 1,2,3 → CF grouping in Stage 4
- Stages 1,2,4 → C-matrix computation
- All stages → CVM entropy calculation

---

## Optimization Justification

**Why remove Stages 2 and 3 from serialization?**

```
Stage 2 outputs: ordClusCoordList, ordClusOrbitList, ordClusRcList
  → Only `mh` and `lc` scalars used at runtime ✓ Kept
  → Raw coordinates never accessed ✗ Removed

Stage 3 outputs: disCFCoordList, disCFOrbitList, disCFData
  → Used only for grouping in Stage 4 ✓ Done at generation
  → Never accessed at runtime ✗ Removed
```

**Why keep Stages 1, 4, and C-Matrix:**

```
Stage 1: kbdis, nijTable, kbCoefficients
  → Essential for entropy calculation ✓ Keep

Stage 4: cfCoordList, lcf, tcf, mcf
  → Defines runtime feature space ✓ Keep

C-Matrix: cmat, lcv, wcv, cv
  → Exponentially expensive to recompute ✓ Keep
  → Required at every CVM iteration ✓ Keep
```

---

## Data Dependencies Graph

```
User Input
  ↓
Type-1a: Cluster Identification
  ├─ Stage 1 ─────→ kbdis ──┐
  │    ↓                     │
  ├─ Stage 2 ─────→ mh, lc ─┤
  │    ↓                     │
  ├─ Stage 3 ─────→ (discard)
  │    ↓                     │
  ├─ Stage 4 ─────→ cfData ─┤
  │    ↓                     │
  └─ C-Matrix ─────→ cmat ──┘
       ↓
  Serialization → cluster_data.json
       ↓
Type-2: CVM Thermodynamics
  ├─ Load: kbdis, mh, lc, cfData, cmat
  ├─ At each T: compute CV[itc][inc][incv]
  └─ Evaluate: Sc = -R * Sum[kbdis * mhdis * mh * wcv * CV * Log[CV]]
       ↓
  EquilibriumState (Gibbs energy, entropy, etc.)
```

---

## Complete Verification Checklist

**After Type-1a completion:**
- [ ] All four stages computed without error
- [ ] `kbdis` computed correctly from `nijTable`
- [ ] `mh = ordClusMList / mhdis` values sensible (typically 0–2)
- [ ] `lc` sums correctly: `sum(lc) = tc` (total clusters)
- [ ] `cfCoordList` properly grouped by cluster class
- [ ] `lcf` sums correctly: `sum(lcf) = tcf` (total CFs)
- [ ] `cmat` dimensions: rows = total clusters, cols = total CFs
- [ ] `lcv` matches `lcf` (CF count per cluster)
- [ ] `wcv` all positive (normalization factors)

**Before running Type-2 (CVM):**
- [ ] `cluster_data.json` loads without error
- [ ] All required fields present: `kbdis`, `mhdis`, `mh`, `lc`, `lcv`, `wcv`
- [ ] Cluster hierarchy indexing consistent: `itc < tcdis`, `inc < lc[itc]`, `incv < lcv[itc][inc]`
- [ ] No missing data for entropy calculation

**During Type-2 (CVM):**
- [ ] `CV[itc][inc][incv]` values in expected range [-1, 1]
- [ ] Entropy computation includes all terms from all cluster types
- [ ] KB coefficients applied with correct signs

---

## References

- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — conceptual architecture
- [CLUSTER_DATA_GENERATION_FLOW.md](CLUSTER_DATA_GENERATION_FLOW.md) — Stage 1-4 implementation
- [CMATRIX_CALCULATION_FLOW.md](CMATRIX_CALCULATION_FLOW.md) — C-matrix computation
- [CVM_ENTROPY_CALCULATION.md](CVM_ENTROPY_CALCULATION.md) — Type-2 usage in entropy
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — serialization choices
- [README_CLUSTER_DATA.md](README_CLUSTER_DATA.md) — overview guide
