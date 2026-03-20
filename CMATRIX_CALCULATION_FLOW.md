# C-Matrix Calculation Flow

## Overview

After the four-layer cluster and correlation function generation, the **C-matrix (correlation matrix)** is computed. This matrix is exponentially expensive to compute and is therefore **always retained in serialized form** (never removed for optimization).

The C-matrix represents the mapping from correlation functions (CFs) to cluster occupation probabilities in the given phase, organized by HSP cluster correspondence.

---

## The C-Matrix Calculation (Mathematica)

### Step 1: Generate Site List for Maximal Cluster

```mathematica
maxClusSiteList = genSiteList[maxClusCoord];
(* Site coordinates of maximal cluster are assigned a site-number *)
(* Example: {site1, site2, ..., siteN} with indices *)
```

**Purpose:** Map each site in the maximal cluster to an index for site operator assignment.

---

### Step 2: Generate CF Site Operators

```mathematica
groupCfCoordList = groupSubClus[maxClusCoord, cfData, basisSymbolList];
(* Generating all possible CFs of maxClusCoord with basisSymbolList and grouping *)

cfSiteOpList = genCfSiteOpList[groupCfCoordList, maxClusSiteList];
(* List of site operators of each CF *)
```

**Purpose:**
- `groupCfCoordList` — all possible CFs within the maximal cluster
- `cfSiteOpList` — map each CF to its constituent site operators

**Example:**
- CF = "pair on sites 1,2" → `σ₁ⁱ σ₂ʲ` (product of site operators)

---

### Step 3: Create Substitution Rules

```mathematica
substituteRules = genSubstituteRules[cfSiteOpList, cfSymbol];
(* Maps CF symbols to site operator products *)
```

**Purpose:** Create symbolic rules like:
```
CF[1][1][1] → σ₁ⁱ σ₂ʲ
CF[1][1][2] → σ₁ⁱ σ₃ᵏ
```

These allow symbolic manipulation and later numerical substitution.

---

### Step 4: Generate Occupation Number Rules

```mathematica
numSites = Length[maxClusSiteList];
pRules = genPRules[numSites, numComp, siteOcSymbol, siteOpSymbol];
(* Conventional basis has been applied here onwards *)
```

**Purpose:**
- Create occupation probability rules for the **conventional basis**
- Maps: site operator `σᵢ` → occupation probability `P(occupation | i)`
- Example: `{σ₁ = 2P₁ - 1, σ₁² = 1, ...}`

**Why conventional basis?** Simplifies energy expressions and makes CVM approximations cleaner.

---

### Step 5: Compute C-Matrix

```mathematica
cMatData = genCV[maxClusSiteList, ordClusData, cfData, substituteRules,
                 pRules, numComp, cfSymbol, siteOcSymbol, siteOpSymbol];

wcv = cMatData[[3]];   (* Weights / scalings *)
lcv = cMatData[[2]];   (* Lengths of correlation vectors per cluster *)
cmat = cMatData[[1]];  (* The actual C-matrix *)
```

**Inputs (dependencies):**
- `maxClusSiteList` — site indices in maximal cluster
- `ordClusData` — given-phase clusters (Stage 2)
- `cfData` — given-phase CFs grouped by HSP clusters (Stage 4)
- `substituteRules` — CF → site operator mappings
- `pRules` — occupation number rules

**Outputs:**
- `cmat` — matrix rows correspond to HSP cluster classes; columns correspond to CFs
- `lcv` — length of correlation vector for each cluster class
- `wcv` — weights or normalization factors

**Computational cost:** Exponential in number of sites and components. This is why it's kept in full (never discarded).

---

### Step 6: Create Correlation Function Symbols

```mathematica
uuList = Table[Table[Table[cfSymbol[i][j][k], {k, 1, lcf[[i]][[j]]}],
                     {j, 1, Length[lcf[[i]]]}], {i, 1, Length[lcf]}];
uList = Flatten[uuList];
(* uList = {u[1,1,1], u[1,1,2], ..., u[i,j,k], ...} *)
(* where lcf[[i]][[j]] is the count of CFs in group (i,j) *)
```

**Purpose:** Create symbolic CF variables `u[i,j,k]` indexed by:
- `i` — HSP cluster class
- `j` — cluster within class
- `k` — CF within cluster

---

### Step 7: Organize C-Matrix by Cluster Class

```mathematica
cv = Table[Table[cmat[[i]][[j]].uList, {j, 1, lc[[i]]}], {i, 1, tcdis}];
(* cv[[i]][[j]] = C-matrix row for cluster j of class i
   dotted with correlation function vector uList *)
```

**Purpose:**
- Group C-matrix rows by HSP cluster class
- Each `cv[[i]][[j]]` is a linear combination of CFs
- Represents the "average correlation" for that cluster in that class

**Output structure:**
```
cv = {
  { cv[1,1], cv[1,2], ..., cv[1,lc[1]] },     (* Class 1 *)
  { cv[2,1], cv[2,2], ..., cv[2,lc[2]] },     (* Class 2 *)
  ...
  { cv[tcdis,1], ..., cv[tcdis,lc[tcdis]] }   (* Class tcdis *)
}
```

---

## Dependency Graph

```
Stage 1: HSP Clusters
  ↓
Stage 2: Given-phase Clusters → ordClusData, lc
  ↓
Stage 3: HSP CFs
  ↓
Stage 4: Given-phase CFs grouped → cfData, lcf, lc
  ↓
C-MATRIX CALCULATION
  ├─ Step 1: maxClusSiteList from maxClusCoord
  ├─ Step 2-3: cfSiteOpList, substituteRules from cfData
  ├─ Step 4: pRules (multi-component occupation rules)
  ├─ Step 5: genCV(ordClusData, cfData, ...) → cmat, lcv, wcv
  └─ Step 6-7: Organize by cluster class → cv
```

**Critical dependency:** C-matrix calculation requires:
- `ordClusData` from Stage 2
- `cfData` from Stage 4
- All four cluster/CF stages must complete before C-matrix can be computed

---

## What Gets Stored

From `cluster_data.json`:

```
cmatrixResult {
  cmat: [[...], [...], ...]          (* Exponentially large *)
  lcv: [...]                          (* CV vector lengths *)
  wcv: [...]                          (* Weights *)
  cv: [[...], [...], ...]             (* Organized by cluster *)
}
```

**Total size:** ~2,799 bytes in baseline (but exponentially grows with system size)

**Storage decision:** **ALWAYS KEEP FULL CMATRIX**
- Recomputation is exponentially expensive
- Storage is negligible compared to generation cost
- Required for both CVM and MCS at runtime

---

## Runtime Usage

### CVM (`CVMPhaseModel`)

At each temperature point:
1. Read `cmat` rows for each cluster class
2. Evaluate `cv = cmat · uList` where `uList` is the current CF value vector
3. Use `cv` to compute cluster free energy contributions
4. Solve Lagrange multipliers for occupation variables

### MCS (`EmbeddingGenerator`)

At each MC step:
1. Calculate proposed CF changes from trial spin flips
2. Compute energy change using C-matrix projection
3. Apply Metropolis acceptance criterion

Both engines rely on the pre-computed C-matrix; recomputation mid-calculation would be prohibitively expensive.

---

## Workflow Completeness

The Type-1a (cluster identification) workflow is **complete** after:

1. ✓ Four-layer cluster/CF generation
2. ✓ C-matrix computation
3. ✓ All data serialized to `cluster_data.json`

The Type-2 (thermodynamic) workflows then load this data and use:
- Multiplicity + KB coefficients (Stage 1)
- CF scalars (Stage 4)
- C-matrix (C-Matrix stage)

---

## Storage Optimization Notes

| Data | Size | Runtime Use | Can Remove? |
|------|------|-------------|------------|
| `cmat` | ~2,799 B | CVM/MCS | **No** — expensive |
| `lcv`, `wcv` | ~100 B | CVM/MCS | **No** — metadata |
| `cv` | ~100 B | CVM/MCS | **No** — organized |

The C-matrix section is **fully retained** in the optimized `cluster_data.json`.

---

## Verification

After C-matrix calculation, verify:

- [ ] `cmat` dimensions: rows = number of cluster classes × clusters per class
- [ ] `cmat` dimensions: columns = total number of CFs
- [ ] `lcv[i]` = sum of CFs in all clusters of class `i`
- [ ] `wcv` values are positive (weight/normalization)
- [ ] `cv` properly groups C-matrix rows by cluster class
- [ ] C-matrix correctly handles decoration per component

---

## References

- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — four-layer architecture
- [CLUSTER_DATA_GENERATION_FLOW.md](CLUSTER_DATA_GENERATION_FLOW.md) — cluster/CF generation
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — optimization plan (C-matrix kept in full)
