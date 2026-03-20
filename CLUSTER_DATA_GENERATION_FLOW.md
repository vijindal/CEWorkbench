# Cluster Data Generation Flow: Mathematica Implementation

## The Four Stages (with actual code)

### Stage 1: HSP (Highest Symmetric Phase) Clusters

```mathematica
basisSymbolListBin = genBasisSymbolList[numCompBin, siteOpSymbol];
disClusData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolListBin];

disClusCoordList = disClusData[[1]];
mhdis = disClusData[[2]];           (* Multiplicities *)
disClusOrbitList = disClusData[[3]];
rcdis = disClusData[[4]];           (* Sublattice site counts *)
tcdis = disClusData[[5]] - 1;       (* Total cluster types *)

nijTable = getNijTable[disClusCoordList, mhdis, disClusOrbitList];
kbdis = generateKikuchiBakerCoefficients[msdis, nijTable];
```

**Inputs:**
- `disMaxClusCoord` — coordinate list of maximal clusters in HSP
- `disSymOpList` — symmetry operations of HSP
- `basisSymbolListBin` — site occupation operators (binary: only s[1])

**Outputs:**
- `mhdis` — HSP cluster multiplicities
- `nijTable` — Kikuchi-Baker nij table (invariant across all phases)
- `kbdis` — Kikuchi-Baker coefficients for HSP

**Purpose:** Foundation layer. `nijTable` and `kbdis` are derived once from the highest-symmetry phase and remain valid for all lower phases.

---

### Stage 2: Given Phase Clusters

```mathematica
(* Generate clusters in the given phase *)
clusData = genClusCoordList[maxClusCoord, symOpList, basisSymbolListBin];
clusCoordList = clusData[[1]];
clusMList = clusData[[2]];
clusOrbitList = clusData[[3]];

(* Transform to HSP reference frame *)
clusCoordList = ordToDisordCoord[rotateMat, translateMat, clusCoordList];

(* Classify given-phase clusters as members of HSP clusters *)
ordClusData = transClusCoordList[disClusData, clusData, clusCoordList];
ordClusCoordList = ordClusData[[1]];
ordClusMList = ordClusData[[2]];
ordClusOrbitList = ordClusData[[3]];
ordClusRcList = ordClusData[[4]];

lc = Map[Length, ordClusCoordList];  (* Length of each cluster class *)
tc = clusData[[5]] - 1;              (* Total cluster types in given phase *)
nxc = lc[[tcdis]];                   (* Point clusters *)
nc = tc - nxc;                        (* Pair+ clusters *)
rc = ordClusRcList;
mh = ordClusMList / mhdis;            (* Ratio: given-phase / HSP multiplicities *)
```

**Inputs:**
- `maxClusCoord` — coordinate list of maximal clusters in given phase
- `symOpList` — symmetry operations of given phase
- `rotateMat`, `translateMat` — transformation from given phase to HSP frame

**Outputs:**
- `lc` — list of cluster class lengths (grouped by HSP clusters)
- `tc` — total cluster types in given phase
- `nxc`, `nc` — point vs. pair+ cluster counts
- `rc` — sublattice site counts
- `mh` — multiplicities relative to HSP

**Purpose:** Validation layer. Ensures given-phase clusters map correctly to HSP structure. Outputs are grouped *by HSP clusters*.

---

### Stage 3: HSP Correlation Functions

```mathematica
(* Generate CFs in HSP with full decoration *)
basisSymbolList = genBasisSymbolList[numComp, siteOpSymbol];
disCFData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolList];

disCFCoordList = disCFData[[1]];
mcfdis = disCFData[[2]];       (* HSP CF multiplicities *)
disCFOrbitList = disCFData[[3]];
rcfdis = disCFData[[4]];
tcfdis = disCFData[[5]] - 1;   (* Total CF types in HSP *)
```

**Inputs:**
- `basisSymbolList` — site operators with full decoration (multi-component)
- Same coordinate structure as Stage 1, but with expanded basis

**Outputs:**
- `tcfdis` — total CF types in HSP
- `mcfdis` — HSP CF multiplicities
- `rcfdis` — sublattice site counts for CFs

**Purpose:** Reference basis layer. Defines the CF feature space structure before applying to the given phase.

---

### Stage 4: Given Phase Correlation Functions (with grouping)

```mathematica
(* Generate CFs in the given phase *)
CFData = genClusCoordList[maxClusCoord, symOpList, basisSymbolList];
CFCoordList = CFData[[1]];

(* Transform to HSP reference frame *)
CFCoordList = ordToDisordCoord[rotateMat, translateMat, CFCoordList];

(* Classify given-phase CFs as members of HSP CFs *)
ordCFData = transClusCoordList[disCFData, CFData, CFCoordList];

(* *** KEY STEP: GROUP CFs BY HSP CLUSTERS *** *)
cfData = groupCFData[disClusData, disCFData, ordCFData, rotateMat, translateMat, basisSymbolListBin];

cfCoordList = cfData[[1]];      (* CFs grouped by HSP clusters *)
lcf = readLength[cfCoordList];  (* CF count per group *)
tcf = Sum[Sum[lcf[[i]][[j]], {j, 1, Length[lcf[[i]]]}], {i, 1, tcdis}];  (* Total CFs *)
nxcf = Sum[lcf[[tcdis]][[j]], {j, 1, Length[lcf[[tcdis]]]}];             (* Point CFs *)
ncf = tcf - nxcf;               (* Pair+ CFs *)
mcf = cfData[[2]];              (* CF multiplicities *)
cfOrbitList = cfData[[3]];
rcf = cfData[[4]];              (* Sublattice site counts *)
```

**Inputs:**
- `disClusData` — HSP clusters (Stage 1)
- `disCFData` — HSP CFs (Stage 3)
- `ordCFData` — given-phase CFs (before grouping)

**Outputs:**
- `cfCoordList` — CFs grouped by HSP cluster correspondence
- `lcf` — CF count per group
- `tcf` — total CFs in given phase
- `nxcf`, `ncf` — point vs. pair+ CF counts
- `mcf` — CF multiplicities
- `rcf` — sublattice site counts

**Purpose:** Final layer. **This is where all four stages come together.** The `groupCFData` function takes:
1. HSP clusters (Stage 1)
2. HSP CFs (Stage 3)
3. Given-phase CFs (from Stage 4 logic)

...and outputs CFs organized by HSP cluster correspondence.

---

## Why All Four Stages Are Essential

| Stage | Function | Required by Next Stage | Runtime Use |
|-------|----------|------------------------|------------|
| 1: HSP Clusters | Foundation: multiplicities, `nijTable`, KB coeff | Stage 2 (validation), Stage 4 (grouping) | Multiplicity, KB |
| 2: Given-phase Clusters | Validation + structure mapping | Stage 4 (grouping reference) | *(consumed by Stage 4)* |
| 3: HSP CFs | Reference basis with full decoration | Stage 4 (grouping reference) | *(input to grouping)* |
| 4: Given-phase CFs | **Final feature space** | *(final output)* | **CF scalars: lcf, tcf, ncf** |

### Critical: The `groupCFData` Call

```mathematica
cfData = groupCFData[disClusData, disCFData, ordCFData, ...]
```

This single function call depends on **all three previous stages**:
- **`disClusData`** (Stage 1) — reference cluster structure
- **`disCFData`** (Stage 3) — reference CF basis
- **`ordCFData`** (Stage 4 intermediate) — given-phase CFs to organize

Without completing Stages 1–3 correctly, Stage 4's output (`cfCoordList` grouped by HSP clusters) is incorrect.

---

## Data Size Implications

From `cluster_data.json` (345 KB baseline):

| Stage | Output Stored | Size | Runtime Read | Can Remove? |
|-------|---|---|---|---|
| 1 | `mhdis`, `nijTable`, `kbdis` | ~70 KB | ✓ Multiplicity, KB coeff | No (needed at runtime) |
| 2 | `ordClusCoordList`, `ordClusMList`, etc. | ~73 KB | ✗ Never | **Yes** — transient |
| 3 | `disCFData`, `mcfdis`, `rcfdis` | ~98 KB | ✗ Never | **Yes** — transient |
| 4 | `cfCoordList`, `lcf`, `tcf`, `mcf` | ~98 KB | ✓ All CF scalars | No (needed at runtime) |

**Optimization:** Remove Stages 2–3 from serialization (but still generate them during workflow).

Result: 345 KB → ~13.6 KB (**96% reduction**)

---

## Implementation in Java

The Mathematica flow maps to Java as:

```
ClusterIdentificationWorkflow
├─ Stage 1: genClusterData(hsp)
│  └─ Output: ClusterIdentificationResult.getDisorderedClusterResult()
│     (contains mhdis, nijTable, kbCoeff)
├─ Stage 2: genClusterData(givenPhase)
│  └─ Output: ClusterIdentificationResult.getOrderedClusterResult()
│     [REMOVED FROM SERIALIZATION via @JsonIgnore]
├─ Stage 3: genCFData(hsp, numComp)
│  └─ Output: CFIdentificationResult.getDisorderedCFResult()
│     [REMOVED FROM SERIALIZATION via @JsonIgnore]
└─ Stage 4: genCFData(givenPhase, numComp) + groupCFData(...)
   └─ Output: CFIdentificationResult.getOrderedCFResult()
      (contains cfCoordList, lcf, tcf, ncf, mcf)
      [KEPT IN SERIALIZATION — runtime feature space]
```

---

## Verification Checklist

When implementing the storage optimization:

- [ ] Stage 1 output (`mhdis`, `nijTable`, `kbdis`) remains in serialized data
- [ ] Stage 2 output marked with `@JsonIgnore` on `getOrderedClusterResult()`
- [ ] Stage 3 output marked with `@JsonIgnore` on `getDisorderedCFResult()` OR `getPhaseCFData()` etc.
- [ ] Stage 4 output (`cfCoordList`, `lcf`, `tcf`, `mcf`) remains in serialized data
- [ ] Backward compatibility: `@JsonIgnoreProperties(ignoreUnknown = true)` on affected classes
- [ ] Test: load old JSON with all four stages → deserializes correctly
- [ ] Test: run CVM/MCS on new optimized JSON → same thermodynamic results

---

## References

- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — conceptual layers
- [CLUSTER_DATA_STORAGE_OPTIMIZATION.md](CLUSTER_DATA_STORAGE_OPTIMIZATION.md) — Java implementation details
- [README_CLUSTER_DATA.md](README_CLUSTER_DATA.md) — overview guide
