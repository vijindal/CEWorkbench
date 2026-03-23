# Cluster Data Reference: Generation → Thermodynamics

## Architecture Overview

Cluster data is produced by the **Type-1a (ClusterIdentificationWorkflow)** in four stages, then consumed by the **Type-2 (CVM/MCS)** thermodynamic engines.

```
Type-1a                                   Type-2
─────────────────────────────────         ──────────────────────────
Stage 1: HSP Clusters → kbdis, mhdis  ──► CVM entropy: kbdis × mhdis
Stage 2: Phase Clusters → mh, lc      ──► CVM entropy: mh, lc
Stage 3: HSP CFs (transient)          ─── Not stored (grouping input only)
Stage 4: Phase CFs → cfCoordList,lcf  ──► Feature space (ncf, lcf)
C-Matrix: cmat, lcv, wcv, cv          ──► CVM/MCS projection at runtime
```

---

## Stage 1 — HSP Clusters

**Mathematica:**
```mathematica
disClusData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolListBin];
mhdis = disClusData[[2]];           (* HSP multiplicities *)
tcdis = disClusData[[5]] - 1;       (* Total cluster types *)
nijTable = getNijTable[disClusCoordList, mhdis, disClusOrbitList];
kbdis = generateKikuchiBakerCoefficients[msdis, nijTable];
```

**Outputs (serialized):** `mhdis[itc]`, `kbdis[itc]`, `nijTable`

KB coefficients are invariants derived from the highest-symmetry phase; valid for all lower-symmetry phases in the system.

---

## Stage 2 — Given-Phase Clusters

**Mathematica:**
```mathematica
clusData = genClusCoordList[maxClusCoord, symOpList, basisSymbolListBin];
clusCoordList = ordToDisordCoord[rotateMat, translateMat, clusCoordList];
ordClusData = transClusCoordList[disClusData, clusData, clusCoordList];
mh = ordClusMList / mhdis;          (* Multiplicity ratio: phase / HSP *)
lc = Map[Length, ordClusCoordList]; (* Clusters per HSP class *)
```

**Outputs (serialized as scalars):** `mh[itc][inc]`, `lc[itc]`
Raw `ordClusCoordList` etc. are **not serialized** (`@JsonIgnore`).

---

## Stage 3 — HSP Correlation Functions *(transient)*

```mathematica
disCFData = genClusCoordList[disMaxClusCoord, disSymOpList, basisSymbolList];
tcfdis = disCFData[[5]] - 1;
```

Used only as reference input to Stage 4 grouping. **Not serialized.**

---

## Stage 4 — Given-Phase CFs (grouped by HSP cluster)

```mathematica
CFData = genClusCoordList[maxClusCoord, symOpList, basisSymbolList];
CFCoordList = ordToDisordCoord[rotateMat, translateMat, CFCoordList];
ordCFData = transClusCoordList[disCFData, CFData, CFCoordList];
cfData = groupCFData[disClusData, disCFData, ordCFData, rotateMat, translateMat, basisSymbolListBin];
lcf = readLength[cfCoordList];      (* CF count per cluster group *)
tcf = Sum[Sum[lcf[[i]][[j]], ...]]  (* Total CFs *)
```

**Outputs (serialized):** `cfCoordList`, `lcf[itc][inc]`, `tcf`, `mcf[itc][inc][incv]`

The `groupCFData` call depends on all three previous stages.

---

## C-Matrix Stage

```mathematica
maxClusSiteList = genSiteList[maxClusCoord];
groupCfCoordList = groupSubClus[maxClusCoord, cfData, basisSymbolList];
cfSiteOpList    = genCfSiteOpList[groupCfCoordList, maxClusSiteList];
substituteRules = genSubstituteRules[cfSiteOpList, cfSymbol];
pRules          = genPRules[numSites, numComp, siteOcSymbol, siteOpSymbol];
cMatData = genCV[maxClusSiteList, ordClusData, cfData, substituteRules, pRules, ...];
cmat = cMatData[[1]];   lcv = cMatData[[2]];   wcv = cMatData[[3]];
cv = Table[Table[cmat[[i]][[j]].uList, {j,1,lc[[i]]}], {i,1,tcdis}];
```

**Outputs (always serialized in full):** `cmat`, `lcv[itc][inc]`, `wcv[itc][inc][incv]`, `cv[itc][inc]`

Recomputation is exponentially expensive; storage cost is negligible by comparison.

---

## Serialized `cluster_data.json` Structure

```
disorderedClusterResult
  ├── mhdis[itc]          Stage 1 ✓
  ├── kbCoefficients[itc] Stage 1 ✓
  ├── nijTable            Stage 1 ✓
  ├── mh[itc][inc]        Stage 2 scalars ✓
  └── lc[itc]             Stage 2 scalars ✓

orderedCFResult
  ├── cfCoordList         Stage 4 ✓
  ├── lcf[itc][inc]       Stage 4 ✓
  ├── tcf, ncf, mcf       Stage 4 ✓
  └── rcf[itc][inc]       Stage 4 ✓

cmatrixResult
  ├── cmat                C-Matrix ✓
  ├── lcv[itc][inc]       C-Matrix ✓
  ├── wcv[itc][inc][incv] C-Matrix ✓
  └── cfBasisIndices      C-Matrix ✓  (CF site-operator decoration indices)

(Stages 2 raw coords, Stage 3: removed via @JsonIgnore → ~96% file size reduction)

**Note:** `cv[itc][inc]` is **NOT** a stored field. It is computed at runtime as the dot product:
`cv[itc][inc] = cmat[itc][inc] · uList` during each CVM N-R iteration.
```

---

## Type-2 Runtime: CVM Thermodynamics

### Configurational Entropy

```mathematica
Sc = -R * Sum[
  kbdis[[itc]] * mhdis[[itc]] *
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

**CV[itc][inc][incv]** — correlation values computed at runtime from occupation variables via C-matrix projection: `CV = cmat · occupations`.

### Configurational Enthalpy

```mathematica
Hc = Sum[mhdis[[i]] * eList[[i]] * (uList[[i]] - (1-xB)*uA[[i]] - xB*uB[[i]]), {i, tcf}];
Hc = Hc /. {u[tcdis][1][1] -> (2*xB - 1)} /. {u[tcdis+1][1][1] -> 1};
```

- `eList[i]` — ECI from CE Hamiltonian (external fit)
- `uList[i]` — CF values from CVM solver
- `uA[i]`, `uB[i]` — pure-phase reference CFs (linear interpolation baseline)
- Point correlation substitution: `u = 2xB - 1` enforces composition constraint

### Gibbs Free Energy

```
Gc = Hc - T · Sc
```

---

## Variable Index

| Symbol | Source | Index | Role |
|--------|--------|-------|------|
| `tcdis` | Stage 1 | scalar | outer loop limit |
| `kbdis[itc]` | Stage 1 | per cluster type | KB coefficient |
| `mhdis[itc]` | Stage 1 | per cluster type | HSP multiplicity |
| `lc[itc]` | Stage 2 | per cluster type | cluster count per class |
| `mh[itc][inc]` | Stage 2 | per cluster | multiplicity ratio (phase/HSP) |
| `lcv[itc][inc]` | C-Matrix | per cluster | CF count per cluster |
| `wcv[itc][inc][incv]` | C-Matrix | per CF | CF weight/normalization |
| `cfBasisIndices[l]` | C-Matrix | per CF | Site-operator basis indices `[k₁,k₂,...,kₙ]`; routes point CFs & builds N-R initial guess |
| `cmat` | C-Matrix | matrix | occupation→CF projection |
| `CV[itc][inc][incv]` | Runtime | per CF | CF value (from cmat·occ) |
| `tcf` | Stage 4 | scalar | total CF count |
| `lcf[itc][inc]` | Stage 4 | per cluster | CF count per cluster group |
| `eList[i]` | Hamiltonian | per CF | ECI (fitted) |
| `uList[i]` | Runtime | per CF | CF value (flat) |
| `uA[i]`, `uB[i]` | Runtime | per CF | pure-phase reference CFs |

---

## Java Class Mapping

```
ClusterIdentificationWorkflow
├── Stage 1 → disorderedClusterResult  (AllClusterData.getDisorderedClusterResult())
├── Stage 2 → [transient ordClusData]  (@JsonIgnore on getOrderedClusterResult())
├── Stage 3 → [transient disCFData]    (@JsonIgnore on getDisorderedCFResult())
├── Stage 4 → orderedCFResult          (AllClusterData.getDisorderedCFResult().getNcf(), etc.)
└── C-Matrix→ cmatrixResult            (AllClusterData.getCMatrixResult())

CVMEngine.compute(ThermodynamicInput)
├── cvmInput = new CVMInput(disorderedClusterResult, disorderedCFResult, cmatrixResult, ...)
├── eci[] = evaluateECI(cec, T)        (eci[l] = a[l] + b[l]*T)
└── CVMPhaseModel.create(cvmInput, eci, T, composition) → EquilibriumState
    ├── cfBasisIndices → ClusterVariableEvaluator.computeRandomCFs()
    │                    (N-R initial guess: uRand[l] = Π pointCFs[kᵢ-1])
    ├── cfBasisIndices → ClusterVariableEvaluator.buildFullCFVector()
    │                    (routes point CFs to correct columns for ternary K≥3)
    └── cfBasisIndices → NewtonRaphsonSolverSimple.CVMData.cfRank[]
                         (cluster rank = cfBasisIndices[l].length)

MCSEngine.compute(ThermodynamicInput)
├── clusterData = disorderedClusterResult.getDisClusterData()
├── eci[] = extractEci(cec.cecTerms, T)
└── MCSRunner.build().run() → MCResult → EquilibriumState
```
