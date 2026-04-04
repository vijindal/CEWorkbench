# TYPE-1a: Cluster & Correlation Function Identification - Complete Dataflow

## Overview

TYPE-1a is **Stage 1-3 (Cluster Identification → CF Identification → C-Matrix Building)** of the CVM pipeline.
This is where crystal geometry becomes mathematical data structures for thermodynamic calculations.

---

## Entry Point
**File:** [Main.java:38](app/src/main/java/org/ce/ui/cli/Main.java#L38)
**Entry Method:** `main(args)` with mode="type1a"

---

## Complete Execution Flow

```
╔════════════════════════════════════════════════════════════════════════════════╗
║                           MAIN ENTRY (CLI)                                     ║
╚════════════════════════════════════════════════════════════════════════════════╝
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Main.java:38 main()                                                            │
│                                                                                │
│ Read command-line arguments:                                                   │
│  mode   = "type1a"                                                             │
│  elements = "A-B" (default)                                                    │
│  structure = "BCC_B2" (default)                                                │
│  model = "T" (default)                                                         │
│                                                                                │
│ Outputs:                                                                       │
│  CLUSTER_ID = "A-B_BCC_B2_T_bin"  (or custom)                                 │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Main.java:81-87                                                                │
│ Build ClusterIdentificationRequest                                            │
│                                                                                │
│ Configuration:                                                                 │
│  .disorderedClusterFile("clus/BCC_A2-T.txt")   ◄─ Maximal clusters (HSP)     │
│  .orderedClusterFile("clus/BCC_B2-T.txt")      ◄─ Maximal clusters (phase)   │
│  .disorderedSymmetryGroup("BCC_A2-SG")         ◄─ Symmetry ops (HSP)         │
│  .orderedSymmetryGroup("BCC_B2-SG")            ◄─ Symmetry ops (phase)       │
│  .numComponents(2)                             ◄─ K = 2 (binary)              │
│                                                                                │
│ Output: ClusterIdentificationRequest object                                   │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
╔════════════════════════════════════════════════════════════════════════════════╗
║                     STAGE 1: CLUSTER IDENTIFICATION                            ║
║                                                                                ║
║ Main.java:89                                                                   ║
║ AllClusterData clusterData = ClusterIdentificationWorkflow.identify(config);  ║
╚════════════════════════════════════════════════════════════════════════════════╝
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ ClusterIdentificationWorkflow.identify() [ClusterIdentificationWorkflow.java] │
│                                                                                │
│ This is a 5-step workflow:                                                    │
│                                                                                │
│ ══════════════════════════════════════════════════════════════════════════    │
│ STEP 1: Load Disordered Phase (HSP) Geometry & Symmetry                       │
│ ══════════════════════════════════════════════════════════════════════════    │
│                                                                                │
│ INPUT FILES:                                                                   │
│  • clus/BCC_A2-T.txt         ◄─ List of maximal clusters (atoms coordinates) │
│  • symmetry/BCC_A2-SG        ◄─ Space group operations (rotations, trans.)  │
│                                                                                │
│ PROCESSING:                                                                    │
│  InputLoader.parseClusterFile("clus/BCC_A2-T.txt")                           │
│    └─ Parses cluster coordinates from text file                              │
│                                                                                │
│  InputLoader.parseSpaceGroup("BCC_A2-SG")                                    │
│    └─ Parses space group symmetry operations                                 │
│                                                                                │
│ OUTPUT:                                                                        │
│  • List<Cluster> disorderedClusters                                          │
│  • SpaceGroup disorderedSpaceGroup                                           │
│  • List<SymmetryOperation> disorderedSymOps                                  │
│                                                                                │
│ ══════════════════════════════════════════════════════════════════════════    │
│ STEP 2: Load Ordered Phase Geometry & Symmetry                               │
│ ══════════════════════════════════════════════════════════════════════════    │
│                                                                                │
│ INPUT FILES:                                                                   │
│  • clus/BCC_B2-T.txt         ◄─ Maximal clusters (ordered phase)             │
│  • symmetry/BCC_B2-SG        ◄─ Space group operations                       │
│                                                                                │
│ PROCESSING: Same as STEP 1, but for ordered phase                            │
│                                                                                │
│ OUTPUT:                                                                        │
│  • List<Cluster> orderedClusters                                             │
│  • SpaceGroup orderedSpaceGroup                                              │
│  • List<SymmetryOperation> orderedSymOps                                     │
│  • Transformation matrix & translation vector (from space group files)        │
│                                                                                │
│ ══════════════════════════════════════════════════════════════════════════    │
│ STEP 3: Cluster Identification (Stage 1)                                      │
│ ══════════════════════════════════════════════════════════════════════════    │
│                                                                                │
│ INPUT:                                                                         │
│  • disorderedClusters, disorderedSymOps (HSP geometry)                       │
│  • orderedClusters, orderedSymOps (phase geometry)                           │
│  • Transformation matrix (ordered → HSP frame)                                │
│                                                                                │
│ CALLS:                                                                         │
│  ClusterIdentifier.identify(                                                  │
│      disMaxClusCoord, disSymOps,                                              │
│      maxClusCoord, symOps,                                                    │
│      rotateMat, translateMat)                                                 │
│                                                                                │
│ OUTPUT: ClusterIdentificationResult                                           │
│  ├─ tcdis          (number of HSP cluster types, excl. empty cluster)        │
│  ├─ tc             (total clusters including point/empty)                     │
│  ├─ mhdis[]        (HSP cluster multiplicities: mh[0], mh[1], ...)          │
│  ├─ kb[]           (Kikuchi-Baker entropy coefficients)                      │
│  ├─ mh[][]         (normalized multiplicities: mh[t][j])                     │
│  ├─ lc[]           (ordered clusters per HSP type: lc[t])                    │
│  ├─ disClusterData (HSP cluster orbits & info)                               │
│  └─ ordClusterData (ordered-phase clusters classified to HSP types)          │
│                                                                                │
│ ══════════════════════════════════════════════════════════════════════════    │
│ STEP 4: Correlation Function Identification (Stage 2)                        │
│ ══════════════════════════════════════════════════════════════════════════    │
│                                                                                │
│ INPUT:                                                                         │
│  • ClusterIdentificationResult (from STEP 3)                                 │
│  • numComponents = 2 (binary)                                                │
│  • Same geometry & symmetry data as STEP 3                                   │
│                                                                                │
│ CALLS:                                                                         │
│  CFIdentifier.identify(                                                       │
│      clusterResult,                                                           │
│      disMaxClusCoord, disSymOps,                                              │
│      maxClusCoord, symOps,                                                    │
│      rotateMat, translateMat,                                                 │
│      numComp=2)                                                               │
│                                                                                │
│ OUTPUT: CFIdentificationResult                                                │
│  ├─ tcdis, tcfdis   (HSP cluster & CF type counts)                           │
│  ├─ tcf, ncf        (total & non-point CF counts)                            │
│  ├─ nxcf            (point CF count = K-1)                                   │
│  ├─ lcf[][]         (CF count per (type, group): lcf[t][j])                  │
│  ├─ disClusterData  (HSP clusters with full decoration)                      │
│  ├─ disOrbitList    (HSP CF orbits)                                          │
│  ├─ ordClusterData  (ordered-phase CFs classified to HSP)                    │
│  └─ GroupedCFResult (CFs grouped by type and cluster)                        │
│                                                                                │
│ ══════════════════════════════════════════════════════════════════════════    │
│ STEP 5: Build C-Matrix (Stage 3)                                              │
│ ══════════════════════════════════════════════════════════════════════════    │
│                                                                                │
│ INPUT:                                                                         │
│  • ClusterIdentificationResult                                                │
│  • CFIdentificationResult                                                     │
│  • maxClusters (ordered-phase clusters as geometry reference)                │
│  • numComponents = 2                                                          │
│                                                                                │
│ CALLS:                                                                         │
│  CMatrixBuilder.build(                                                        │
│      clusterResult,                                                           │
│      cfResult,                                                                │
│      orderedClusters,  ◄─ used as reference for site enumeration            │
│      numComponents)                                                           │
│                                                                                │
│ CONSTRUCTION STEPS:                                                            │
│  1. SiteListBuilder.buildSiteList(maxClusters)                               │
│     └─ Extract all unique atomic sites from ordered clusters                 │
│     OUTPUT: List<Position> siteList (e.g., 2 sites for B2)                  │
│                                                                                │
│  2. PRulesBuilder.build(siteList.size(), numElements)                       │
│     └─ Build p-rules (substitution rules for occupations)                    │
│     OUTPUT: PRules (defines site occupancies)                                │
│                                                                                │
│  3. CFSiteOpListBuilder.build(cfSiteOpList, siteList)                       │
│     └─ Map each CF to site-operator sequences                               │
│                                                                                │
│  4. SubstituteRulesBuilder.build(cfSiteOpList, siteList)                    │
│     └─ Build polynomial substitution rules                                   │
│     OUTPUT: SubstituteRules                                                  │
│                                                                                │
│  5. Extract cfBasisIndices (basis-index decorations per CF)                   │
│     └─ Specifies which basis power each CF carries                           │
│     OUTPUT: int[][] cfBasisIndices[tcf][rank]                               │
│                                                                                │
│  6. For each ordered-phase cluster type t & group j:                         │
│     a) Generate all site configurations (C(nsites, K) choose rules)          │
│     b) For each config, compute CV row via polynomial substitution           │
│     c) Group equivalent rows (same polynomial = same CV value)               │
│     d) Assign weights (how many configs map to same CV)                      │
│     OUTPUT: cmat[t][j][v][k+1]  (C-matrix)                                  │
│             wcv[t][j][v]         (CV weights)                                │
│             lcv[t][j]            (CV count per group)                        │
│                                                                                │
│ FINAL OUTPUT: CMatrixResult                                                   │
│  ├─ cmat[tcdis][lc[t]][lcv[t][j]][tcf+1]   (the C-matrix)                  │
│  ├─ lcv[tcdis][lc[t]]                       (CV counts)                      │
│  ├─ wcv[tcdis][lc[t]][lcv[t][j]]            (CV weights)                     │
│  ├─ cfBasisIndices[tcf][]                   (CF basis decorations)           │
│  └─ other topology data                                                       │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Return AllClusterData from ClusterIdentificationWorkflow.identify()           │
│                                                                                │
│ AllClusterData contains:                                                      │
│  ├─ disClusterResult    (ClusterIdentificationResult from Stage 1)           │
│  ├─ ordClusterResult    (same as above, for symmetry in API)                │
│  ├─ disCFResult         (CFIdentificationResult from Stage 2)                │
│  ├─ ordCFResult         (same as above)                                      │
│  └─ cmatrix             (CMatrixResult from Stage 3)                         │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Main.java:108                                                                  │
│ clusterStore.save(CLUSTER_ID, fullData)                                       │
│                                                                                │
│ SAVE TO DISK:                                                                 │
│  ~/CEWorkbench/cluster-data/A-B_BCC_B2_T_bin/cluster_data.json               │
│                                                                                │
│ JSON STRUCTURE:                                                                │
│  {                                                                             │
│    "clusterId": "A-B_BCC_B2_T_bin",                                          │
│    "timestamp": "2026-03-26T...",                                             │
│    "disClusterResult": { tcdis, tc, mhdis, kb, mh, lc, ... },              │
│    "disCFResult": { tcf, ncf, nxcf, lcf, ... },                             │
│    "cmatrixResult": {                                                         │
│      "cmat": [ [ [ [c00, c01, ...], ... ], ... ], ... ],                    │
│      "lcv": [ [2, 3, ...], [1, 2, ...], ... ],                              │
│      "wcv": [ [[1, 1, 1, 2], ...], ... ],                                   │
│      "cfBasisIndices": [ [1], [1,1], [2], [1,2], ... ]                     │
│    }                                                                          │
│  }                                                                             │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────────────┐
│ Print Summary to Console                                                      │
│                                                                                │
│ Cluster Identification complete: tcdis=5, tc=6, mhdis=[12, 24, 8, 6, 1]     │
│ C-matrix built: 5 cluster types                                              │
│ Saved: ~/CEWorkbench/cluster-data/A-B_BCC_B2_T_bin/cluster_data.json         │
│                                                                                │
└────────────────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Component Breakdown

### STAGE 1: Cluster Identification

**File:** [ClusterIdentifier.java:107](app/src/main/java/org/ce/domain/cluster/ClusterIdentifier.java#L107)

#### 1a: HSP Cluster Enumeration

```
INPUTS:
  • disorderedClusters:   max clusters with HSP geometry
  • disorderedSymOps:     symmetry operations of HSP
  • basisBin = ["s1"]     binary basis (one site operator)

PROCESSING:
  ClusCoordListGenerator.generate(clusters, symOps, basisBin)
    ├─ Enumerate all clusters by site decorations
    ├─ Apply symmetry operations to find orbits
    └─ Return distinct cluster types

OUTPUTS:
  • ClusCoordListResult (disordered)
    ├─ clusList:        cluster coordinates by size
    ├─ orbitList:       clusters grouped by symmetry orbits
    └─ multiplicities:  |orbit| for each cluster type

ALGEBRA INVOLVED:
  • Cluster enumeration using generating functions
  • Symmetry orbit computation (applies rotation matrices)
  • Multiplicity calculation (|orbit| = group order / stabilizer order)
```

#### 1b: Kikuchi-Baker Coefficients

```
INPUTS:
  • Nij table: containment relationships (which clusters contain which)
  • Multiplicities from 1a

PROCESSING:
  NijTableCalculator.compute(disClusterList, disOrbitList)
    └─ For each pair (i,j): count how many orbits of type j
       are contained in clusters of type i

  KikuchiBakerCalculator.computeKikuchiBaker(mhdis, nijTable)
    └─ Solve recurrence:
       kb[t] = kb[t] - Σ_{s<t} (kb[s] · Nij[s][t])

OUTPUTS:
  • kb[tcdis]:  entropy coefficients
  • mhdis[]:    multiplicities
  • lc[]:       number of ordered clusters per HSP type
  • mh[][]      normalized multiplicities

ALGEBRA:
  • Matrix inversion (Nij inverse relationship)
  • Weighted sum computation
```

### STAGE 2: Correlation Function Identification

**File:** [CFIdentifier.java:92](app/src/main/java/org/ce/domain/cluster/CFIdentifier.java#L92)

#### 2a: HSP CF Enumeration (Full Decoration)

```
INPUTS:
  • disorderedClusters:        HSP geometry
  • disorderedSymOps:          HSP symmetry
  • basisN = ["s1", "s2"]      K-component basis (K-1 independent ops)

PROCESSING:
  ClusCoordListGenerator.generate(clusters, symOps, basisN)
    └─ Enumerate clusters with all K-component decorations
       (each site can have s1, s2, ... decorations)

OUTPUTS:
  • ClusCoordListResult (for CF decoration)
    ├─ tcfdis:   total CF types in HSP
    ├─ mcfdis[]: CF multiplicities

ALGEBRA:
  • Combinatorial enumeration of decorations
  • Symmetry-orbit classification
```

#### 2b: Ordered-Phase CF Grouping

```
INPUTS:
  • ClusterIdentificationResult (from Stage 1)
  • CF data for ordered phase (from 2a)
  • Transformation matrix (ordered → HSP frame)

PROCESSING:
  1. Generate ordered-phase CFs with full K-component basis
  2. Transform coordinates from ordered to HSP frame:
     OrderedToDisorderedTransformer.transform(coords, rotateMat, translateMat)
  3. Classify transformed CFs into HSP types:
     OrderedClusterClassifier.classify(transformed, hspClusters)
  4. Group CFs by cluster type and group:
     CFGroupGenerator.groupCFData(...)

OUTPUTS:
  • CFIdentificationResult
    ├─ tcf:      total CFs (including point CFs)
    ├─ ncf:      non-point CFs (optimization variables)
    ├─ nxcf:     point CFs (= K-1 for K components)
    ├─ lcf[][]   CFs per (type, group)
    └─ GroupedCFResult: CFs organized by cluster

ALGEBRA:
  • 3×3 matrix multiplication (rotation):
    coords_new = rotateMat · coords_old
  • Vector addition (translation):
    coords_new += translateMat
  • Classification via geometric distance/tolerance
```

### STAGE 3: C-Matrix Building

**File:** [CMatrixBuilder.java:18](app/src/main/java/org/ce\domain\cluster\CMatrixBuilder.java#L18)

#### Core Algorithm

```
FOR each ordered-phase cluster type t:
  FOR each group j in type t:
    cluster = ordClustersByType[t][j]

    1. Extract sites from cluster:
       siteIndices = flattenSiteIndices(cluster, siteList)
       (List<Integer> of site positions)

    2. Generate all possible occupancies (K choices per site):
       configs = generateConfigs(nsites, K)
       (e.g., for 3 sites & K=2: [A,A,A], [A,A,B], [A,B,A], ...)

    3. For each configuration:
       a) Compute CV via polynomial evaluation:
          row = computeCvRow(siteIndices, config, pRules, substituteRules, ...)

          HOW:
          • Use PRules to map config → polynomial
          • Use SubstituteRules to evaluate polynomial as CF linear combination
          • Result: row[0..tcf-1] = CF coefficients, row[tcf] = constant

       b) Create polynomial signature:
          key = PolynomialKey(row)  ◄─ canonical form of polynomial

       c) Count configurations mapping to same polynomial:
          countMap[key]++

       d) Store the polynomial:
          rowMap[key] = row

    4. Build C-matrix for this group:
       ncv = countMap.size()  ◄─ number of distinct CVs
       cmat[t][j] = double[ncv][tcf+1]  ◄─ each row is a polynomial
       wcv[t][j]  = int[ncv]             ◄─ weights (counts)

       FOR each unique polynomial in countMap:
         cmat[t][j][v] = polynomial row
         wcv[t][j][v]  = count of configs → this polynomial
```

#### Polynomial Computation Details

```
INPUTS (to computeCvRow):
  • siteIndices:     which sites are in this cluster
  • config:          occupation {A,B,A,...}
  • pRules:          p[i][s] = polynomial for site i with occupation s
  • substituteRules: how to multiply polynomials
  • cfColumn:        which CF index each polynomial maps to

ALGEBRA:
  1. Expand polynomial (p-rules):
     For each site i in cluster with occupation config[i]:
       get p[i][config[i]] as polynomial

     Multiply all these polynomials together:
       poly = p[site1] × p[site2] × ... × p[siteN]

  2. Simplify using substitution rules:
     Reduce cross-products using:
       s[a] · s[b] = substitution_rule[a][b]
       s[a]² = 1  (idempotent)

  3. Express as CF linear combination:
     poly = c0·CF0 + c1·CF1 + ... + ctcf·CFtcf

     Extract:
       row[k] = coefficient of CFk
       row[tcf] = constant term

OUTPUT:
  • row[0..tcf]: CF coefficients (the C-matrix row)
```

---

## Data Structures (Input & Output)

### ClusterIdentificationResult (Stage 1 Output)

```java
├─ tcdis (int)                      // # HSP cluster types (excl. point)
├─ tc (int)                         // # total clusters (incl. point)
├─ mhdis (List<Double>)             // HSP multiplicities [12, 24, 8, 6, 1, 0]
├─ kb (double[])                    // Kikuchi-Baker [1.0, 0.5, 0.25, ...]
├─ mh (double[][])                  // mh[t][j] - normalized mults
├─ lc (int[])                       // lc[t] - # ordered clusters per HSP type
├─ disClusterData (ClusterData)
│  ├─ clusList (List<Cluster>)
│  ├─ orbitList (List<List<Cluster>>)
│  └─ multiplicities (List<Double>)
└─ ordClusterData (ClassifiedClusterResult)
   ├─ clusters (List<List<Cluster>>)  // [t][j] organization
   └─ multiplicities (List<Double>)
```

### CFIdentificationResult (Stage 2 Output)

```java
├─ tcfdis (int)                     // # CF types in HSP
├─ tcf (int)                        // # total CFs (point + non-point)
├─ ncf (int)                        // # non-point CFs (optimization vars)
├─ nxcf (int)                       // # point CFs (= K-1)
├─ lcf (int[][])                    // lcf[t][j] - CFs per (type, group)
├─ disClusterData (ClusterData)
├─ disOrbitList (List<List<Cluster>>)
├─ ordClusterData (ClassifiedClusterResult)
└─ GroupedCFResult
   └─ groupedData[t][j]             // CFs for each (type, group)
```

### CMatrixResult (Stage 3 Output)

```java
├─ cmat (List<List<double[][]>>)    // cmat[t][j][v][k] - the C-matrix
│  where:
│    t = HSP cluster type (0..tcdis-1)
│    j = ordered cluster group (0..lc[t]-1)
│    v = cluster variable (0..lcv[t][j]-1)
│    k = correlation function (0..tcf)
│  Last column [k=tcf] is the constant term
│
├─ lcv (int[][])                    // lcv[t][j] - # CVs per group
├─ wcv (List<List<int[]>>)          // wcv[t][j][v] - CV weights
├─ cfBasisIndices (int[][])         // cfBasisIndices[col][rank] - decorations
└─ ... other topology data
```

---

## Key Algebra Classes

| Class | Purpose |
|-------|---------|
| **ClusCoordListGenerator** | Generate cluster orbits under symmetry |
| **NijTableCalculator** | Build containment table (clusters ⊂ clusters) |
| **KikuchiBakerCalculator** | Compute entropy coefficients from Nij |
| **CFGroupGenerator** | Organize CFs by cluster type & group |
| **SiteListBuilder** | Extract unique sites from cluster geometry |
| **PRulesBuilder** | Build polynomial rules for site occupancies |
| **CFSiteOpListBuilder** | Map CFs to site-operator sequences |
| **SubstituteRulesBuilder** | Build polynomial multiplication rules |
| **OrderedToDisorderedTransformer** | Apply transformation matrix to coordinates |
| **OrderedClusterClassifier** | Classify transformed clusters to HSP types |
| **PolynomialKey** | Canonical representation of polynomial for deduplication |

---

## Computational Complexity

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| **Cluster enumeration** | O(clusters × symmetry ops) | Manageable: ~10-100 clusters |
| **Nij table** | O(tcdis²) | Small: tcdis ~ 5-10 |
| **KB coefficients** | O(tcdis²) | Matrix inversion |
| **CF enumeration** | O(clusters × K^sites) | K^sites can be large! (2^4 = 16 for 4 sites) |
| **CV generation** | O(configs × polynomial_degree) | Configs = K^nsites per group |
| **Polynomial deduplication** | O(configs × log(configs)) | Hashing & grouping |
| **C-matrix final** | O(tcdis·lc·lcv·tcf) | Usually < 1 second |
| **Total Stage 1-3** | **O(K^max_sites)** | Dominated by exponential CV generation |

---

## Files to Study (Priority Order)

| Priority | File | Key Methods | Purpose |
|----------|------|-----------|---------|
| **1** | [Main.java](app/src/main/java/org/ce/ui/cli/Main.java) | `main():38` | CLI entry, TYPE-1a logic |
| **2** | [ClusterIdentificationWorkflow.java](app/src/main/java/org/ce/workflow/ClusterIdentificationWorkflow.java) | `identify():38` | Orchestrates all 3 stages |
| **3** | [ClusterIdentifier.java](app/src/main/java/org/ce/domain/cluster/ClusterIdentifier.java) | `identify():107` | Stage 1: cluster algebra |
| **4** | [CFIdentifier.java](app/src/main/java/org/ce/domain/cluster/CFIdentifier.java) | `identify():92` | Stage 2: CF algebra |
| **5** | [CMatrixBuilder.java](app/src/main/java/org/ce/domain/cluster/CMatrixBuilder.java) | `build():18` | Stage 3: polynomial algebra |
| **6** | [InputLoader.java](app/src/main/java/org/ce/storage/InputLoader.java) | `parseClusterFile()` | File I/O for clusters |
| **7** | [ClusCoordListGenerator.java](app/src/main/java/org/ce/domain/cluster/ClusCoordListGenerator.java) | `generate()` | Cluster orbit enumeration |
| **8** | [NijTableCalculator.java](app/src/main/java/org/ce/domain/cluster/NijTableCalculator.java) | `compute()` | Containment table |
| **9** | [KikuchiBakerCalculator.java](app/src/main/java/org/ce/domain/cluster/KikuchiBakerCalculator.java) | `computeKikuchiBaker()` | Entropy coefficients |
| **10** | [SubstituteRulesBuilder.java](app/src/main/java/org/ce/domain/cluster/SubstituteRulesBuilder.java) | `build()` | Polynomial multiplication rules |

---

## Execution Characteristics

### Typical Execution for BCC_A2/B2 (Binary)

```
Input files (2 files):
  ├─ clus/BCC_A2-T.txt       (20 KB, ~20 clusters)
  ├─ clus/BCC_B2-T.txt       (20 KB, ~20 clusters)
  ├─ symmetry/BCC_A2-SG      (parsed from configuration)
  └─ symmetry/BCC_B2-SG      (parsed from configuration)

Stage 1: ~100 ms
  • Cluster enumeration
  • Nij table (5×5)
  • KB coefficients

Stage 2: ~300 ms
  • CF enumeration (K^nsites exponential)
  • Classification & grouping

Stage 3: ~200 ms
  • Polynomial deduplication
  • C-matrix assembly

Total: ~600 ms, produces ~850 KB JSON file

Output (cluster_data.json):
  • tcdis = 5 (A2 clusters)
  • tcf = 16 (total CFs)
  • ncf = 15 (non-point CFs)
  • nxcf = 1 (point CFs, = K-1)
```

### Bottleneck: CF Enumeration

The exponential CF enumeration can be slow for:
- **Large clusters** (many sites per cluster)
- **Many components** (K=3, K=4 → K^nsites grows rapidly)

Example:
```
Ternary (K=3) with 4-site cluster:
  K^nsites = 3^4 = 81 configurations per group

Quaternary (K=4) with 4-site cluster:
  K^nsites = 4^4 = 256 configurations per group
```

Most time is spent in polynomial deduplication hashing.

---

## Summary: What TYPE-1a Produces

**For Use in TYPE-2 (Thermodynamic Calculation):**

The cluster_data.json file (AllClusterData) contains everything needed to:
1. ✅ Define the free-energy functional geometry
2. ✅ Compute Kikuchi-Baker entropy weights
3. ✅ Build the C-matrix for CV evaluation
4. ✅ Map correlation functions to chemical components
5. ✅ Specify cluster variables for the minimization loop

**All at zero thermodynamic cost** — this is pure geometry and combinatorics!
