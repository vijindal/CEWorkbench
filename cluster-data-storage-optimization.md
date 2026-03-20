# Plan: Cluster Data Storage Optimization

## Baseline (BCC_B2_T_bin/cluster_data.json)

| Section | Bytes | % |
|---|---|---|
| `disorderedClusterResult` | 73,383 | 21.2% |
| `orderedClusterResult` | 73,383 | 21.2% |
| `disorderedCFResult` | 98,102 | 28.4% |
| `orderedCFResult` | 98,102 | 28.4% |
| `cmatrixResult` | 2,799 | 0.8% |
| **TOTAL** | **345,888** | |

Python confirms: `orderedClusterResult == disorderedClusterResult` and `orderedCFResult == disorderedCFResult` — exact duplicate objects stored twice.

---

## What each engine actually reads at runtime

### CVM (`CVMPhaseModel`)
```
disorderedClusterResult.disClusterData.multiplicities   ← 31 bytes
disorderedClusterResult.kbCoefficients                  ← 27 bytes
disorderedClusterResult.nijTable                        ← 85 bytes  [kept for KB recomputation]
disorderedClusterResult.mh / lc / tcdis / nxcdis / tc  ← ~69 bytes
disorderedCFResult.lcf / tcf / nxcf / ncf / tcfdis     ← ~28 bytes
cmatrixResult  (all 4 fields)                           ← 2,799 bytes
```

### MCS (`EmbeddingGenerator`)
```
disorderedClusterResult.disClusterData.orbitList        ← 21,234 bytes (raw)
disorderedClusterResult.disClusterData.tc               ← 1 byte
  — from each Cluster: sublattices → sites → symbol + position
```

### Neither engine reads
- `orderedClusterResult` — entire section
- `orderedCFResult` — entire section
- `disorderedClusterResult.phaseClusterData` — raw ordered clusters before CF identification
- `disorderedClusterResult.ordClusterData` — classified ordered clusters (CVM only uses derived `mh`/`lc`)
- `disorderedCFResult.disCFData` — CF cluster geometry (CVM only uses scalar summaries)
- `disorderedCFResult.phaseCFData`, `ordCFData`, `groupedCFData` — CMatrix build-time intermediates
- `disorderedClusterResult.disClusterData.clusCoordList` — canonical cluster reps (neither engine)
- `disorderedClusterResult.disClusterData.rcList` — sublattice site counts (neither engine)
- `disorderedClusterResult.disClusterData.numPointSubClusFound` — CF normalization count

---

## Bug found: `Cluster.getAllSites()` serialized by accident

`getAllSites()` is a computed method — it builds a flat `List<Site>` by iterating `sublattices`:

```java
public List<Site> getAllSites() {
    List<Site> all = new ArrayList<>();
    for (Sublattice sub : sublattices) all.addAll(sub.getSites());
    return all;
}
```

Because it follows the JavaBean `getXxx()` naming convention, Jackson serializes it as `"allSites"`.
It appears **844 times** in the file, duplicating all site data already present in `"sublattices"`.

`allSites` in kept sections (disClusterData.orbitList):
- Raw `orbitList`: 21,234 bytes
- `allSites` fraction ≈ 50% → ~10,617 bytes of pure waste
- After removing: `orbitList` shrinks to ~10,617 bytes

Fix: add `@JsonIgnore` to `getAllSites()`. The method continues to work at runtime; it just isn't persisted.
Also add `@JsonIgnoreProperties(ignoreUnknown = true)` to `Cluster` so old JSON files with `"allSites"` keys load without error.

---

## Complete changes

### A — Stop serializing duplicate top-level objects (add `@JsonIgnore` to getter)

| File | Getter | Why |
|---|---|---|
| `AllClusterData.java` | `getOrderedClusterResult()` | Always same instance as disordered; never accessed at runtime |
| `AllClusterData.java` | `getOrderedCFResult()` | Same |

### B — Stop serializing workflow-only nested fields (add `@JsonIgnore` to getter)

| File | Getter | Why |
|---|---|---|
| `ClusterIdentificationResult.java` | `getPhaseClusterData()` | Raw ordered clusters before CF identification — workflow only |
| `ClusterIdentificationResult.java` | `getOrdClusterData()` | Classified ordered clusters — CMatrix build only; CVM uses `mh`/`lc` |
| `CFIdentificationResult.java` | `getDisCFData()` | CF cluster geometry — CVM needs only scalars (tcf/ncf/lcf) |
| `CFIdentificationResult.java` | `getPhaseCFData()` | Workflow only |
| `CFIdentificationResult.java` | `getOrdCFData()` | Workflow only |
| `CFIdentificationResult.java` | `getGroupedCFData()` | CMatrix build only |
| `ClusCoordListResult.java` | `getClusCoordList()` | Neither CVM nor MCS uses canonical reps at runtime |
| `ClusCoordListResult.java` | `getRcList()` | Neither engine uses sublattice site counts |
| `ClusCoordListResult.java` | `getNumPointSubClusFound()` | CF normalization count — workflow only |

### C — Fix accidental getter serialization (add `@JsonIgnore` to getter)

| File | Getter | Why |
|---|---|---|
| `Cluster.java` | `getAllSites()` | Computed method serialized by accident; 844 occurrences in file |

### D — Add `@JsonIgnoreProperties(ignoreUnknown = true)` for backward compatibility

All classes that will have `@JsonIgnore` getters (so old JSON files with those fields still load):

`Cluster`, `ClusterIdentificationResult`, `CFIdentificationResult`, `ClusCoordListResult`, `ClassifiedClusterResult`, `GroupedCFResult`, `CMatrixResult`

### E — `@JsonProperty` constructor params remain unchanged (backward compat read path)

For every dropped field, **keep** its `@JsonProperty` annotation on the `@JsonCreator` constructor parameter. Jackson uses constructors for deserialization — old JSON files with those fields will continue to deserialize normally. New JSON files simply won't have those fields, so the constructor receives `null` for those params (which is fine since they're never accessed at runtime).

---

## Projected file sizes

| Section | Current | After |
|---|---|---|
| `disorderedClusterResult` | 73,383 | ~10,830 |
| `orderedClusterResult` | 73,383 | **0** (dropped) |
| `disorderedCFResult` | 98,102 | ~28 |
| `orderedCFResult` | 98,102 | **0** (dropped) |
| `cmatrixResult` | 2,799 | 2,799 |
| **TOTAL** | **345,888** | **~13,657** |

**Reduction: ~96% (25× smaller)**

---

## Notes

- **KB coefficients**: `nijTable` is kept. CVM reads pre-computed `kbCoefficients` directly (`CVMPhaseModel:248`); it never re-derives them from `nijTable`. `nijTable` is retained as the source material needed to recompute KB coefficients without a full type-1a re-run.
- **C-Matrix**: `cmatrixResult` is exponentially expensive to recompute — kept in full.
- **No logic changes**: All changes are annotations only. No computation paths or data structures change.

---

## Verification

1. `./gradlew compileJava` — `BUILD SUCCESSFUL`
2. Run type-1a (cluster identification) → inspect new `cluster_data.json` size
3. CVM single-point calculation end-to-end → non-zero Gibbs energy
4. MCS single-point calculation → sweep progress visible
5. Load an old `cluster_data.json` (with the extra fields) → deserializes without error
