# Plan: Implement CVCF Basis for MCS ΔE Calculation

## Context

MCS currently computes ΔE using orthogonal site-operator products `u[t]` indexed by cluster type `t`, with `eci[t]` loaded positionally from the CEC file. The CEC file stores CVCF ECIs (e4AB, e3AB, e22AB, e21AB) by name, not by cluster-type position — so the current `extractEci()` mapping is **wrong**.

The correct approach: keep CVCF ECIs intact, and at ΔE computation time transform the orthogonal CF changes `Δu[t]` to CVCF CF changes `Δv[l]` using `BINARY_T_INV`, then dot with CVCF ECIs.

## Matrix Direction

- **T** (`BINARY_T`): rows=ortho `t`, cols=CVCF `j` → `u[t] = Σ_j T[t][j] · v[j]`  (CVCF → ortho)
- **Tinv** (`BINARY_T_INV`): rows=CVCF `l`, cols=ortho `t` → `v[l] = Σ_t Tinv[l][t] · u[t]`  (ortho → CVCF)

**Use Tinv** to go from orthogonal CFs to CVCF CFs.

## Formula

At each MC step when site `i` changes, the current code accumulates per-cluster-type contributions:
```
Δu[t] = Σ_{embeddings e of type t containing i}  (φ_new − φ_old) · Π_{k≠i} φ_k
```

Transform to CVCF CFs:
```
Δv[l] = Σ_t  Tinv[l][t] · Δu[t]      (l = 0..numNonPointCfs−1)
```
(`Tinv` = `BINARY_T_INV`, rows=CVCF, cols=ortho cluster type)

Then:
```
ΔE = Σ_l  eci_cvcf[l] · Δv[l]
```

This keeps the CVCF ECIs unchanged and makes CVCF CFs observable (can be logged/validated against CVM CFs).

## Implementation Steps

### Step 1: Verify `CvCfBasis` exposes `Tinv`

Read [CvCfBasis.java](app/src/main/java/org/ce/domain/cluster/cvcf/CvCfBasis.java).
If `Tinv` (`double[][]`) is not a public field, add it and populate it in `BccA2TModelCvCfTransformations.basisForNumComponents()` from `BINARY_T_INV` (already defined there).

### Step 2: Update `MCSEngine.extractEci()` → load CVCF ECIs by CF name

In [MCSEngine.java](app/src/main/java/org/ce/domain/engine/mcs/MCSEngine.java), change `extractEci()` to return an array indexed by CF index `l` (0..numNonPointCfs−1), not by cluster type:

```java
private static double[] extractEciCvcf(CECTerm[] terms, double temperature, CvCfBasis basis) {
    int ncf = basis.numNonPointCfs;
    double[] eci = new double[ncf];
    for (CECTerm term : terms) {
        if (term.name == null || !term.name.startsWith("e")) continue;
        String cfName = "v" + term.name.substring(1);   // e4AB → v4AB
        int idx = basis.cfNames.indexOf(cfName);
        if (idx >= 0 && idx < ncf) {
            eci[idx] = term.a + term.b * temperature;
        }
    }
    return eci;
}
```

Pass `eci_cvcf` and `basis.Tinv` (and `basis.numNonPointCfs`) down to `MCSRunner.Builder`.

### Step 3: Extend `MCSRunner.Builder` and constructor

Add fields:
```java
private double[]   eciCvcf;   // indexed by CF (l)
private double[][] tinv;      // [numCFs][tc]
private int        numNonPointCfs;
```
Add builder setters. In `MCSRunner.run()`, pass these to `MCEngine`/`ExchangeStep` in place of the old `eci[t]` array.

### Step 4: Modify `LocalEnergyCalc.deltaESingleSite()` (and `totalEnergy()`)

Accumulate `Δu[t]` per cluster type (exactly what it computes now), then:

```java
// transform ortho → CVCF
double[] deltaV = new double[numNonPointCfs];
for (int l = 0; l < numNonPointCfs; l++) {
    for (int t = 0; t < tc; t++) {
        deltaV[l] += tinv[l][t] * deltaU[t];
    }
}
// dot with CVCF ECIs
double dE = 0.0;
for (int l = 0; l < numNonPointCfs; l++) {
    dE += eciCvcf[l] * deltaV[l];
}
return dE;
```

Same pattern for `totalEnergy()`.

### Step 5: Verify cluster-type ordering matches BINARY_T rows

The transform assumes MCS cluster type `t=0` ↔ orthogonal CF row 0 (u[1]=tetra), `t=1` ↔ row 1 (u[2]=tri), etc.
Add a one-time `LOG.info` in `MCSEngine.compute()` that prints the orbit cluster sizes from `ClusCoordListResult.getOrbitList()` — confirm descending order matches: 4-body, 3-body, 2nd-pair, 1st-pair, point, point.

## Files to Modify

| File | Change |
|------|--------|
| [CvCfBasis.java](app/src/main/java/org/ce/domain/cluster/cvcf/CvCfBasis.java) | Add `Tinv` field if missing |
| [BccA2TModelCvCfTransformations.java](app/src/main/java/org/ce/domain/cluster/cvcf/BccA2TModelCvCfTransformations.java) | Pass `BINARY_T_INV` into `CvCfBasis` constructor |
| [MCSEngine.java](app/src/main/java/org/ce/domain/engine/mcs/MCSEngine.java) | Replace `extractEci()` with `extractEciCvcf()`; pass tinv + numNonPointCfs to builder |
| [MCSRunner.java](app/src/main/java/org/ce/domain/engine/mcs/MCSRunner.java) | Add `eciCvcf`, `tinv`, `numNonPointCfs` fields; propagate to MCEngine/ExchangeStep |
| [LocalEnergyCalc.java](app/src/main/java/org/ce/domain/engine/mcs/LocalEnergyCalc.java) | Accumulate `Δu[t]`, transform → `Δv[l]`, dot with `eciCvcf[l]` |
| [ExchangeStep.java](app/src/main/java/org/ce/domain/engine/mcs/ExchangeStep.java) | Accept updated LocalEnergyCalc signature |

## Verification

1. **Ordering log**: Confirm cluster-type ordering in log output matches BINARY_T row order.
2. **CVM cross-check**: Same composition/temperature (e.g. Nb-Ti, x=0.5, T=1000 K) — MCS `⟨H⟩/site` should match CVM `H/site` within statistical error.
3. **Zero-ECI**: With all ECIs = 0, both CVM and MCS return H = 0.
4. **CF output**: Log `Δv[l]` averages during MCS run; compare against CVM equilibrium `v[l]` values.
