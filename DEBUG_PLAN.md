# Debug Plan: RandomStateEvaluator — Wrong Results vs TernaryCMatrixRunner

## Context

`TernaryCMatrixRunner` and `RandomStateEvaluator` both compute the disordered (random)
state CFs for Nb-Ti-V / BCC_A2 / T. `TernaryCMatrixRunner` passes (all CVs ∈ (0,1),
S ≈ ideal). `RandomStateEvaluator` produces wrong results: large negative CVCF CFs
→ negative cluster variables → entropy blow-up (S = −1,304,093 J/(mol·K)).

The two runners take different code paths through the same underlying pipeline. This
plan traces those paths step-by-step, adds printouts to expose where values diverge,
and identifies the root cause.

---

## Root Cause Hypothesis

The two runners differ in **two key ways**:

### Difference 1 — C-matrix used for CV evaluation

| Runner | C-matrix source | What it is |
|--------|----------------|------------|
| `TernaryCMatrixRunner` | `cMatData` from `CMatrixPipeline.run()` → passed directly to `verifyRandomCVs()` | **Orthogonal** basis C-matrix |
| `RandomStateEvaluator` | `this.cmat = mdCvcf.getCmat()` inside `CVMGibbsModel.initialize()` | **CVCF-transformed** C-matrix (`C_orth · T`) |

`verifyRandomCVs()` (`CMatrixPipeline.java:818`) calls `evaluateCVs(uFull, cmatData.cmat, ...)`
with the **orthogonal** C-matrix.
`evaluateClusterVariables()` (`CVMGibbsModel.java:454`) calls `evaluateCVs(uFull, this.cmat, ...)`
with the **CVCF** C-matrix.

### Difference 2 — CF vector passed to the C-matrix

| Runner | CF vector builder | What it produces |
|--------|------------------|------------------|
| `TernaryCMatrixRunner` | `pipelineResult.computeRandomCFs()` | **Orthogonal** CFs: products of point φₖ values |
| `RandomStateEvaluator` | `basis.computeRandomStateVectors()` → `rv.uCvcfNonPoint` → `buildFullCVCFVector()` | **CVCF** CFs: T⁻¹·uOrthFull, then mole fractions appended |

**Expected pairing:** orthogonal C-matrix × orthogonal CF vector = CVs.
Also CVCF C-matrix × CVCF CF vector = CVs (equivalent via T transformation).
`TernaryCMatrixRunner` uses the orthogonal pair. `RandomStateEvaluator` is supposed
to use the CVCF pair — but if T or T⁻¹ is wrong, the CVCF CFs will be wrong and
produce garbage CVs.

### Difference 3 — Point vars in the full vector

`pipelineResult.computeRandomCFs()` appends **orthogonal point CFs** (φₖ values) at
positions `[ncf..ncf+K-2]` and `1.0` at `[ncf+K-1]` — total length `ncf + K`.

`buildFullCVCFVector()` appends **mole fractions** (x₀, x₁, x₂) at positions
`[ncf..ncf+K-1]` — total length `ncf + K`.

This is by design if the C-matrix columns are aligned accordingly. **If the column
count of `this.cmat` doesn't match the vector length, the matrix multiply silently
produces garbage.**

---

## Debugging Steps

### Step 1 — Add cross-check block in `RandomStateEvaluator`

**File:** `src/main/java/org/ce/debug/RandomStateEvaluator.java`

After the existing STEP 6 CV block, add a **CROSS-CHECK** section that recomputes CVs
via the orthogonal path (as `TernaryCMatrixRunner` does) and prints both side by side:

- Run `CMatrixPipeline.run()` independently to get `mdOrth`
- Call `pipelineResult.computeRandomCFs(x)` to get the orthogonal full CF vector
- Call `CMatrixPipeline.evaluateCVs(uOrthFull, mdOrth.getCmat(), ...)` to get CVs
- Print both sets of CVs (CVCF path vs orthogonal path) side by side

### Step 2 — Print C-matrix and vector dimensions

In the cross-check block, print:
- Column count of `this.cmat` (CVCF path)
- Column count of `mdOrth.getCmat()` (orthogonal path)
- Length of `uFull` vector in each path

If column counts differ from vector lengths → matrix/vector alignment bug.

### Step 3 — Compare CF vectors

Print `rv.uCvcfNonPoint` (CVCF path) alongside the non-point slice of
`pipelineResult.computeRandomCFs(x)` (orthogonal path).

If T is correct, applying T⁻¹ to the orthogonal CFs should give back the CVCF CFs,
and both paths should produce identical CVs. If they differ, T is wrong.

---

## Files to Modify

| File | Change |
|------|--------|
| `src/main/java/org/ce/debug/RandomStateEvaluator.java` | Add CROSS-CHECK block: orthogonal-path CV recomputation, dimension printouts, side-by-side CF comparison |
| `src/main/java/org/ce/model/cvm/CVMGibbsModel.java` | Expose `getPipelineResult()` accessor (store `PipelineResult pr` as a field during `initialize()`) |

---

## Verification

```bash
# Run the debug evaluator — look for CROSS-CHECK section output
./gradlew randomEval

# If CVs match in cross-check but not in STEP 6 → bug is in CVCF CF vector (T⁻¹ wrong)
# If CVs also wrong in cross-check → bug is upstream (pipeline, cfBasisIndices, composition)
# If dimensions differ → matrix/vector alignment bug

# TernaryCMatrixRunner for reference (should produce all CVs ∈ (0,1))
./gradlew runTernaryDebug
```
