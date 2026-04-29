# AlloyMC Engine Implementation Plan

Implementation of a high-performance Metropolis Monte Carlo engine for Cluster Expansion models in CEWorkbench.

## Current Status: BLOCKED (Ternary Mapping)

> [!CAUTION]
> **Ternary Verification Discrepancy**
> While the engine is 100% accurate for Binary systems (Nb-Ti), the Ternary system (Nb-Ti-V) exhibits a mapping error at the pure Nb endpoint.
> - **Expected**: All CVCF non-point terms = 0.0 for pure Nb.
> - **Actual**: `v3AB`, `v22AB`, and `v21AB` result in `-2.0`.
> - **Impact**: Resulting mixing energy is physically incorrect (-2.4 MJ/mol).

## Root Cause Analysis
The discrepancy indicates a mismatch between the **Orthogonal CF Vector (u)** measured from the lattice and the **Transformation Matrix (T^-1)** columns.
- Possible Column Shift: The point functions or constant term indices in the transformation matrix may not align with the `[NonPoint | Point | Constant]` layout in `AlloyMC`.
- Basis Inconsistency: `ClusterMath.buildBasis(3)` values `[1, 0, -1]` must be confirmed against the `CvCfBasis` derivation assumptions.

## Proposed Resolution Steps

### 1. Diagnostic Mapping Trace
Implement a detailed trace in `AlloyMC` to print the individual contributions of each orthogonal term to the problematic CVCF functions (`v3AB`, `v21AB`).
- Identify which orthogonal column index is introducing the `-2.0` weight.

### 2. Standardize Basis Layout
Ensure the `correlationFunctions` vector in `AlloyMC` strictly follows the `totalCfs` ordering defined by the `PipelineResult`:
- Verify if `ncf` includes or excludes point orbits in the `Tinv` matrix derivation.

### 3. Implement Robust Delta-E (After Fix)
Once mapping is verified:
- Implement `calculateDeltaE(int i, int j)` using `siteToCfIndex` CSR structure for $O(1)$ updates.
- Implement `runSweep()` with local energy tracking.

## Verification Plan
### Automated Tests
- `AlloyMCTest`: Must produce 0.0 J for pure Nb in the ternary Nb-Ti-V system.
- `AlloyMCTest`: Verify randomized ternary energy is plausible.

### Manual Verification
- Trace the $T^{-1}$ matrix construction for a ternary tetrahedron to ensure linear independence and correct endpoint centering.
