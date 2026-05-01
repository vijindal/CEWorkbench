# AlloyMC Engine Implementation Plan

Implementation of a high-performance Metropolis Monte Carlo engine for Cluster Expansion models in CEWorkbench.

## Current Status: ACTIVE

> [!NOTE]
> **Ternary Mapping Resolved**
> The previously identified -2.0 energy discrepancy in the Nb-Ti-V system has been resolved. The root cause was a basis function indexing mismatch between the engine's state vector population and the CVCF transformation matrix.

## Key Resolutions
- **Dynamic Basis Mapping**: Refactored `AlloyMC` to use `MCSGeometry.getCfBasisIndices()` metadata for point CF column assignments.
- **Full Vector Transformation**: Updated `MCSGeometry` to transform the entire measured orthogonal vector ($T^{-1}u$) rather than hard-appending compositions.
- **Symbolic Verification**: Integrated symbolic definitions ($p[A]*p[B]$) directly into the `CvCfBasis` and `AlloyMCTest` diagnostic suite.

## Roadmap
### 1. Finalize Hamiltonian Integration
- Integrate the verified CVCF state vector into the `AlloyMC` energy calculator.
- Verify mixing energy calculation for the Nb-Ti-V equiatomic disordered state.

### 2. Metropolis Sampling Loop
- Implement the standard Metropolis-Hastings acceptance criterion.
- Optimize site-swap logic for local energy updates.

### 3. Verification & Benchmarking
- Validate against reference CVM results for simple BCC alloys.
- Benchmark site-swap performance for large supercells ($L > 20$).

## Verification Plan
### Automated Tests
- `AlloyMCTest`: Must produce 0.0 J for pure Nb in the ternary Nb-Ti-V system.
- `AlloyMCTest`: Verify randomized ternary energy is plausible.

### Manual Verification
- Trace the $T^{-1}$ matrix construction for a ternary tetrahedron to ensure linear independence and correct endpoint centering.
