# STAGE 4: Engine Setup Logging Detail

## STAGE 4 Components

Complete logging for CVMEngine.compute() with all substeps showing data extraction from each Stage.

---

## STAGE 4a: Validate Input

```
[INFO] === CVMEngine.compute() START ===
[INFO] Input parameters:
[INFO]   systemId: BCC_B2_T_bin
[INFO]   systemName: A-B_A2
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.500000, 0.500000]
[INFO]   numComponents: 2

[FINE] Validating C-Matrix...
[FINE] ✓ C-Matrix found

[FINE] Validating composition...
[FINE] ✓ Composition valid (sum=1.000000000)

[FINE] Validating temperature...
[FINE] ✓ Temperature valid
```

---

## STAGE 4b: Extract Stage 1-3 Data from AllClusterData

### Stage 1: Cluster Identification
```
[FINE] STAGE 4a: Extract Stage 1 (Cluster Identification)...
[FINE]   ✓ tcdis=5 (cluster types)
[FINE]   ✓ kb coefficients, mh multiplicities, lc counts loaded
```

**Data extracted:**
- `tcdis` = 5 (number of HSP cluster types)
- `kb[tcdis]` = Kikuchi-Baker entropy coefficients (5 values)
- `mhdis[tcdis]` = HSP cluster multiplicities (5 values)
- `mh[t][j]` = Normalized multiplicities per type and group
- `lc[t]` = Cluster count per HSP type (5 values)

### Stage 2: CF Identification
```
[FINE] STAGE 4b: Extract Stage 2 (CF Identification)...
[FINE]   ✓ tcf=6, ncf=4 (CF counts)
[FINE]   ✓ lcf array, CF basis indices loaded
```

**Data extracted:**
- `tcf` = 6 (total CFs: 4 non-point + 2 point)
- `ncf` = 4 (non-point CFs for optimization)
- `lcf[tcf]` = CF multiplicities per type
- `cfBasisIndices` = mapping for orthogonal basis (legacy, now null for CVCF)

### Stage 3: C-Matrix (CVCF Basis)
```
[FINE] STAGE 4c: Extract Stage 3 (C-Matrix - CVCF basis)...
[FINE]   ✓ cmat[t][j][v][k] transformation matrix loaded
[FINE]   ✓ lcv (CV counts), wcv (CV weights) loaded
[FINE]   ✓ cfBasisIndices: null (CVCF)
```

**Data extracted:**
- `cmat[t][j]` = C-matrix for each (cluster type, cluster group) pair
  - Dimensions: `[nv][tcf+1]` where last column is constant term
  - Shape: `cmat[5][variable][variable][7]` for BCC_B2
- `lcv[t][j]` = CV counts per (type, group)
- `wcv[t][j][v]` = CV weight per variable per (type, group)
- `cfBasisIndices` = **null** (indicates CVCF basis, not orthogonal)

---

## STAGE 4c: Create CVMInput Bundle

```
[FINE] STAGE 4d: Create CVMInput bundle...
[FINE]   ✓ CVMInput created with numComponents=2
```

**CVMInput wraps:**
- Stage 1 result (cluster topology)
- Stage 2 result (CF identification)
- Stage 3 result (C-matrix CVCF basis)
- System metadata (systemId, systemName)
- CVCF basis object (BccA2 binary basis)

---

## STAGE 4d: Evaluate ECI at Temperature

```
[FINE] Evaluating ECI at T=1000.0 K...
[FINE]   Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
[FINE]   CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
[FINE]   Processing 4 CEC terms...
[FINER]     eci[2] (v22AB) = -260.0 + 0.0*1000.0 = -2.6000e+02
[FINER]     eci[3] (v21AB) = -390.0 + 0.0*1000.0 = -3.9000e+02
[FINE]   ECI array: [[2]=-2.60e+02, [3]=-3.90e+02]
[FINE] ✓ ECI evaluated (4 non-point terms)
```

**ECI Evaluation:**
- For each CEC term with name like `e4AB`, `e3AB`, etc.
- Replace `e` with `v` to get CF name: `v4AB`, `v3AB`, etc.
- Lookup CF name in registry to get column index
- Evaluate: `eci[idx] = a + b * T`
  - `e22AB`: `a = -260.0`, `b = 0.0` → `eci[2] = -260.0 + 0.0*1000 = -260.0`
  - `e21AB`: `a = -390.0`, `b = 0.0` → `eci[3] = -390.0 + 0.0*1000 = -390.0`
- Missing terms default to 0 (direct inheritance property)

---

## STAGE 4e: Create CVMPhaseModel and Minimize

```
[INFO] Running CVM N-R minimization...
[INFO] CVMPhaseModel.ensureMinimized — STARTING minimization
[INFO]   systemId: BCC_B2_T_bin
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.5, 0.5]
[INFO]   ncf: 4, tcf: 6
[FINE]   tolerance: 1.0e-6
[FINE]   ECI: [[2]=-2.60e+02, [3]=-3.90e+02]

[INFO] NewtonRaphsonSolverSimple.solve — ENTER
[FINE]   ncf=4, tcf=6, tcdis=5, T=1000.0
[FINE]   moleFractions=[0.5, 0.5]
[FINE]   maxIter=400, tolerance=1.0e-6

[FINE] NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1000.0, xB=0.5, tolerance=1.0e-6
[FINE]   Initial CF (random state): [0.0625, 0.0, 0.25, 0.25]

[FINEST] NewtonRaphsonSolverSimple — iter   1: G=2.88145994e+03 H=0.00000000e+00 S=-2.88145994e+00 ||Gu||=1.87e-10

[FINE] NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=1, ||Gu||=1.87e-10

[INFO] CVMPhaseModel.minimize — SUCCESS: T=1000.0 K, x=[0.5, 0.5], iterations=1, ||dG||=1.87e-10, G=2881.4599 J/mol, elapsed=10 ms

[INFO] CVMPhaseModel.ensureMinimized — COMPLETE
```

---

## Complete STAGE 4 Output Example

```
[INFO] === CVMEngine.compute() START ===
[INFO] Input parameters:
[INFO]   systemId: BCC_B2_T_bin
[INFO]   systemName: A-B_A2
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.500000, 0.500000]
[INFO]   numComponents: 2

[FINE] Validating C-Matrix...
[FINE] ✓ C-Matrix found

[FINE] Validating composition...
[FINE] ✓ Composition valid (sum=1.000000000)

[FINE] Validating temperature...
[FINE] ✓ Temperature valid

[FINE] STAGE 4a: Extract Stage 1 (Cluster Identification)...
[FINE]   ✓ tcdis=5 (cluster types)
[FINE]   ✓ kb coefficients, mh multiplicities, lc counts loaded

[FINE] STAGE 4b: Extract Stage 2 (CF Identification)...
[FINE]   ✓ tcf=6, ncf=4 (CF counts)
[FINE]   ✓ lcf array, CF basis indices loaded

[FINE] STAGE 4c: Extract Stage 3 (C-Matrix - CVCF basis)...
[FINE]   ✓ cmat[t][j][v][k] transformation matrix loaded
[FINE]   ✓ lcv (CV counts), wcv (CV weights) loaded
[FINE]   ✓ cfBasisIndices: null (CVCF)

[FINE] STAGE 4d: Create CVMInput bundle...
[FINE]   ✓ CVMInput created with numComponents=2

[FINE] Evaluating ECI at T=1000.0 K...
[FINE]   Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
[FINE]   CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
[FINE]   Processing 4 CEC terms...
[FINER]     eci[2] (v22AB) = -260.0 + 0.0*1000.0 = -2.6000e+02
[FINER]     eci[3] (v21AB) = -390.0 + 0.0*1000.0 = -3.9000e+02
[FINE]   ECI array: [[2]=-2.60e+02, [3]=-3.90e+02]
[FINE] ✓ ECI evaluated (4 non-point terms)

[INFO] Running CVM N-R minimization...
[INFO] CVMPhaseModel.ensureMinimized — STARTING minimization
[INFO]   systemId: BCC_B2_T_bin
[INFO]   temperature: 1000.0 K
[INFO]   composition: [0.5, 0.5]
[INFO]   ncf: 4, tcf: 6
[FINE]   tolerance: 1.0e-6
[FINE]   ECI: [[2]=-2.60e+02, [3]=-3.90e+02]

[INFO] NewtonRaphsonSolverSimple.solve — ENTER
[FINE]   ncf=4, tcf=6, tcdis=5, T=1000.0
[FINE]   moleFractions=[0.5, 0.5]
[FINE]   maxIter=400, tolerance=1.0e-6

[FINE] NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1000.0, xB=0.5, tolerance=1.0e-6
[FINE]   Initial CF (random state): [0.0625, 0.0, 0.25, 0.25]

[FINEST] NewtonRaphsonSolverSimple — iter   1: G=2.88145994e+03 H=0.00000000e+00 S=-2.88145994e+00 ||Gu||=1.87e-10

[FINE] NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=1, ||Gu||=1.87e-10
[FINE]   Newton-Raphson solver returned
[FINE]   Cached 1 iteration snapshots
[FINE]   Solver converged!
[FINE]   Equilibrium CFs (non-point): [0.0625, 0.0, 0.25, 0.25]
[FINE]   Evaluating CVMFreeEnergy at equilibrium CFs...
[FINER] CVMFreeEnergy.evaluate — ENTER: T=1000.0, ncf=4, tcf=6
[FINER] CVMFreeEnergy.evaluate — EXIT: G=2.88145994e+03, H=0.00000000e+00, S=-2.88145994e+00
[FINE]   Thermodynamic values: G=2.88145994e+03, H=0.00000000e+00, S=-2.88145994e+00

[INFO] CVMPhaseModel.minimize — SUCCESS: T=1000.0 K, x=[0.5, 0.5], iterations=1, ||dG||=1.87e-10, G=2881.4599 J/mol, elapsed=10 ms

[INFO] CVMPhaseModel.ensureMinimized — COMPLETE

[INFO] ✓ CVM minimization converged in 1 iterations (||Gu||=1.87e-10)
[INFO]   G(eq) = 2.88145994e+03 J/mol
[INFO]   H(eq) = 0.00000000e+00 J/mol
[INFO]   S(eq) = -2.88145994e+00 J/(mol·K)

[INFO] === CVMEngine.compute() SUCCESS ===
```

---

## Key Logging Points in STAGE 4

| Substep | What's Logged | Level | Purpose |
|---------|---------------|-------|---------|
| **Validation** | C-matrix, composition, temperature checks | FINE | Verify input integrity |
| **Stage 1 Extract** | `tcdis`, KB coefficients, multiplicities | FINE | Show cluster topology data |
| **Stage 2 Extract** | `tcf`, `ncf`, CF basis info | FINE | Show CF identification data |
| **Stage 3 Extract** | C-matrix dimensions, CV counts, cfBasisIndices | FINE | Show CVCF transformation status |
| **CVMInput Create** | Bundle completion, numComponents | FINE | Verify input object created |
| **ECI Evaluation** | Basis info, CF names, per-term calculation | FINE/FINER | Trace ECI computation at temperature |
| **N-R Setup** | Solver parameters, initial state | FINE | Show minimization initialization |
| **N-R Iterations** | Per-iteration G, H, S, gradient norm | FINEST | Track convergence progress |
| **Final Results** | Equilibrium G, H, S | INFO | Report equilibrium thermodynamics |

---

## Summary

**STAGE 4 Logging Coverage: ✅ Complete**

- ✅ Input validation with checkmarks
- ✅ Stage 1 (Cluster) data extraction
- ✅ Stage 2 (CF) data extraction
- ✅ Stage 3 (C-matrix CVCF) data extraction
- ✅ CVMInput bundle creation
- ✅ ECI evaluation with per-term details
- ✅ N-R minimization setup and convergence
- ✅ Final thermodynamic results

All 8 substeps have detailed entry/exit logging with input and output parameters visible at FINE level.
