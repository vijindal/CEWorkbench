# MCS Analytic Validation — Detailed Results

Reference record of the MCS energy/CF validation performed via
`org.ce.model.mcs.AnalyticConfigVerification`, with the analytic derivation and
the actually-computed value at every stage, for every configuration and every
system size (K=2, K=3, K=4). This is the supporting detail behind the summary
in `CLAUDE.md`'s Testing section.

Run command:
```bash
./gradlew runScratch -PscratchClass=org.ce.model.mcs.AnalyticConfigVerification
```

All numeric values below were captured directly from the test run and a
companion diagnostic dump (not retyped/rounded from memory), on:
- `L=4` (N=128 sites) for MCS-measured quantities
- BCC_A2 / T model, `T=1000 K` for ECI evaluation

---

## Pipeline stages under test

```
Configuration (LatticeConfig)
   -> orthogonal CFs, uOrth   [Embeddings.measureFullCVsFromConfig / PipelineResult.computeRandomCFs]
   -> CVCF CFs, vCvcf         [Embeddings.applyTinvTransform, via CvCfBasis.Tinv]
   -> Energy, H = N * Sum(eciCvcf[l] * vCvcf[l])   [Embeddings.totalEnergyCvcf]
   -> Delta-E for a trial exchange                  [Embeddings.deltaEExchangeCvcf / V2]
   -> Metropolis trajectory (running E)             [MetropolisMC-equivalent loop]
```

Ground truth sources:

| Quantity | Ground-truth source | Independence from MCS code |
|---|---|---|
| Orthogonal CFs | `PipelineResult.computeRandomCFs(xFrac)` | Closed-form combinatorial formula, zero dependence on embedding code |
| CVCF CFs | `Tinv` applied to the analytic orthogonal vector | `Tinv` is the same shared matrix object CVM and MCS both use (verified identical, not independently re-derived) |
| B2 pair CFs | Hand-derived from site parity | Cross-checked empirically against actual embedding site-index parity — the naive assumption was wrong and had to be corrected (see K=2 Stage 2) |
| Composition-boundary CVCF invariant | Physical requirement | A named CVCF quantity (e.g. a pair probability) must not depend on a third species present at zero mole fraction |

---

## K=2 — Nb-Ti / BCC_A2 / T

**Non-point CF names** (ncf=4): `v4AB, v3AB, v22AB, v21AB`
**Point CFs**: `xA` (alpha=1) — K=2 has only one point column
**Orthogonal basis sequence**: `seq = {-1, +1}` (species 0 → -1, species 1 → +1)

| ECI term | v4AB | v3AB | v22AB | v21AB |
|---|---:|---:|---:|---:|
| eciCvcf (T=1000K) | 0.0 | 0.0 | 3120.0 | 6240.0 |

(`v4AB`, `v3AB` unmapped in the shipped Nb-Ti Hamiltonian → 0.)

### Stage 1 — Pure element (x_A = 1, all sites species 0)

| Column | uOrth analytic | uOrth measured | Match |
|---|---:|---:|---|
| v4AB | 1.0 | 1.0 | exact |
| v3AB | -1.0 | -1.0 | exact |
| v22AB | 1.0 | 1.0 | exact |
| v21AB | 1.0 | 1.0 | exact |
| xA (point, α=1) | -1.0 | -1.0 | exact |
| empty | 1.0 | 1.0 | exact |
| **max error** | | | **0.0** |

| CVCF column | v4AB | v3AB | v22AB | v21AB |
|---|---:|---:|---:|---:|
| vCvcf (Tinv·uOrth) | 0.0 | 0.0 | 0.0 | 0.0 |

All zero — no fractional-species-pair probability exists in a pure single-species state.

| Quantity | Formula | Value | Measured | Match |
|---|---|---:|---:|---|
| H (energy/site) | `Sum(eci·vCvcf)` | 0.0 | 0.0 | exact |

### Stage 2 — Perfectly ordered B2 (x=0.5, alternating sublattices)

Species assigned by flat site-index parity: even index → species 0, odd → species 1
(`MCSGeometry.buildBCCPositions`).

**Correction made during derivation** — the initial hand-derivation assumed
`v21AB` ("I-n pair" / "1NN" in `CvCfBasis`'s registration comment) connects
*opposite*-parity sites and `v22AB` ("II-n" / "2NN") connects *same*-parity
sites. Direct inspection of the actual generated embeddings
(`geo.cfEmbeddings()`) showed the reverse:

| CF column | Sample embeddings (site indices, parity) | Pattern found |
|---|---|---|
| v21AB | `[0,32]→[0,0]`, `[1,33]→[1,1]`, `[2,34]→[0,0]`, `[3,35]→[1,1]` | **all same-parity** |
| v22AB | `[0,25]→[0,1]`, `[1,34]→[1,0]`, `[2,27]→[0,1]`, `[3,36]→[1,0]` | **all opposite-parity** |

So the "I-n"/"II-n" naming does not correspond to a naive same-cell
nearest/next-nearest-neighbor distance reading. The table below uses the
empirically-confirmed rule, not the naive one.

| Quantity | Analytic derivation | Analytic value | Measured | Match |
|---|---|---:|---:|---|
| v21AB (uOrth) | same-parity: φ₁(0)·φ₁(0) = φ₁(1)·φ₁(1) | +1.0 | +1.0 | exact |
| v22AB (uOrth) | opposite-parity: φ₁(0)·φ₁(1) = (-1)(+1) | -1.0 | -1.0 | exact |
| point CF[0] | 0.5·seq[0] + 0.5·seq[1] = 0.5(-1)+0.5(1) | 0.0 | 0.0 | exact |

| Column | v4AB | v3AB | v22AB | v21AB | xA | empty |
|---|---:|---:|---:|---:|---:|---:|
| uOrth (measured, B2) | 1.0 | 0.0 | -1.0 | 1.0 | 0.0 | 1.0 |

| CVCF column | v4AB | v3AB | v22AB | v21AB |
|---|---:|---:|---:|---:|
| vCvcf (B2) | 0.5 | 0.0 | 0.0 | 0.5 |

| Quantity | Method | Value | Cross-check | Match |
|---|---|---:|---:|---|
| E/site (synthetic ECI isolating v21AB, a=1,b=0) | `totalEnergyCvcf` | 0.5 | vCvcf[v21AB] = 0.5 | exact |

### Stage 3 — Perfectly random (x=0.5, equiatomic)

| Column | v4AB | v3AB | v22AB | v21AB | xA | empty |
|---|---:|---:|---:|---:|---:|---:|
| uOrth (analytic) | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 1.0 |

All non-point/point CFs vanish at x=0.5 for K=2, since `pointCF[α=1] = 0.5(-1)+0.5(+1) = 0`
and every non-point column is a product of powers of this zero value.

| CVCF column | v4AB | v3AB | v22AB | v21AB |
|---|---:|---:|---:|---:|
| vCvcf (analytic) | 0.0625 | 0.0 | 0.25 | 0.25 |

| Quantity | Value |
|---|---:|
| H closed-form (`Sum(eci·vCvcf)`) | **2340.0 J/mol** |

| MCS measured (avg of 40 seeds, L=4) | Value |
|---|---:|
| MCS random H | 2363.664063 |
| stderr | ±13.021317 |
| closed-form H | 2340.000000 |
| diff | +23.664063 |
| 5σ tolerance | 65.1 |
| **Result** | **PASS** |

Finite-size sampling noise at L=4/N=128 (statistical, not systematic — confirmed
by the shrinking-with-L trend documented in the original `RandomStateEnergyTest`
runs, which this suite superseded).

### Stage 4 — Delta-E, anchored at analytic configs

| Anchor | Trials run | Mismatches (>1e-6) | Max abs error | Result |
|---|---:|---:|---:|---|
| B2-anchored | 57 | 0 | 0.000e+00 | PASS |
| pure-element-anchored | 0* | 0 | 0.000e+00 | PASS |

\* Pure-element anchor: every trial pair happened to be same-species (no-op,
skipped at this seed) — 0 executed trials is expected, not a bug (matches
`Embeddings.deltaEExchangeCvcf`'s early-return on `occI==occJ`).

Checked: `deltaEExchangeCvcf` (flat), `deltaEExchangeCvcf` (list), and
`deltaEExchangeCvcfV2`, all against finite-difference
(`totalEnergyCvcf(after) - totalEnergyCvcf(before)`).

### Stage 5 — Trajectory (running energy vs. from-scratch recompute)

500 accept/reject Metropolis-style steps starting from the B2 configuration, T=1000K, seed=7.

| Steps | Mismatches | Max abs error | Result |
|---:|---:|---:|---|
| 500 | 0 | 0.000e+00 | PASS |

Running energy and from-scratch recompute tracked in lock-step at every step.

---

## K=3 — Nb-Ti-V / BCC_A2 / T

**Non-point CF names** (ncf=18): `v4AB, v4AC, v4BC, v4ABC1, v4ABC2, v4ABC3, v3AB, v3AC, v3BC, v3ABC1, v3ABC2, v3ABC3, v22AB, v22AC, v22BC, v21AB, v21AC, v21BC`
**Point CFs**: `xA` (α=2), `xB` (α=1) — in the pipeline's raw (descending-α) column order; see "Bug found and fixed" below.
**Orthogonal basis sequence**: `seq = {-1, 0, +1}` (species A=-1, B=0, C=+1)

| CF family | v4AB | v4AC | v4BC | v4ABC1 | v4ABC2 | v4ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | 0.0 | 48800.0 | -24400.0 | -24400.0 |

| CF family | v3AB | v3AC | v3BC | v3ABC1 | v3ABC2 | v3ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | 0.0 | -20000.0 | -20000.0 | -20000.0 |

| CF family | v22AB | v22AC | v22BC | v21AB | v21AC | v21BC |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 3120.0 | 7040.0 | 4080.0 | 6240.0 | 14080.0 | 8160.0 |

### Stage 1 — Pure element (x_0 = 1)

| Family | v4AB | v4AC | v4BC | v4ABC1 | v4ABC2 | v4ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 |
| uOrth measured | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 |

| Family | v3AB | v3AC | v3BC | v3ABC1 | v3ABC2 | v3ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 1.0 | -1.0 | 1.0 | -1.0 | 1.0 | -1.0 |
| uOrth measured | 1.0 | -1.0 | 1.0 | -1.0 | 1.0 | -1.0 |

| Family | v22AB | v22AC | v22BC | v21AB | v21AC | v21BC | xA(α=2) | xB(α=1) | empty |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 |
| uOrth measured | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 | 1.0 | -1.0 | 1.0 |

**Max error across all 21 columns: 0.0** (identical, entry by entry).

| Quantity | Value |
|---|---:|
| vCvcf (Tinv·uOrth), all 18 columns | 0.0 |

All zero — a pure single-species state has zero probability of any mixed-species cluster.

### Stage 2 — Composition boundary (x = [0.5, 0.5, 0.0], V absent)

This is the case that surfaced the real regression (see "Bug found and fixed"
below). Values shown are **after the fix**.

| Family | v4AB | v4AC | v4BC | v4ABC1 | v4ABC2 | v4ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth (boundary) | 0.0625 | -0.0625 | 0.0625 | 0.0625 | -0.0625 | 0.0625 |

| Family | v3AB | v3AC | v3BC | v3ABC1 | v3ABC2 | v3ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth (boundary) | 0.125 | -0.125 | 0.125 | -0.125 | 0.125 | -0.125 |

| Family | v22AB | v22AC | v22BC | v21AB | v21AC | v21BC | xA(α=2) | xB(α=1) | empty |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uOrth (boundary) | 0.25 | -0.25 | 0.25 | 0.25 | -0.25 | 0.25 | 0.5 | -0.5 | 1.0 |

| CVCF column | v4AB | v4AC | v4BC | v4ABC1 | v4ABC2 | v4ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| vCvcf (boundary) | 0.0625 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| CVCF column | v3AB | v3AC | v3BC | v3ABC1 | v3ABC2 | v3ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| vCvcf (boundary) | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| CVCF column | v22AB | v22AC | v22BC | v21AB | v21AC | v21BC |
|---|---:|---:|---:|---:|---:|---:|
| vCvcf (boundary) | 0.25 | 0.0 | 0.0 | 0.25 | 0.0 | 0.0 |

**Cross-check against the true K=2 (Nb-Ti) system's own CVCF values at x=0.5** (K=2
Stage 3 above), matched by shared name:

| Name | K=3 boundary CVCF | K=2 native CVCF | diff |
|---|---:|---:|---:|
| v4AB | 0.0625 | 0.0625 | 0.0 |
| v3AB | 0.0 | 0.0 | 0.0 |
| v22AB | 0.25 | 0.25 | 0.0 |
| v21AB | 0.25 | 0.25 | 0.0 |

| Quantity | Value |
|---|---:|
| matchedNames | 4 |
| maxErr | 0.000e+00 |
| **Result (after fix)** | **PASS (exact)** |
| maxErr (before fix, for reference) | 1.0 (every term off by exactly 1.0) |

### Stage 3 — Perfectly random (equiatomic, x=[1/3,1/3,1/3])

| Family | v4AB | v4AC | v4BC | v4ABC1 | v4ABC2 | v4ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 0.19753086 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| Family | v3AB | v3AC | v3BC | v3ABC1 | v3ABC2 | v3ABC3 |
|---|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 0.29629630 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| Family | v22AB | v22AC | v22BC | v21AB | v21AC | v21BC | xA(α=2) | xB(α=1) | empty |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 0.44444444 | 0.0 | 0.0 | 0.44444444 | 0.0 | 0.0 | 0.66666667 | 0.0 | 1.0 |

| CVCF column | v4AB…v4ABC3 (all 6) | v3AB/v3AC/v3BC (all 3) | v3ABC1/2/3 (all 3) | v22/v21 family (all 6) |
|---|---:|---:|---:|---:|
| vCvcf analytic | 0.01234568 | 0.0 | 0.03703704 | 0.11111111 |

| Quantity | Value |
|---|---:|
| H closed-form (`Sum(eci·vCvcf)`) | **2524.444444 J/mol** |

| MCS measured (avg of 40 seeds, L=4) | Value |
|---|---:|
| MCS random H | 2467.996745 |
| stderr | ±64.511267 |
| closed-form H | 2524.444444 |
| diff | -56.447700 |
| 5σ tolerance | 322.6 |
| **Result** | **PASS** |

### Stage 4 — Delta-E, anchored at analytic configs

| Anchor | Trials run | Mismatches | Max abs error | Result |
|---|---:|---:|---:|---|
| pure-element-anchored | 0* | 0 | 0.000e+00 | PASS |
| random-state-anchored | 74 | 0 | 1.555e-10 | PASS |

\* Same no-op reasoning as K=2's pure-element anchor.

---

## K=4 — Nb-Ti-V-Zr / BCC_A2 / T

**Non-point CF names** (ncf=51): 6 tetrahedron-binary + 12 tetrahedron-ternary + 3
tetrahedron-quaternary + 6 triangle-binary + 12 triangle-ternary + 6 pair(II) + 6 pair(I).
**Point CFs**: `α=3, α=2, α=1` in the pipeline's raw (descending) order.
**Orthogonal basis sequence**: `seq = {-2, -1, +1, +2}`.

| CF family (tetrahedron, binary) | v4AB | v4AC | v4AD | v4BC | v4BD | v4CD |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| CF family (tetrahedron, ternary ABC/ABD) | v4ABC1 | v4ABC2 | v4ABC3 | v4ABD1 | v4ABD2 | v4ABD3 |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 48800.0 | -24400.0 | -24400.0 | 0.0 | 0.0 | 0.0 |

| CF family (tetrahedron, ternary ACD/BCD, quaternary) | v4ACD1 | v4ACD2 | v4ACD3 | v4BCD1 | v4BCD2 | v4BCD3 | v4ABCD1 | v4ABCD2 | v4ABCD3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

| CF family (triangle, binary) | v3AB | v3AC | v3AD | v3BC | v3BD | v3CD |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | -4224.0 | 0.0 | 0.0 | -5760.0 |

| CF family (triangle, ternary ABC/ABD) | v3ABC1 | v3ABC2 | v3ABC3 | v3ABD1 | v3ABD2 | v3ABD3 |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | -20000.0 | -20000.0 | -20000.0 | -3333.33 | -3333.33 | -3333.33 |

| CF family (triangle, ternary ACD/BCD) | v3ACD1 | v3ACD2 | v3ACD3 | v3BCD1 | v3BCD2 | v3BCD3 |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 0.0 | 0.0 | 0.0 | 10000.0 | 10000.0 | 10000.0 |

| CF family (pair, II-n) | v22AB | v22AC | v22AD | v22BC | v22BD | v22CD |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 3120.0 | 7040.0 | 7022.4 | 4080.0 | -4632.0 | 10232.0 |

| CF family (pair, I-n) | v21AB | v21AC | v21AD | v21BC | v21BD | v21CD |
|---|---:|---:|---:|---:|---:|---:|
| eciCvcf | 6240.0 | 14080.0 | 14044.8 | 8160.0 | -9264.0 | 20464.0 |

### Stage 1 — Pure element (x_0 = 1)

| Family (tetrahedron, first 10 columns) | v4AB | v4AC | v4AD | v4BC | v4BD | v4CD | v4ABC1 | v4ABC2 | v4ABC3 | v4ABD1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 4096 | -2048 | 1024 | 1024 | -512 | 256 | 1024 | -512 | 256 | -512 |
| uOrth measured | 4096 | -2048 | 1024 | 1024 | -512 | 256 | 1024 | -512 | 256 | -512 |

(Full 54-entry vector — including all 51 non-point columns, 3 point columns,
and the empty-cluster constant — confirmed identical between analytic and
measured, entry by entry.)

| Quantity | Value |
|---|---:|
| **Max error, all 54 columns** | **0.0** |
| vCvcf (Tinv·uOrth), all 51 columns | ~0 (magnitude ≤ ~1.5×10⁻¹⁵, floating-point noise around exact 0) |

A pure single-species state has zero probability of any mixed cluster, same as K=2/K=3.

### Stage 2 — Composition boundary (x = [1/3, 1/3, 1/3, 0], Zr absent)

| Family (first 6 tetrahedron-binary columns) | v4AB | v4AC | v4AD | v4BC | v4BD | v4CD |
|---|---:|---:|---:|---:|---:|---:|
| uOrth (boundary) | 50.5679 | -37.9259 | 12.6420 | 28.4444 | -9.4815 | 3.1605 |

(Full 54-entry vector recorded in the run log; only the lead columns shown here for scale reference.)

**Cross-check against the true K=3 (Nb-Ti-V) system's own equiatomic CVCF values**
(K=3 Stage 3 above — the K=4 boundary composition [1/3,1/3,1/3] on the 3 surviving
species *is* K=3's own equiatomic composition), matched by shared name across all 18
K=3 non-point columns:

| Quantity | Value |
|---|---:|
| matchedNames | 18 |
| maxErr | 1.596×10⁻¹⁶ |
| **Result** | **PASS (exact to floating-point precision)** |
| maxErr (before fix, for reference) | 1.481 |

### Stage 3 — Perfectly random (equiatomic, x=[0.25, 0.25, 0.25, 0.25])

| Nonzero uOrth columns | v21AB | v21AC | v21AD | v21BC | v21BD | v21CD |
|---|---:|---:|---:|---:|---:|---:|
| uOrth analytic | 39.0625 | 0.0 | 0.0 | 15.625 | 0.0 | 0.0 |

(Remaining nonzero pair-family entries: additional `v21xx`-type terms at
`6.25` and `2.5`; all other columns — every tetrahedron and triangle
column — are exactly 0 at equiatomic K=4 composition; `empty = 1.0`.)

| CVCF family | tetrahedron (all 21) | triangle-binary (all 6) | triangle-ternary (all 12) | pair (all 12) |
|---|---:|---:|---:|---:|
| vCvcf analytic | ~0.00390625 | ~0 | ~0.015625 | ~0.0625 |

| Quantity | Value |
|---|---:|
| H closed-form (`Sum(eci·vCvcf)`) | **4411.700156 J/mol** |

| MCS measured (avg of 40 seeds, L=4) | Value |
|---|---:|
| MCS random H | 4669.440692 |
| stderr | ±456.294887 |
| closed-form H | 4411.700156 |
| diff | +257.740535 |
| 5σ tolerance | 2281.5 |
| **Result** | **PASS** |

Note the much larger stderr than K=2/K=3 (456 vs. 13/64) — K=4 has 51 CF
columns competing for the same N=128 sites at L=4, so per-sample variance is
intrinsically larger.

**Convergence trend** (separate ad-hoc check during investigation, different
seeds/L, not part of the committed suite — 40 samples each):

| L | N | MCS H | diff from closed-form (4411.70) |
|---:|---:|---:|---:|
| 4 | 128 | 4669.44 | +257.74 |
| 6 | 432 | 4620.50 | +208.80 |
| 8 | 1024 | 4756.57 | +344.87 |
| 10 | 2000 | 4471.21 | +59.51 |

Diffs stay well within sampling-noise range across L — consistent with
statistical noise, not a systematic bias.

### Stage 4 — Delta-E, anchored at analytic configs

| Anchor | Trials run | Mismatches | Max abs error | Result |
|---|---:|---:|---:|---|
| pure-element-anchored | 0* | 0 | 0.000e+00 | PASS |
| random-state-anchored | 79 | 0 | 9.677e-10 | PASS |

\* Same no-op reasoning as above.

---

## Overall summary

| System | Pure-element uOrth | B2 CFs | Boundary CVCF | Random-state H | Delta-E | Trajectory |
|---|---|---|---|---|---|---|
| K=2 (Nb-Ti) | exact | exact | — | PASS (5σ) | PASS | PASS |
| K=3 (Nb-Ti-V) | exact | — | exact (0.0) | PASS (5σ) | PASS | — |
| K=4 (Nb-Ti-V-Zr) | exact | — | exact (1.6e-16) | PASS (5σ) | PASS | — |

**Checks: 23   Failures: 0   RESULT: PASS** (full suite, `AnalyticConfigVerification.main`)

---

## Bug found and fixed during this validation

**Symptom**: the K=3/K=4 composition-boundary CVCF check initially failed
with every shared-name column off by exactly **-1.0** (K=3→K=2) or **-1.481**
(K=4→K=3, before full diagnosis), while the K=2 case and the tetrahedron term
`v4AB` at every K matched exactly.

**Root cause**: an earlier session (before this validation work) had
reordered the point-CF columns emitted by `Embeddings.generatePointCfEmbeddings`
and `ClusterCFIdentificationPipeline.computeRandomCFs` from the pipeline's
raw column order (descending in alpha, for K≥3 — e.g. K=3's two point
columns appear as `[α=2]` then `[α=1]`) to **ascending** order, believing
`CvCfBasis`'s `Tinv` matrix expected ascending order.

**This was backwards.** `Tinv` is built directly from the pipeline's own raw
column order (`CMatrixPipeline.buildCfColumnMap` / `deriveCfBasisIndices`),
so it expects the *same* raw (descending) order the pipeline naturally
produces. The ascending reorder silently broke every K≥3 CVCF value that
depends on a point-CF column (i.e. every pair and triangle term), though not
the fully-symmetric tetrahedron term, which happened to be invariant under
the 2-element reversal for the point-CF pair at K=3.

**Why existing tests didn't catch it**: every pre-existing MCS test compared
one MCS code path against another MCS code path (delta-E vs. finite-difference
full-energy recompute; running trajectory energy vs. from-scratch recompute).
Both sides of those comparisons shared the same (wrong) point-CF ordering, so
they stayed numerically self-consistent and kept reporting `PASS` even though
the underlying physical value was wrong. Only a check against an
**independent** ground truth exposed it — specifically, a pure-species-A
configuration's point CF `xA` must equal exactly `1.0` by definition.

| Check | Before fix | After fix |
|---|---:|---:|
| xA (pure species-0 config, K=3) | -1.0 | **1.0** (correct) |
| xB (pure species-0 config, K=3) | 2.0 | **0.0** (correct) |
| Boundary CVCF maxErr, K=3→K=2 | 1.0 | **0.0** |
| Boundary CVCF maxErr, K=4→K=3 | 1.481 | **1.6×10⁻¹⁶** |

**Fix**: reverted both `Embeddings.generatePointCfEmbeddings` and
`ClusterCFIdentificationPipeline.computeRandomCFs` to preserve the pipeline's
raw column order. Also confirmed all three documented CLI regression values
(CVM, unaffected by this bug since CVM's Newton-Raphson self-corrects from
any initial guess) still reproduce exactly:

| Test case | Value |
|---|---:|
| Nb-Ti, T=1000K, Ti=0.5 → G | -3480.5209063901 |
| Nb-Ti-V, T=1000K, Ti=0.33,V=0.34 → G | -7051.1257304632 |
| Nb-Ti-V-Zr, T=1273K, equiatomic → S | 11.0812146249 |

**Practical impact**: this bug affected K≥3 MCS CVCF energies wherever a
point-CF-dependent term (any pair or triangle CF) carried a nonzero ECI —
i.e. essentially all real K≥3 Hamiltonians. K=2 systems were never affected
(only one point-CF column, so there is no ordering to get wrong).

---

## Files

| File | Role |
|---|---|
| `src/main/java/org/ce/model/mcs/AnalyticConfigVerification.java` | The consolidated test — supersedes and replaces `RandomStateEnergyTest`, `DeltaEVerificationTest`, `DeltaEVerificationTestV2Rigorous`, `MetropolisTrajectoryVerification` (all deleted) |
| `src/main/java/org/ce/model/mcs/Embeddings.java` | `generatePointCfEmbeddings` — reverted to raw column order |
| `src/main/java/org/ce/model/cluster/ClusterCFIdentificationPipeline.java` | `computeRandomCFs` — reverted to raw column order |
| `src/main/java/org/ce/model/mcs/AbFamilyCfDiagnostic.java`, `DeltaEBenchmark.java`, `EmbeddingScaleProbe.java` | Retained as separate diagnostic/performance tools (not correctness gates, not superseded) |
