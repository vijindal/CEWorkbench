# CVCF Basis Derivation — Status

Tracks progress on `CvCfBasis.java` registrations (Stage 4 of cluster/CF identification) across structures and component counts. See `CLAUDE.md` for the pipeline architecture and session contract.

## Current Status: ACTIVE — BCC_B2 binary complete, ternary in progress

## Verified and working (regression-tested, 0.00e+00 self-test diff)

| Structure | Model | K | Status |
|---|---|---|---|
| BCC_A2 | T | 2, 3, 4 | Verified |
| FCC_A1 | T | 2, 3, 4 | Verified |
| HCP_A3 | T | 2 (binary) | Verified — rebuilt from Jindal & Lele (2025) CALPHAD paper, Appendix 2, eq 59-64 |
| HCP_A3 | T | 3 (ternary) | Verified — rebuilt from paper eq 65+ with corrected site numbering (math site→our site: 1→1, 2→4, 3→2, 4→3, 5→5, 6→6); tetrahedron ABC-mixed terms are literal Mathematica translations, not paper-inferred |
| BCC_B2 | T | 2 (binary) | Verified — derived by hand (see below), rank 9/9, exact 0.00e+00 diff, η=0 at equiatomic |

## Blocked / not started

| Structure | Model | K | Status |
|---|---|---|---|
| HCP_A3 | T | 4 (quaternary) | Blocked — registration exists but M-matrix is singular (rank 30/31... actually pipeline reports tcf=85, needs 86-wide basis, but our 84-CV registration is short by 1). Table 19 says 83+1=84 should be right; real pipeline disagrees. Unresolved discrepancy, same status as before the HCP rework. |
| BCC_B2 | T | 3 (ternary) | In progress — see below |
| BCC_B2 | T | 4 (quaternary) | Not started |

## BCC_B2 binary — derivation summary

B2 has two sublattices: α = {p1, p4}, β = {p2, p3} (same 4-site tetrahedron as A2, `BCC_A2-T.txt`'s coordinates, just partitioned per `BCC_B2-T.txt`'s `{{p1,p4}},{{p2,p3}}` grouping).

Orthogonal-CF orbit structure (tcf=8, from pipeline Stage 2b):
- tetrahedron (t=0): unsplit, ncv=9
- triangle (t=1): splits into {1,2,3} and {1,2,4}
- II-n pair (t=2, mh=4 in A2): unsplit as **one** pipeline orbit, but α-α ({1,4}) and β-β ({2,3}) are physically distinct once sublattice occupation matters — picked as **two** independent CVCF variables anyway (see below)
- I-n pair (t=3, mh=3 in A2): two orbit *entries* in Stage 2 output, but same physical pair type — **one** CVCF variable, unchanged from A2
- point (t=4): splits into α and β

Final 9-variable CVCF basis:

| name | formula | ECI shared with A2's... |
|---|---|---|
| V4AB | `p1[A]p2[A]p3[B]p4[B]` | new — self-paired tetrahedron orbit (AABB≡ABAB≡BABA≡BBAA under (1↔4)(2↔3); **no antisymmetric partner exists**, must use plain product) |
| V31AB | `p1[A]p2[B]p3[B] − p1[B]p2[A]p3[A]` | v3AB |
| V32AB | `p2[A]p1[B]p4[B] − p2[B]p1[A]p4[A]` | v3AB (same ECI as V31AB) |
| V221AB | `p1[A]p4[B]` | v22AB |
| V222AB | `p2[A]p3[B]` | v22AB (same ECI as V221AB) |
| V21AB | `p1[A]p2[B]` | v21AB (unchanged) |
| xA, xB | `p1[A]`, `p1[B]` | — |
| eta | `p1[A] − p2[A]` (LRO parameter, no ECI) | — |

**Pitfall found the hard way:** the antisymmetric combination `p1[A]p4[B] − p1[B]p4[A]` for the II-n pair evaluates to an *exact zero row* in the real M-matrix — sites 1 and 4 (and separately 2, 3) have no AB/BA distinction at that orbit under B2's true space group. Only the **plain product** works. Verify any new antisymmetric pick against the real M-matrix (rank check), not just plausibility — see `AGENT_NOTES` pattern below.

## Shared-code fixes (apply to all structures, not just B2)

1. **`CvCfBasis.Definition`** gained an explicit `numPointCfs` field (default = `numComponents`, matching prior behavior for A2/FCC/HCP). B2 registers `numPointCfs=3` (xA, xB, eta) since it has one more point-like quantity than `numComponents`. Use the `register(..., numPointCfs)` overload for any future ordered structure with order parameters.
2. **`ClusterCFIdentificationPipeline.computeRandomCFs`** previously hardcoded `pointCfCount = K-1`, silently leaving split point-orbit columns at 0 instead of the correct random-state value (this produced a spurious η=0.5 at equiatomic instead of the correct η=0). Fixed to use the tracked `nxcf` field. This affects the "random CF at composition" self-test for any ordered structure with more than `K-1` point orbits.
3. **`ClusterCFIdentificationPipeline.runFullWorkflow`** now checks `CvCfBasis.isSupported(structurePhase, ...)` directly before falling back to `resolveParentStructure(structurePhase)` — previously an ordered structure with its own registered CVCF definition (like the new `BCC_B2_T_2`) was silently resolved to its disordered parent's definition instead.

## Next session: BCC_B2 ternary (K=3)

Confirmed via pipeline: tcf=35 (needs 36-variable basis). Per-cluster-type CV counts (from `classified[t][j][k]` Stage-3 breakdown, counting `(j,k)` entries):

| Cluster type | # CVCF variables |
|---|---|
| Tetrahedron | 9 |
| Triangle | 12 |
| II-n pair (t=2, A2 mh=4) | 4 |
| I-n pair (t=3, A2 mh=3) | 6 |
| Point (incl. order parameters) | 4 |
| **Total** | 35 |

Point CFs confirmed via literature research (Sanchez–Ducastelle–Gratias generalized cluster description): **xA, xB, xC, ηB, ηC** where ηB = xB(α)−xB(β), ηC = xC(α)−xC(β), and ηA = −(ηB+ηC) is dependent (traceless constraint). This matches the pipeline's nxcf=4 (independent point DOF) + xA as a conventional redundant 5th.

**Open question, not yet resolved:** does the ternary I-n pair (t=3) split into 6 CVCF variables (2 per species-pair AB/AC/BC), or does it stay unsplit like binary's I-n pair (which needed only 1 variable, `V21AB`, despite Stage 2 showing two `j`-orbit entries)? The `classified[3][j][k]` data shows k=0,1 for all three `j` (species-pairs), suggesting a split — but binary's analogous `classified[3][0][k=0,1]` *also* showed a k-split and yet the correct answer (verified by rank + self-test) was **no split needed**. This "classified k-split ≠ CVCF split" trap means per-cluster-type CV counts from `classified[]` should be treated as an upper bound / sanity check, not a direct prescription — always verify via M-matrix rank, same as the II-n antisymmetric-cancels-to-zero pitfall above.

Next step: propose I-n pair as 3 unsplit variables (AB, AC, BC, matching A2 exactly) first — cheaper to test — and only fall back to 6 if rank comes up short.

## Verification pattern (use for any new registration)

1. Get real orbit structure: `./gradlew run --args="type1a <elements> <structure> T --verbose"`, grep `t=X j=Y:` lines for `ncv`/`wcv`, and `[RESULT] lcv:` for orbit counts.
2. Propose CVCF picks, wire into `CvCfBasis.java` via `register(...)`.
3. Compile, run the same command — watch for "Matrix is singular" errors.
4. If singular: temporarily add M-matrix row dump (`for i: emit("MROW "+i+" ("+name+"): "+row)`gated behind an env var, e.g. `DEBUG_MMATRIX`), pull into Python/numpy, check `np.linalg.matrix_rank` and SVD null vector to find which CV/orbit direction is missing or redundant. Remove the debug block once resolved — don't leave it in.
5. Once rank is full: check `[SELF-TEST] CV VERIFICATION` output — every cluster type should show `Diff: 0.00e+00` (or ~1e-16 floating noise) at equiatomic composition. Any order parameter (like η) should read exactly 0 at equiatomic/random state.
6. Regression-check all previously-verified structures before calling it done — the shared-code fixes above (points #1-3) can silently affect unrelated structures.
