# Plan: Move Model Building into Model Layer — MCSRunner accepts ModelSession

## Context

The user wants the calculation layer to pass model specifications to the model layer, which then
takes full responsibility for building the model (cluster algebra, CEC reading, initial config,
ECI evaluation and transform). Steps 2 and 3 (Hamiltonian evaluation, ECI transform) should
happen together with Step 4 (initialization) inside the model layer — not scattered across
ThermodynamicWorkflow and MCSRunner.

**Current problem:**
- `ThermodynamicWorkflow.runMcs()` does Step 2 (CECEvaluator.evaluate → eci[]) and validates
  C-matrix before handing off to MCSRunner — model-layer concerns leaked into the calculation layer
- `MCSRunner.Builder` requires the caller to decompose `ModelSession` into 6 raw fields
  (clusterData, eci, basis, matrixData, lcf, numComp) — the session was already available
- `MCSRunner.run()` does Step 3 (ECI transform via Tinv) — fine, stays there

**Target:**
- `ThermodynamicWorkflow.runMcs()` becomes a thin dispatcher — passes `ModelSession` directly
- `MCSRunner.Builder` gains `session(ModelSession, double temperature)` — extracts all
  model fields internally, performs CECEvaluator.evaluate(), validates C-matrix, builds ECIs
- Steps 2, 3, 4 all visible and sequenced inside the model layer

---

## What Each Layer Does After the Change

```
UI Layer
  → ModelSpecifications (elements, structure, model, engine)
  → CalculationSpecifications (T, x, L, nEquil, nAvg)
  → CalculationService.execute()

Calculation Layer (ThermodynamicWorkflow.runMcs)
  → retrieves pre-built ModelSession                (unchanged)
  → wires progress callbacks                        (unchanged)
  → calls MCSRunner.builder().session(session, T)   ← NEW: passes session directly
  → calls builder.build().run()                     (unchanged)
  → runs MCSStatisticsProcessor on MCResult         (unchanged)

Model Layer (MCSRunner.Builder)
  Step 2: CECEvaluator.evaluate(session.cecEntry, T, session.cvcfBasis) → eci (CVCF basis)
  Step 2b: validate C-matrix dimensions            ← moved from ThermodynamicWorkflow
  Step 3: ECI transform via Tinv (stays in MCSRunner.run() as today)
  Step 4: build lattice, embeddings, LatticeConfig (stays in MCSRunner.run() as today)
```

---

## Precise Changes

### 1. `MCSRunner.Builder` — add `session()` method, remove `eci()` requirement

**New import in MCSRunner.java:**
```java
import org.ce.model.ModelSession;
import org.ce.model.hamiltonian.CECEvaluator;
import org.ce.model.PhysicsConstants;
```

**New field in Builder:**
```java
private ModelSession session = null;   // set by session() convenience method
```

**New Builder method:**
```java
public Builder session(ModelSession s, double temperature) {
    this.session = s;
    this.T       = temperature;
    return this;
}
```

**Modified `build()` in Builder** — add block that auto-fills from session if provided:
```java
public MCSRunner build() {
    if (session != null) {
        // Step 2: Evaluate Hamiltonian — ECI in CVCF basis at temperature T
        var matData = session.clusterData.getMatrixData();
        int cmatCols = /* extract from matData */;
        if (cmatCols != session.cvcfBasis.totalCfs())
            throw new IllegalStateException("C-matrix dimension mismatch ...");
        this.eci       = CECEvaluator.evaluate(session.cecEntry, T, session.cvcfBasis, "MCS");
        // Fill remaining fields from session
        this.clusterData = session.clusterData.getDisorderedClusterResult().getDisClusterData();
        this.numComp     = session.numComponents();
        this.basis       = session.cvcfBasis;
        this.matrixData  = matData;
        this.lcf         = session.clusterData.getLcf();
        if (R <= 0) this.R = PhysicsConstants.R_GAS;
    }
    // existing validation follows unchanged
    if (clusterData == null) throw new IllegalStateException("clusterData required");
    if (eci == null)         throw new IllegalStateException("eci required");
    if (T <= 0)              throw new IllegalStateException("T must be > 0");
    ...
}
```

### 2. `ThermodynamicWorkflow.runMcs()` — remove Steps 2 & 3, use `session()` builder method

**Remove from runMcs():**
- `CECEvaluator.evaluate(...)` call (Step 2 — moves to MCSRunner.Builder.build())
- C-matrix dimension validation block (Step 2b — moves to MCSRunner.Builder.build())
- All 13 individual `.clusterData()`, `.eci()`, `.numComp()`, `.basis()`, `.matrixData()`, `.lcf()` builder calls

**Replace with:**
```java
MCSRunner.Builder builder = MCSRunner.builder()
    .session(session, request.temperature)   // model building now in model layer
    .composition(request.composition)
    .nEquil(nEquil)
    .nAvg(nAvg)
    .L(L)
    .seed(System.currentTimeMillis())
    .cancellationCheck(Thread.currentThread()::isInterrupted);
```

**Remove import from ThermodynamicWorkflow.java:**
```java
// remove:  import org.ce.model.hamiltonian.CECEvaluator;
// remove:  import org.ce.model.cluster.ClusterCFIdentificationPipeline.ClusCoordListData;
```

### 3. No other file changes

- `ModelSession` — unchanged (already holds everything needed)
- `MetropolisMC` — unchanged (already receives transformed ECIs)
- `CalculationService` — unchanged
- `ThermodynamicRequest` — unchanged (T still passed from calculation layer — that's correct)
- `LocalEnergyCalcTest` — unchanged

---

## Step Visibility After Change

### Inside `MCSRunner.Builder.build()`:
```
Step 2: CECEvaluator.evaluate(session.cecEntry, T, basis) → eci (CVCF basis)
Step 2b: validate C-matrix dims  
```

### Inside `MCSRunner.run()` (unchanged body, now fully self-contained):
```
Step 3: buildEciByOrbitType(eci_cvcf, tc, basis) → eci_orth (orthogonal basis)
Step 4: buildBCCPositions(L) → positions
        Embeddings.generate(positions, clusterData, L) → emb
        LatticeConfig(N, numComp); config.randomise(xFrac, rng)
        → MetropolisMC(emb, eci_orth, ...) engine
        → engine.run(config, sampler) → MCResult
```

---

## Critical Files

- [MCSRunner.java](src/main/java/org/ce/model/mcs/MCSRunner.java) — Builder gains `session()` + `build()` fills fields from session
- [ThermodynamicWorkflow.java](src/main/java/org/ce/calculation/workflow/thermo/ThermodynamicWorkflow.java) — `runMcs()` simplified to ~10 builder lines
- [ModelSession.java](src/main/java/org/ce/model/ModelSession.java) — read-only, no changes
- [MetropolisMC.java](src/main/java/org/ce/model/mcs/MetropolisMC.java) — no changes
- [LocalEnergyCalcTest.java](src/test/java/org/ce/model/mcs/LocalEnergyCalcTest.java) — no changes

---

## Verification

```bash
./gradlew compileJava compileTestJava   # must BUILD SUCCESSFUL
./gradlew test                          # LocalEnergyCalcTest must pass (3 tests)
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 0.5 --verbose"
```
