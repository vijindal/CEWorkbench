# Architecture Review: Clean Layer Contracts

**Date:** 2026-04-11  
**Status:** ✅ CONTRACTS VERIFIED AND CLEAN

---

## Executive Summary

The codebase now has **clean, testable contracts at layer boundaries**:

- **Model Layer Contract:** Given (T, x, ECI), return equilibrium state with thermodynamic values
- **Calculation Layer Contract:** Given model outputs, compute scientific quantities (error bars, derived properties) the user requested
- **UI Layer Contract:** Present results and accept user input

All layer boundaries are unidirectional: `model ← calculation ← ui`. No upward dependencies.

---

## Model Layer Contract

### Physics Evaluators (Stateless)

#### CVMGibbsModel.evaluate()

```java
ModelResult evaluate(
    double[] u,                  // correlation functions (free variables)
    double[] moleFractions,       // composition
    double temperature,           // K
    double[] eci                  // ECI parameters (a + b·T)
) → ModelResult {
    G, H, S,                     // free energy, enthalpy, entropy
    Gu, Guu,                     // gradient, Hessian of G
    Hu, Hu,                      // derivatives of H, S
    Suu
}
```

**Contract:** Pure function. No side effects. Given physical parameters, compute thermodynamic quantities and their derivatives.

**Testable:** ✅
```java
double[] u = {0.5, 0.0};
CVMGibbsModel.ModelResult result = model.evaluate(u, x, T, eci);
assert result.G < result.H * T;  // G = H - TS at equilibrium
```

---

#### MCSampler.sample() / getSeriesXxx()

```java
void sample(
    LatticeConfig config,
    EmbeddingData emb,
    double currentEnergy
)

List<Double> getSeriesHmix()     // Raw Hmix time series (J)
List<Double> getSeriesE()         // Raw E time series (J)
List<Double>[] getSeriesCF()      // Raw CF time series
```

**Contract:** Accumulate raw time series. No statistics. Pure accumulation.

**Testable:** ✅
```java
MCSampler sampler = new MCSampler(...);
for (sweep in sweeps) sampler.sample(config, emb, E);
List<Double> series = sampler.getSeriesHmix();
assert series.size() == sweeps;
```

---

### Algorithm Drivers (Own Loops & Convergence)

#### CVMSolver.minimize()

```java
EquilibriumResult minimize(
    CVMGibbsModel model,
    double[] moleFractions,
    double temperature,
    double[] eci,
    double tolerance,
    Consumer<String> progressSink,      // optional logging
    Consumer<ProgressEvent> eventSink   // optional events
) → EquilibriumResult {
    u[],                    // optimized correlation functions
    ModelResult,            // G/H/S at equilibrium
    iterations,             // convergence count
    finalGradientNorm,      // ||∇G||
    converged,              // success flag
    trace                   // iteration snapshots
}
```

**Contract:** Newton-Raphson minimization. Input: model, conditions, tolerance. Output: equilibrium state + convergence info.

**Testable:** ✅
```java
CVMSolver solver = new CVMSolver();
CVMSolver.EquilibriumResult result = solver.minimize(
    model, x, T, eci, 1e-5);
assert result.converged;
assert result.finalGradientNorm < 1e-5;
assert result.trace.size() == result.iterations;
```

---

#### MCEngine.run()

```java
MCResult run(
    LatticeConfig config,
    MCSampler sampler
) → MCResult {
    temperature, composition,
    acceptRate,
    meanHmixPerSite,                // ⟨Hmix⟩/N
    meanEnergyPerSite,              // ⟨E⟩/N
    meanCFs[],                      // average correlation functions
    heatCapacity,                   // crude ⟨E²⟩-⟨E⟩² estimate
    // Statistics are NaN (post-processing handled in calculation layer)
    tauInt = NaN,
    stdEnergyPerSite = NaN,
    stdHmixPerSite = NaN,
    stdCFs[] = null,
    cvJackknife = NaN,
    cvStdErr = NaN
}
```

**Contract:** Run equilibration then averaging sweeps. Accumulate means in sampler. Return result with **only raw averages**. Statistics post-processing is not model responsibility.

**Testable:** ✅
```java
MCEngine engine = new MCEngine(...);
MCResult result = engine.run(config, sampler);
assert result.getHmixPerSite() != Double.NaN;
assert result.getTauInt() == Double.NaN;  // model doesn't compute
List<Double> rawSeries = sampler.getSeriesHmix();
assert rawSeries.size() == nAvg;
```

---

## Calculation Layer Contract

### Orchestrators

#### ThermodynamicWorkflow.runCvm()

```
Input:  ModelSession + ThermodynamicRequest
        (session: cluster data, Hamiltonian, CVCF basis)
        (request: T, x, progress sinks)

Output: ThermodynamicResult {
    T, x,
    gibbsEnergy, enthalpy, entropy,
    optimizedCFs[] at equilibrium,
    (stdEnthalpy, heatCapacity are NaN)
}
```

**Steps:**
1. Validate inputs
2. Build CVMGibbsModel from session
3. Evaluate ECI at T
4. **Call model-layer optimizer:** CVMSolver.minimize() → EquilibriumResult
5. Extract and return ThermodynamicResult

**No post-processing.** CVM returns exact values (G, H, S); no statistics needed.

**Testable:** ✅
```java
ThermodynamicWorkflow wf = new ThermodynamicWorkflow();
ThermodynamicResult result = wf.runCalculation(session, request);
assert result.temperature == request.temperature;
assert result.gibbsEnergy != Double.NaN;
assert result.enthalpy != Double.NaN;
```

---

#### ThermodynamicWorkflow.runMcs()

```
Input:  ModelSession + ThermodynamicRequest
        (session: cluster data, Hamiltonian, CVCF basis)
        (request: T, x, L, nEquil, nAvg, progress sinks)

Output: ThermodynamicResult {
    T, x,
    gibbsEnergy = NaN,              // not available in canonical ensemble
    enthalpy = ⟨H⟩/site,             // with error bar stdHmixPerSite
    entropy = NaN,
    heatCapacity = Cv from jackknife, // with error bar cvStdErr
    avgCFs[] with stdCFs[]
}
```

**Steps:**
1. Validate C-matrix dimensions
2. Evaluate ECI at T
3. **Call model-layer optimizer:** MCSRunner.run() → MCSRunResult
4. **NEW: Call calculation-layer statistics:** MCSStatisticsProcessor
5. Rebuild MCResult with computed statistics
6. Return ThermodynamicResult

**Key:** Separates simulation (model) from statistics (calculation).

**Testable:** ✅
```java
ThermodynamicWorkflow wf = new ThermodynamicWorkflow();
ThermodynamicResult result = wf.runCalculation(session, request);
assert result.temperature == request.temperature;
assert result.enthalpy != Double.NaN;
assert result.getCvJackknife() != Double.NaN;
assert result.getCvStdErr() > 0;  // has error bar
```

---

### Statistical Post-Processor

#### MCSStatisticsProcessor

```java
MCSStatisticsProcessor(
    N, R, T,
    List<Double> seriesHmix,       // raw time series from MCSampler
    List<Double> seriesE,
    List<Double>[] seriesCF
)

void computeStatistics()

// Returns:
getTauInt()                        // integrated autocorrelation time (Sokal)
getStatInefficiency()              // s = 1 + 2·τ_int
getNEff()                          // effective independent samples
getBlockSizeUsed()                 // auto-determined
getNBlocks()
getStdEnergyPerSite()              // SEM of ⟨E⟩/N
getStdHmixPerSite()                // SEM of ⟨H⟩/N
getStdCFs()[]                      // SEM of each CVCF CF
getCvJackknife()                   // unbiased Cv from leave-one-out
getCvStdErr()                      // jackknife SEM of Cv
```

**Contract:** Given raw time series and physical parameters (N, R, T), compute publication-quality statistics.

**Responsibilities:**
- Sokal automatic windowing for τ_int
- Block averaging with proper SEM
- Jackknife resampling for Cv
- No physics; pure statistics

**Testable:** ✅
```java
List<Double> series = Arrays.asList(1.0, 1.1, 0.9, 1.2, ...);
MCSStatisticsProcessor proc = new MCSStatisticsProcessor(N, R, T, series, series, cfs);
proc.computeStatistics();
assert proc.getTauInt() >= 0.5;
assert proc.getStdHmixPerSite() > 0;
assert proc.getNBlocks() >= 2;
```

---

## UI Layer Contract

### CalculationPanel

```
Input:  User selects:
        - System (elements, structure, model)
        - Engine (CVM or MCS)
        - Temperature, composition
        - MCS parameters (L, nEquil, nAvg)

Output: Display ThermodynamicResult to OutputPanel
        - G/H/S (CVM) or H/Cv (MCS)
        - Error bars (MCS only)
        - Correlation functions
```

**Responsibilities:**
- Collect user inputs
- Validate system identity
- Trigger session rebuild
- Call CalculationService.runCalculation()
- Stream progress to OutputPanel

**Does NOT:** Perform physics or statistics; pure UI.

---

## Dependency Flow

```
┌─ Model Layer ─────────────────────────────────────────┐
│  Physics Evaluators      Algorithm Drivers            │
│  • CVMGibbsModel.evaluate()                           │
│  • LocalEnergyCalc (static)           • CVMSolver     │
│  • SiteOperatorBasis                  • MCEngine      │
│  • CvCfEvaluator                      • MCSRunner     │
│  • EmbeddingData                      • MCSampler     │
│                                                        │
│  ProgressEvent (produced here)                        │
└────────────────────────┬──────────────────────────────┘
                         ↑
                         │ (calls optimizer)
                         │
┌─ Calculation Layer ────┼──────────────────────────────┐
│  Orchestrators & Statistics Processors                │
│  • ThermodynamicWorkflow                              │
│  • MCSStatisticsProcessor                             │
│  • CalculationService                                 │
│                                                        │
│  Consumes ProgressEvent for feedback                  │
└────────────────────────┬──────────────────────────────┘
                         ↑
                         │ (returns ThermodynamicResult)
                         │
┌─ UI Layer ────────────┼───────────────────────────────┐
│  Input & Display                                       │
│  • CalculationPanel                                   │
│  • OutputPanel (consumes ProgressEvent)               │
│  • MainWindow, WorkbenchContext                       │
└────────────────────────────────────────────────────────┘
```

**Direction:** Unidirectional downward. `ui` → `calculation` → `model`. Never reverse.

---

## Contract Verification Checklist

### Model Layer ✅

- [x] CVMSolver: accepts (model, T, x, ECI) → EquilibriumResult
- [x] MCEngine: accepts (config, sampler) → MCResult (means only)
- [x] MCSampler: accumulates raw time series, no statistics
- [x] No imports from calculation or ui layers (except ProgressEvent which it produces)
- [x] All evaluators are stateless functions
- [x] All optimizers own their loops and convergence
- [x] ProgressEvent moved to model layer (produced here, consumed elsewhere)

### Calculation Layer ✅

- [x] ThermodynamicWorkflow: receives (ModelSession, ThermodynamicRequest)
- [x] Dispatches to model-layer optimizers
- [x] MCSStatisticsProcessor: receives raw time series, computes statistics
- [x] Rebuilds MCResult with computed statistics
- [x] No upward dependencies to ui
- [x] Clear separation: model simulates, calculation analyzes

### UI Layer ✅

- [x] CalculationPanel: collects input, triggers workflows
- [x] Consumes ProgressEvent from model via calculation layer
- [x] No physics or statistics computation
- [x] Depends only on calculation layer API

---

## Testing Implications

### Model Layer Tests (Pure Unit Tests)

```java
@Test
void testCVMGibbsModel() {
    CVMGibbsModel model = ...;
    ModelResult result = model.evaluate(u, x, T, eci);
    assertEquals(...);  // No mocks; pure math
}

@Test
void testCVMSolverConvergence() {
    CVMSolver solver = new CVMSolver();
    EquilibriumResult result = solver.minimize(model, x, T, eci, 1e-5);
    assertTrue(result.converged);
}

@Test
void testMCSamplerAccumulation() {
    MCSampler sampler = new MCSampler(...);
    sampler.sample(config, emb, E);
    sampler.sample(config, emb, E);
    assertEquals(2, sampler.getSampleCount());
}
```

**Characteristics:**
- No mocks (model operates on concrete data)
- No external services
- Deterministic inputs → deterministic outputs
- Can be run in CI/CD in seconds

### Calculation Layer Tests (Integration Tests)

```java
@Test
void testMCSStatisticsProcessor() {
    List<Double> series = generateTimeSeries(...);
    MCSStatisticsProcessor proc = new MCSStatisticsProcessor(N, R, T, series, ...);
    proc.computeStatistics();
    assertTrue(proc.getTauInt() > 0);
    assertTrue(proc.getNBlocks() >= 2);
}

@Test
void testThermodynamicWorkflowCVM() {
    ModelSession session = buildMockSession();
    ThermodynamicRequest req = new ThermodynamicRequest(1000, x, null, null);
    ThermodynamicResult result = workflow.runCvm(session, req);
    assertNotNull(result.gibbsEnergy);
}
```

**Characteristics:**
- May mock ModelSession (data; not behavior)
- Can mock progress sinks
- Tests orchestration, not physics
- Fast; no long simulations

### UI Layer Tests (System Tests)

```java
// GUI tests would verify:
// - Button clicks trigger workflow
// - Results display correctly
// - Progress updates appear
```

---

## Summary: Clean Contracts Achieved

| Aspect | Status | Evidence |
|--------|--------|----------|
| Model → no upward deps | ✅ | No calculation/ui imports except ProgressEvent (which model produces) |
| Calculation → model only | ✅ | Imports and calls model-layer optimizers directly |
| UI → calculation only | ✅ | Calls CalculationService, receives ThermodynamicResult |
| Physics vs. Statistics separated | ✅ | Model returns NaN for statistics; Calculation computes them |
| Optimizers own loops | ✅ | CVMSolver owns NR loop; MCEngine owns sweep loop |
| Evaluators are stateless | ✅ | CVMGibbsModel.evaluate() is pure function |
| Time series accumulation only | ✅ | MCSampler.sample() + getSeriesXxx() |
| Clean input → output | ✅ | All public methods have clear signatures |
| Testable contracts | ✅ | All layers can be unit-tested independently |

**Result:** Architecture is **clean, maintainable, and testable**. Layer responsibilities are explicit. Contracts are unambiguous.

---

## Related Commits

- `9035e85` — Separate MCS statistical post-processing from model layer simulation
- `0ffb04f` — Simplify calculation layer: eliminate redundant ThermodynamicInput intermediary
- `368be39` — Resolve layer violation: move ProgressEvent from calculation to model layer (this review)

