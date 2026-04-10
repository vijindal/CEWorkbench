# CEWorkbench Three-Layer Restructuring Plan

## Context

The user wants a clear three-layer architecture:
- **UI Layer** — collects system config (elements, structure, model type) and calculation params, passes to Model layer
- **Model Layer** — a persistent object (singleton per session) that holds resolved cluster data + Hamiltonian; created once per system config, reused across all calculations
- **Calculation Layer** — uses the persistent model; results stream to UI real-time or batch

**Core problem being solved:** Both `CVMGibbsModel.fromThermodynamicInput()` and `MCSEngine.compute()` currently call `ClusterIdentificationWorkflow.identify()` on **every single calculation**, even repeated T-scans on the same system. Cluster identification is expensive and should only run when the system (elements, structure, model) actually changes.

## Approach: Introduce a Persistent Model Registry

Surgical changes only — no rewrites. 4 new files, 3 modified files, 1 one-liner addition.

---

## New Classes (4 files)

### 1. `org/ce/model/SystemKey.java` — Cache key
Immutable value object. Identifies a unique system configuration.
```java
public final class SystemKey {
    public final String structurePhase;  // from cec.structurePhase
    public final String model;           // from cec.model
    public final int numComponents;
    public final String resolvedHamiltonianId; // after CVCF resolution
    // equals() + hashCode() on all 4 fields
}
```

### 2. `org/ce/model/CVMModel.java` — Persistent CVM model
Created once; reused for all CVM calculations on the same system.
```java
public final class CVMModel {
    public final SystemKey key;
    public final AllClusterData clusterData;
    public final CVMGibbsModel gibbsModel;   // built at creation, not per-calc
    public final CvCfBasis basis;
    public final CECEntry cec;

    public static CVMModel create(SystemKey key, AllClusterData clusterData, CECEntry cec) {
        CvCfBasis basis = CvCfBasis.Registry.INSTANCE.get(
                cec.structurePhase, cec.model, numComponents);
        CVMGibbsModel gibbsModel = new CVMGibbsModel(
                clusterData.getDisorderedClusterResult(),
                clusterData.getDisorderedCFResult(),
                clusterData.getCMatrixResult(),
                basis);
        return new CVMModel(key, clusterData, gibbsModel, basis, cec);
    }
}
```

### 3. `org/ce/model/MCSModel.java` — Persistent MCS model
```java
public final class MCSModel {
    public final SystemKey key;
    public final AllClusterData clusterData;
    public final CvCfBasis basis;
    public final CMatrix.Result cmatResult;  // pre-extracted
    public final CECEntry cec;

    public static MCSModel create(SystemKey key, AllClusterData clusterData, CECEntry cec) { ... }
}
```

### 4. `org/ce/model/ActiveModelRegistry.java` — Lifecycle manager

The heart of the change. Holds the active `CVMModel` and `MCSModel`; rebuilds them only when `SystemKey` changes.

```java
public class ActiveModelRegistry {
    private volatile SystemKey  activeKey;
    private volatile CVMModel   cvmModel;
    private volatile MCSModel   mcsModel;
    private final CECManagementWorkflow cecWorkflow;

    public synchronized CVMModel getCVMModel(
            String resolvedHamiltonianId, CECEntry cec,
            int numComponents, Consumer<String> progressSink) {
        SystemKey key = new SystemKey(cec.structurePhase, cec.model,
                                      numComponents, resolvedHamiltonianId);
        if (key.equals(activeKey) && cvmModel != null) return cvmModel;  // CACHE HIT
        // CACHE MISS: run cluster identification once
        AllClusterData data = ClusterIdentificationWorkflow.identify(
                buildRequest(cec, numComponents), progressSink);
        cvmModel = CVMModel.create(key, data, cec);
        mcsModel = null;   // invalidate MCS (different model, same cluster data)
        activeKey = key;
        return cvmModel;
    }

    public synchronized MCSModel getMCSModel(
            String resolvedHamiltonianId, CECEntry cec,
            int numComponents, Consumer<String> progressSink) {
        SystemKey key = new SystemKey(cec.structurePhase, cec.model,
                                      numComponents, resolvedHamiltonianId);
        if (key.equals(activeKey) && mcsModel != null) return mcsModel;  // CACHE HIT
        // If activeKey matches but mcsModel is null (cvmModel was built first),
        // reuse the already-identified clusterData from cvmModel
        AllClusterData data = (key.equals(activeKey) && cvmModel != null)
                ? cvmModel.clusterData
                : ClusterIdentificationWorkflow.identify(buildRequest(cec, numComponents), progressSink);
        mcsModel = MCSModel.create(key, data, cec);
        activeKey = key;
        return mcsModel;
    }

    public synchronized void invalidate() {
        activeKey = null; cvmModel = null; mcsModel = null;
    }
}
```

---

## Modified Classes (3 files + 1 one-liner)

### 5. `CVMEngine.java` — Add model-aware overload

Add alongside existing `compute(ThermodynamicInput)` (keep it unchanged):
```java
public EquilibriumState compute(CVMModel model, double temperature, double[] composition,
                                 Consumer<String> progressSink, Consumer<ProgressEvent> eventSink) {
    // Same pipeline, but uses model.gibbsModel, model.basis, model.cec directly
    // No cluster ID call. No getBasis() call.
    validateInputs(temperature, composition);
    double[] eci = CECEvaluator.evaluate(model.cec, temperature, model.basis, "CVM");
    CVMSolver.EquilibriumResult result = new CVMSolver().minimize(
            model.gibbsModel, composition, temperature, eci, 1.0e-5, progressSink, eventSink);
    validateConvergence(result, progressSink);
    return buildEquilibriumState(temperature, composition, result);
}
```

### 6. `MCSEngine.java` — Add model-aware overload

Add alongside existing `compute(ThermodynamicInput)` (keep it unchanged):
```java
public EquilibriumState compute(MCSModel model, double temperature, double[] composition,
                                 int mcsL, int mcsNEquil, int mcsNAvg,
                                 Consumer<String> progressSink, Consumer<ProgressEvent> eventSink) {
    // Same pipeline, but uses model.clusterData, model.basis, model.cec directly
    // No ClusterIdentificationWorkflow.identify() call.
}
```

### 7. `ThermodynamicWorkflow.java` — Route through registry

Add `ActiveModelRegistry` as constructor parameter (old constructor kept for tests):
```java
public ThermodynamicWorkflow(CECManagementWorkflow cecWorkflow,
                              ThermodynamicEngine cvmEngine,
                              ThermodynamicEngine mcsEngine,
                              ActiveModelRegistry modelRegistry) { ... }
```

In `runCalculation()`, after loading `CECEntry cec`, replace the current `engine.compute(input)` call:
```java
int numComponents = request.composition.length;
EquilibriumState state;
if ("CVM".equals(request.engineType)) {
    CVMModel model = modelRegistry.getCVMModel(
            resolvedHamiltonianId, cec, numComponents, request.progressSink);
    state = cvmEngine.compute(model, request.temperature, request.composition,
                              request.progressSink, request.eventSink);
} else {
    MCSModel model = modelRegistry.getMCSModel(
            resolvedHamiltonianId, cec, numComponents, request.progressSink);
    state = mcsEngine.compute(model, request.temperature, request.composition,
                              request.mcsL, request.mcsNEquil, request.mcsNAvg,
                              request.progressSink, request.eventSink);
}
```
Note: `cvmEngine` and `mcsEngine` here cast to `CVMEngine`/`MCSEngine` to access the new overloads, or the interface is left as-is and the new overloads are called directly.

### 8. `CEWorkbenchContext.java` — Wire registry

```java
// Add field:
private final ActiveModelRegistry modelRegistry;

// In constructor:
this.modelRegistry = new ActiveModelRegistry(cecWorkflow);
ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
        cecWorkflow, cvmEngine, mcsEngine, modelRegistry);  // pass registry

// Add accessor:
public ActiveModelRegistry getModelRegistry() { return modelRegistry; }
```

### 9. `CalculationPanel.java` — Invalidate on system change (one-liner)

Inside the `WorkbenchContext` change listener already registered in `CalculationPanel`:
```java
context.addChangeListener(() -> {
    appCtx.getModelRegistry().invalidate();   // ADD THIS LINE
    if (context.hasSystem()) { ... }          // existing code
});
```

---

## Files NOT Changing
- `ThermodynamicEngine.java` (interface) — untouched
- `ThermodynamicInput.java` — untouched  
- `ThermodynamicRequest.java` — untouched  
- `CalculationService.java` — untouched  
- `LineScanWorkflow.java`, `GridScanWorkflow.java`, `FiniteSizeScanWorkflow.java` — untouched (caching is transparent)
- `WorkbenchContext.java` — untouched
- `ClusterIdentificationWorkflow.java` — untouched (still called by registry)
- `Main.java` (CLI) — untouched (gets caching for free via context)

---

## Model Lifecycle

| Event | Result |
|---|---|
| First calculation after startup | MISS → cluster ID runs → model stored |
| Same system, different T or x | HIT → model reused, no cluster ID |
| User switches system in GUI | `CalculationPanel` listener → `registry.invalidate()` → next calc rebuilds |
| Switch CVM→MCS same system | `getMCSModel()` reuses `cvmModel.clusterData`, skips cluster ID |
| T-scan (loop over T) | First point: cluster ID; all others: cache hit |
| CLI invocation | New `CEWorkbenchContext` → new registry → first calc triggers ID; subsequent scan points reuse |

---

## Critical Files

- `src/main/java/org/ce/CEWorkbenchContext.java`
- `src/main/java/org/ce/workflow/thermo/ThermodynamicWorkflow.java`
- `src/main/java/org/ce/domain/engine/cvm/CVMEngine.java`
- `src/main/java/org/ce/domain/engine/mcs/MCSEngine.java`
- `src/main/java/org/ce/ui/gui/CalculationPanel.java`

New files in new package `org/ce/model/`:
- `SystemKey.java`
- `CVMModel.java`
- `MCSModel.java`
- `ActiveModelRegistry.java`

---

## Verification

1. **Single-point CVM**: Run calculation twice on same system → second run should NOT print "Starting Always Fresh Structural Identification" in log
2. **T-scan**: Run temperature scan on same system → cluster ID message appears only once (first point)
3. **System change**: Change elements/structure in GUI → run calculation → cluster ID runs fresh
4. **CVM then MCS**: Set CVM, run, then switch to MCS same system → cluster ID NOT re-run (reuses data from CVMModel)
5. **CLI**: Run `calc_min` → run second `calc_min` in same session (if applicable) → cluster ID on first only
6. **Backward compatibility**: Old `ThermodynamicEngine.compute(ThermodynamicInput)` path still compiles and works for CLI `calc_fixed`
