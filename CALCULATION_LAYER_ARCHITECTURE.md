# Calculation Layer Architecture

## Overview

The calculation layer sits between the **UI layer** (user inputs) and the **model layer** (physics optimizers). It orchestrates how thermodynamic calculations are performed and how results are delivered.

```
┌─────────────────────────────────────────────────────────────┐
│ UI Layer (org.ce.ui)                                        │
│  - CalculationPanel: user inputs (T, composition, MCS params)│
│  - LineScanPanel: parameter scan UI                         │
│  - MapPanel: 2D result visualization                         │
└────────────────────┬────────────────────────────────────────┘
                     │ request (.ThermodynamicRequest)
                     │ single point / scan / special
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ Calculation Layer (org.ce.calculation.workflow)             │
│  - CalculationService: public API router                    │
│  - ThermodynamicWorkflow: single-point CVM/MCS dispatcher   │
│  - ScanWorkflows: GridScan, LineScan, FiniteSizeScan       │
│  - MCSStatisticsProcessor: post-processing for MCS          │
│  - [NEW] Specialized methods for phase diagrams, etc.       │
└────────────────────┬────────────────────────────────────────┘
                     │ result (.ThermodynamicResult)
                     │ progress events (.ProgressEvent)
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ Model Layer (org.ce.model)                                  │
│  - CVMGibbsModel: physics evaluator                         │
│  - CVMSolver: Newton-Raphson optimizer                      │
│  - MCEngine + MCSampler + MCSRunner: MCS algorithm drivers  │
│  - LocalEnergyCalc: energy evaluator (stateless)            │
└─────────────────────────────────────────────────────────────┘
```

## Current Request/Response Contracts

### ThermodynamicRequest (Input)
The parameter bundle sent from UI to calculation layer:

```java
public class ThermodynamicRequest {
    // Calculation parameters (only)
    public final double temperature;           // Temperature in K
    public final double[] composition;         // Mole fractions [x_A, x_B, ...]
    
    // Progress sinks (optional, null allowed)
    public final Consumer<String> progressSink;        // Text output
    public final Consumer<ProgressEvent> eventSink;    // Chart events
    
    // MCS-specific parameters
    public final int mcsL;        // Lattice size (default 4)
    public final int mcsNEquil;   // Equilibration sweeps (default 1000)
    public final int mcsNAvg;     // Averaging sweeps (default 2000)
}
```

### ThermodynamicResult (Output)
The complete thermodynamic state returned to UI:

```java
public class ThermodynamicResult {
    public final double temperature;           // K (echoed from input)
    public final double[] composition;         // mole fractions (echoed)
    
    // Thermodynamic properties
    public final double gibbsEnergy;           // G, J/mol (CVM only)
    public final double enthalpy;              // H, J/mol
    public final double entropy;               // S, J/(mol·K) (CVM only)
    
    // Error bars (from MCS post-processing)
    public final double stdEnthalpy;           // σ(H), J/mol
    public final double heatCapacity;          // C_v, J/(mol·K)
    
    // Correlation functions
    public final double[] optimizedCFs;        // CVM equilibrium CFs (non-point)
    public final double[] avgCFs;              // MCS mean CFs
    public final double[] stdCFs;              // MCS CF standard errors
}
```

## Current Calculation Methods

### 1. Single-Point Calculation
**Entry:** `CalculationService.runSinglePoint(session, request)`

**Path:**
```
UI → CalculationService.runSinglePoint(session, ThermodynamicRequest)
   → ThermodynamicWorkflow.runCalculation(session, request)
   → switch on session.engineConfig.engineType:
       case "CVM" → CVMSolver.minimize()  → CVMEquilibriumState
       case "MCS" → MCSRunner.build().run() → MCResult
                  → MCSStatisticsProcessor.computeStatistics()
   → ThermodynamicResult
```

**Current UI:** CalculationPanel "Calculate" button (T, x at single point)

### 2. Finite-Size Scaling Scan (MCS production)
**Entry:** `CalculationService.runFiniteSizeScan(session, T, composition, nEquil, nAvg, sinks)`

**Path:**
```
UI → FiniteSizeScan.run(L=12, 16, 24)
   → for each L:
       ThermodynamicWorkflow.runCalculation() with mcsL=L
   → Extrapolate to thermodynamic limit (1/N → 0)
   → Final ThermodynamicResult at infinite size
```

**Current UI:** CalculationPanel "Production Run" button (MCS only)

### 3. Line Scan (1D Parameter Scan)
**Entry:** `CalculationService.runLineScanTemperature()` or `.runLineScanComposition()`

**Path:**
```
UI → LineScan.scanTemperature(T_start, T_end, T_step)
   → for each T:
       ThermodynamicWorkflow.runCalculation()
   → List<ThermodynamicResult>
```

**Current UI:** LineScanPanel (temperature or composition scan)

### 4. Grid Scan (2D Parameter Scan)
**Entry:** `CalculationService.runGridScan(tStart, tEnd, tStep, xStart, xEnd, xStep)`

**Path:**
```
UI → GridScan.scanTX()
   → for each (T, x):
       ThermodynamicWorkflow.runCalculation()
   → List<List<ThermodynamicResult>>
```

**Current UI:** MapPanel (T×x grid visualization)

---

## How to Add New Calculation Methods

### Pattern: Add a New Specialized Workflow

**Example: Phase diagram search (spinodal/miscibility gap)**

#### Step 1: Create Request DTO
If the new calculation needs different parameters, create a request class:

```java
package org.ce.calculation.workflow.thermo;

public class PhaseDiagramRequest {
    // Phase diagram search parameters
    public final double temperatureStart;
    public final double temperatureEnd;
    public final double[] compositionStart;
    
    public final Consumer<String> progressSink;
    public final Consumer<ProgressEvent> eventSink;
    
    // CVM/MCS params
    public final int mcsL;
    public final int mcsNEquil;
    public final int mcsNAvg;
    
    // Constructor ...
}
```

#### Step 2: Create Specialized Workflow Class (in same package)
The new workflow calls `ThermodynamicWorkflow` repeatedly and post-processes results:

```java
package org.ce.calculation.workflow.thermo;

public class PhaseDiagramWorkflow {
    
    private final ThermodynamicWorkflow thermoWorkflow;
    
    public PhaseDiagramWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }
    
    /**
     * Searches for spinodal line (limit of metastability).
     * 
     * Model layer contract: given (T, x, engine), returns G/H/S.
     * Calculation layer adds: searches for d²G/dx² = 0.
     */
    public List<SpinodialPoint> findSpinodal(
            ModelSession session,
            PhaseDiagramRequest request) throws Exception {
        
        List<SpinodialPoint> spinodal = new ArrayList<>();
        
        for (double T = request.temperatureStart; 
             T <= request.temperatureEnd; T += tStep) {
            
            // Binary scan: find composition where d²G/dx² = 0
            // Requires: sampling G(T, x) across composition range
            // → calls thermoWorkflow.runCalculation() multiple times per T
            
            double xSpinodal = binarySearchSpinodal(session, T, request);
            spinodal.add(new SpinodialPoint(T, xSpinodal));
        }
        
        return spinodal;
    }
    
    /**
     * Binary search for composition where d²G/dx² = 0 at fixed T.
     * (or equivalently, G_xx = 0)
     */
    private double binarySearchSpinodal(ModelSession session, double T, 
                                        PhaseDiagramRequest request) {
        double xLow = 0.0, xHigh = 1.0;
        
        // Evaluate G at endpoints
        double gxx_low = evaluateGxx(session, T, xLow, request);
        double gxx_high = evaluateGxx(session, T, xHigh, request);
        
        // Binary search for zero crossing
        while (xHigh - xLow > 1e-6) {
            double xMid = (xLow + xHigh) / 2.0;
            double gxx_mid = evaluateGxx(session, T, xMid, request);
            
            if (gxx_mid * gxx_low < 0) {
                xHigh = xMid; gxx_high = gxx_mid;
            } else {
                xLow = xMid; gxx_low = gxx_mid;
            }
        }
        
        return (xLow + xHigh) / 2.0;
    }
    
    /**
     * Evaluates d²G/dx² by finite difference (CVM only).
     * Model layer: calculates G(T, x) and derivatives.
     * Calculation layer: post-processes derivatives.
     */
    private double evaluateGxx(ModelSession session, double T, double x,
                               PhaseDiagramRequest request) {
        double dx = 1e-6;
        
        // Three calls to model layer
        double[] comp_lo = {1 - (x - dx), x - dx};
        double[] comp_mid = {1 - x, x};
        double[] comp_hi = {1 - (x + dx), x + dx};
        
        ThermodynamicResult g_lo  = thermoWorkflow.runCalculation(
            session, new ThermodynamicRequest(T, comp_lo, null, null, ...));
        ThermodynamicResult g_mid = thermoWorkflow.runCalculation(
            session, new ThermodynamicRequest(T, comp_mid, null, null, ...));
        ThermodynamicResult g_hi  = thermoWorkflow.runCalculation(
            session, new ThermodynamicRequest(T, comp_hi, null, null, ...));
        
        // d²G/dx² ≈ (G(x+dx) - 2*G(x) + G(x-dx)) / dx²
        return (g_hi.gibbsEnergy - 2*g_mid.gibbsEnergy + g_lo.gibbsEnergy) 
               / (dx * dx);
    }
}

// Result DTO
public class SpinodialPoint {
    public final double temperature;
    public final double composition;
    
    public SpinodialPoint(double temperature, double composition) {
        this.temperature = temperature;
        this.composition = composition;
    }
}
```

#### Step 3: Add to CalculationService
Register the new workflow in the service:

```java
public class CalculationService {
    
    private final ThermodynamicWorkflow  thermoWorkflow;
    private final PhaseDiagramWorkflow   phaseDiagramWorkflow;
    // ... other workflows ...
    
    public CalculationService(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow           = thermoWorkflow;
        this.phaseDiagramWorkflow     = new PhaseDiagramWorkflow(thermoWorkflow);
        // ...
    }
    
    /**
     * New public entry point for phase diagram calculations.
     */
    public List<SpinodialPoint> findSpinodal(
            ModelSession session,
            PhaseDiagramRequest request) throws Exception {
        return phaseDiagramWorkflow.findSpinodal(session, request);
    }
}
```

#### Step 4: Add UI (e.g., new panel or button)
In the UI layer, add a button/panel that calls the new service method:

```java
public class PhaseDiagramPanel extends JPanel {
    
    private final CalculationService service;
    
    private void runSpinodal() {
        ModelSession session = context.getActiveSession();
        
        PhaseDiagramRequest request = new PhaseDiagramRequest(
            tStart, tEnd, xStart,
            this::onProgress,
            this::onChartEvent,
            mcsL, mcsNEquil, mcsNAvg);
        
        worker = new SwingWorker<List<SpinodialPoint>, Object>() {
            @Override
            protected List<SpinodialPoint> doInBackground() {
                return service.findSpinodal(session, request);
            }
            // ... process results ...
        };
        worker.execute();
    }
}
```

---

## Design Rules

### ✅ DO
- **Each workflow class** owns a specific calculation strategy (scan type, search method, etc.)
- **Workflows call `ThermodynamicWorkflow`** to perform single-point calculations
- **Workflows post-process results** — compute derived quantities, search for features, aggregate data
- **Workflows are stateless** — they accept session + request, return results
- **Request/Result DTOs** hold parameters/output for each calculation type
- **Progress sinks** are optional (null-checked) — allow silent operation for tests/CLI

### ❌ DON'T
- **Workflows should NOT** directly construct model layer objects (CVMSolver, MCEngine, etc.)
  - Always call through `ThermodynamicWorkflow` (which is the dispatcher)
- **Workflows should NOT** re-validate cluster data, Hamiltonian, or basis
  - That is done once by `ModelSession.Builder`
- **Workflows should NOT** assume CVM or MCS
  - Let `ThermodynamicWorkflow` dispatch based on `session.engineConfig.engineType`
- **Don't add calculation logic to result DTOs**
  - `ThermodynamicResult` is immutable data only; post-processing happens in workflows

### Separation of Concerns

| Layer | Responsibility | Example |
|-------|-----------------|---------|
| **Model** | Physics evaluation + optimization | CVMGibbsModel evaluates G/H/S; CVMSolver minimizes |
| **Calculation** | Orchestration + post-processing | LineScan calls ThermodynamicWorkflow N times; aggregates results |
| **UI** | User inputs + visualization | CalculationPanel reads T, x; sends request; displays result |

---

## Reference Implementation

The existing `ScanWorkflows` is the simplest pattern to follow:

```
ScanWorkflows.LineScan
├── scanTemperature() loops T
│   └── for each T: thermoWorkflow.runCalculation() → ThermodynamicResult
├── Result: List<ThermodynamicResult>
└── UI consumes list → charts line plot

ScanWorkflows.GridScan
├── scanTX() nested loop T, x
│   └── for each (T,x): thermoWorkflow.runCalculation()
├── Result: List<List<ThermodynamicResult>>
└── UI consumes grid → heatmap
```

For a phase diagram search, follow this exact pattern: loop the parameter(s) of interest, call `thermoWorkflow.runCalculation()` repeatedly, post-process results (search for features, compute derivatives), return structured output.

---

## Adding Model-Layer Capabilities

If a new calculation method requires **new physics** (not just re-using CVM/MCS), add it to the model layer **first**, then wrap it in a calculation workflow.

Example: adding a new optimization algorithm
1. Implement `NewSolver` in `model/cvm/` or `model/mcs/`
2. Have it return clean contract objects (state + metadata)
3. Update `ThermodynamicWorkflow` to dispatch to it
4. Wrap in calculation workflows as needed

The principle: **physics goes in model, orchestration goes in calculation, UI goes in ui.**
