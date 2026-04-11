# Calculation Methods: Concrete Examples

This document shows how to implement specific thermodynamic calculations using the calculation layer pattern.

---

## Example 1: Spinodal Line (Phase Stability Limit)

**Physics:** Find compositions where d²G/dx² = 0 (limit of metastability).

**Input:** Temperature range, engine (CVM/MCS)

**Output:** List of (T, x_spinodal) points

### Request DTO

```java
package org.ce.calculation.workflow.thermo;

public class SpinodolRequest {
    public final double temperatureStart;
    public final double temperatureEnd;
    public final double temperatureStep;
    public final double compositionMin;
    public final double compositionMax;
    
    public final Consumer<String> progressSink;
    public final Consumer<ProgressEvent> eventSink;
    
    public final int mcsL;
    public final int mcsNEquil;
    public final int mcsNAvg;
    
    public SpinodolRequest(
            double tStart, double tEnd, double tStep,
            double xMin, double xMax,
            Consumer<String> progressSink,
            Consumer<ProgressEvent> eventSink,
            int mcsL, int mcsNEquil, int mcsNAvg) {
        this.temperatureStart = tStart;
        this.temperatureEnd = tEnd;
        this.temperatureStep = tStep;
        this.compositionMin = xMin;
        this.compositionMax = xMax;
        this.progressSink = progressSink;
        this.eventSink = eventSink;
        this.mcsL = mcsL;
        this.mcsNEquil = mcsNEquil;
        this.mcsNAvg = mcsNAvg;
    }
}

// Result DTO
public class SpinodolPoint {
    public final double temperature;
    public final double xAlpha;  // lower phase composition
    public final double xBeta;   // upper phase composition
    
    public SpinodolPoint(double T, double xAlpha, double xBeta) {
        this.temperature = T;
        this.xAlpha = xAlpha;
        this.xBeta = xBeta;
    }
}
```

### Workflow Implementation

```java
package org.ce.calculation.workflow.thermo;

public class SpinodolWorkflow {
    
    private final ThermodynamicWorkflow thermoWorkflow;
    private static final Logger LOG = Logger.getLogger(SpinodolWorkflow.class.getName());
    
    public SpinodolWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }
    
    /**
     * Finds spinodal points by scanning temperature and searching for 
     * composition where d²G/dx² = 0.
     */
    public List<SpinodolPoint> findSpinodal(
            ModelSession session,
            SpinodolRequest request) throws Exception {
        
        LOG.info("SpinodolWorkflow.findSpinodal — ENTER");
        List<SpinodolPoint> spinodal = new ArrayList<>();
        
        int pointCount = 0;
        for (double T = request.temperatureStart; 
             T <= request.temperatureEnd; T += request.temperatureStep) {
            
            if (request.progressSink != null) {
                request.progressSink.accept(String.format(
                    "Spinodal: T=%.1f K — searching composition range [%.3f, %.3f]",
                    T, request.compositionMin, request.compositionMax));
            }
            
            // Binary search for both spinodal points (alpha and beta phases)
            double xAlpha = binarySearchGxx(session, T, 
                request.compositionMin, 0.5, request);
            double xBeta = binarySearchGxx(session, T, 
                0.5, request.compositionMax, request);
            
            if (!Double.isNaN(xAlpha) && !Double.isNaN(xBeta)) {
                spinodal.add(new SpinodolPoint(T, xAlpha, xBeta));
                pointCount++;
                
                if (request.progressSink != null) {
                    request.progressSink.accept(String.format(
                        "  ✓ Spinodal point: T=%.1f  xα=%.4f  xβ=%.4f",
                        T, xAlpha, xBeta));
                }
            } else {
                if (request.progressSink != null) {
                    request.progressSink.accept(String.format(
                        "  ✗ No spinodal found at T=%.1f K", T));
                }
            }
        }
        
        LOG.info("SpinodolWorkflow.findSpinodal — EXIT: " + pointCount + " points");
        return spinodal;
    }
    
    /**
     * Binary search for composition in [xMin, xMax] where d²G/dx² = 0.
     * Returns NaN if sign never changes (no spinodal in interval).
     */
    private double binarySearchGxx(ModelSession session, double T, 
                                    double xMin, double xMax,
                                    SpinodolRequest request) {
        final double TOL = 1e-8;
        final int MAX_ITER = 50;
        
        double gxx_min = evaluateGxx(session, T, xMin, request);
        double gxx_max = evaluateGxx(session, T, xMax, request);
        
        // Check if sign change exists
        if (gxx_min * gxx_max > 0) {
            return Double.NaN;  // No zero crossing
        }
        
        double xLow = xMin, xHigh = xMax;
        
        for (int iter = 0; iter < MAX_ITER; iter++) {
            if (xHigh - xLow < TOL) {
                return (xLow + xHigh) / 2.0;
            }
            
            double xMid = (xLow + xHigh) / 2.0;
            double gxx_mid = evaluateGxx(session, T, xMid, request);
            
            if (gxx_mid == 0) return xMid;
            
            if (gxx_mid * gxx_min < 0) {
                xHigh = xMid;
            } else {
                xLow = xMid;
            }
        }
        
        return (xLow + xHigh) / 2.0;
    }
    
    /**
     * Compute d²G/dx² by finite difference.
     * Model layer evaluates G; calculation layer computes derivative.
     */
    private double evaluateGxx(ModelSession session, double T, double x,
                                SpinodolRequest request) {
        final double DX = 1e-5;
        
        try {
            double[] comp_lo  = {1 - (x - DX), x - DX};
            double[] comp_mid = {1 - x, x};
            double[] comp_hi  = {1 - (x + DX), x + DX};
            
            ThermodynamicRequest req_lo = new ThermodynamicRequest(
                T, comp_lo, null, null, request.mcsL, request.mcsNEquil, request.mcsNAvg);
            ThermodynamicRequest req_mid = new ThermodynamicRequest(
                T, comp_mid, null, null, request.mcsL, request.mcsNEquil, request.mcsNAvg);
            ThermodynamicRequest req_hi = new ThermodynamicRequest(
                T, comp_hi, null, null, request.mcsL, request.mcsNEquil, request.mcsNAvg);
            
            ThermodynamicResult res_lo  = thermoWorkflow.runCalculation(session, req_lo);
            ThermodynamicResult res_mid = thermoWorkflow.runCalculation(session, req_mid);
            ThermodynamicResult res_hi  = thermoWorkflow.runCalculation(session, req_hi);
            
            // d²G/dx² ≈ (G(x+dx) - 2*G(x) + G(x-dx)) / dx²
            double gxx = (res_hi.gibbsEnergy - 2 * res_mid.gibbsEnergy + res_lo.gibbsEnergy)
                       / (DX * DX);
            
            return gxx;
        } catch (Exception e) {
            LOG.warning("evaluateGxx failed at T=" + T + ", x=" + x + ": " + e.getMessage());
            return Double.NaN;
        }
    }
}
```

### Service Registration

```java
public class CalculationService {
    
    private final ThermodynamicWorkflow thermoWorkflow;
    private final SpinodolWorkflow spinodolWorkflow;
    // ... other workflows ...
    
    public CalculationService(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
        this.spinodolWorkflow = new SpinodolWorkflow(thermoWorkflow);
        // ...
    }
    
    public List<SpinodolPoint> findSpinodal(
            ModelSession session,
            SpinodolRequest request) throws Exception {
        return spinodolWorkflow.findSpinodal(session, request);
    }
}
```

---

## Example 2: Equilibrium Phase Diagram

**Physics:** Find coexisting phase compositions at each T by solving G(xα) = G(xβ) and dG/dx|α = dG/dx|β (equal chemical potentials).

**Input:** Temperature range, engine (CVM only)

**Output:** List of (T, xα, xβ) tie lines

### Request & Result DTOs

```java
public class EquilibriumPhaseDiagramRequest {
    public final double temperatureStart;
    public final double temperatureEnd;
    public final double temperatureStep;
    
    public final Consumer<String> progressSink;
    public final Consumer<ProgressEvent> eventSink;
    
    public EquilibriumPhaseDiagramRequest(double tStart, double tEnd, double tStep,
                                          Consumer<String> progressSink,
                                          Consumer<ProgressEvent> eventSink) {
        this.temperatureStart = tStart;
        this.temperatureEnd = tEnd;
        this.temperatureStep = tStep;
        this.progressSink = progressSink;
        this.eventSink = eventSink;
    }
}

public class TieLine {
    public final double temperature;
    public final double xAlpha;      // equilibrium composition α phase
    public final double xBeta;       // equilibrium composition β phase
    public final double gibbsAlpha;  // G of α phase
    public final double gibbsBeta;   // G of β phase
    
    public TieLine(double T, double xA, double xB, double gA, double gB) {
        this.temperature = T;
        this.xAlpha = xA;
        this.xBeta = xB;
        this.gibbsAlpha = gA;
        this.gibbsBeta = gB;
    }
}
```

### Workflow Implementation

```java
public class EquilibriumPhaseDiagramWorkflow {
    
    private final ThermodynamicWorkflow thermoWorkflow;
    private static final Logger LOG = Logger.getLogger(EquilibriumPhaseDiagramWorkflow.class.getName());
    
    public EquilibriumPhaseDiagramWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }
    
    /**
     * Finds equilibrium tie lines by solving:
     *   μ_A(α) = μ_A(β)   ⟹ dG/dx|α = dG/dx|β
     *   G(α) = G(β)        ⟹ lever rule
     * 
     * CVM only (requires Gibbs free energy and derivatives).
     */
    public List<TieLine> findEquilibriumCurve(
            ModelSession session,
            EquilibriumPhaseDiagramRequest request) throws Exception {
        
        if (!session.engineConfig.isCvm()) {
            throw new IllegalArgumentException("Equilibrium phase diagram requires CVM");
        }
        
        LOG.info("EquilibriumPhaseDiagramWorkflow.findEquilibriumCurve — ENTER");
        List<TieLine> diagram = new ArrayList<>();
        
        for (double T = request.temperatureStart;
             T <= request.temperatureEnd; T += request.temperatureStep) {
            
            if (request.progressSink != null) {
                request.progressSink.accept(String.format("Phase diagram: T=%.1f K", T));
            }
            
            // Solve for tie line at this T
            TieLine tieLine = solveTieLine(session, T, request);
            if (tieLine != null) {
                diagram.add(tieLine);
                if (request.progressSink != null) {
                    request.progressSink.accept(String.format(
                        "  ✓ Tie line: xα=%.4f  xβ=%.4f  ΔG=%+.6f J/mol",
                        tieLine.xAlpha, tieLine.xBeta,
                        tieLine.gibbsBeta - tieLine.gibbsAlpha));
                }
            }
        }
        
        LOG.info("EquilibriumPhaseDiagramWorkflow.findEquilibriumCurve — EXIT: " 
                + diagram.size() + " tie lines");
        return diagram;
    }
    
    /**
     * Newton-Raphson solver for tie line at fixed T.
     * Finds (xα, xβ) such that:
     *   f1(xα, xβ) = dG/dx|xα - dG/dx|xβ = 0
     *   f2(xα, xβ) = G(xα) - G(xβ) = 0
     */
    private TieLine solveTieLine(ModelSession session, double T,
                                  EquilibriumPhaseDiagramRequest request) {
        final double TOL = 1e-8;
        final int MAX_ITER = 20;
        
        // Initial guess: spinodal points
        double xAlpha = 0.2;  // Could start from spinodal
        double xBeta = 0.8;
        
        for (int iter = 0; iter < MAX_ITER; iter++) {
            // Evaluate residuals
            double dGdx_alpha = evaluateDGDx(session, T, xAlpha, request);
            double dGdx_beta = evaluateDGDx(session, T, xBeta, request);
            ThermodynamicResult res_a = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, new double[]{1-xAlpha, xAlpha}));
            ThermodynamicResult res_b = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, new double[]{1-xBeta, xBeta}));
            
            double f1 = dGdx_alpha - dGdx_beta;
            double f2 = res_a.gibbsEnergy - res_b.gibbsEnergy;
            
            // Check convergence
            if (Math.abs(f1) < TOL && Math.abs(f2) < TOL) {
                return new TieLine(T, xAlpha, xBeta, 
                                   res_a.gibbsEnergy, res_b.gibbsEnergy);
            }
            
            // Newton step (simplified jacobian)
            double dfdx_alpha = evaluateD2GDx2(session, T, xAlpha, request);
            double dfdx_beta = evaluateD2GDx2(session, T, xBeta, request);
            
            double denom = dfdx_alpha - dfdx_beta;
            if (Math.abs(denom) < 1e-12) {
                return null;  // Jacobian singular
            }
            
            // Simple Newton step (could be improved to 2×2 Jacobian)
            double step = f1 / denom;
            xAlpha -= step * 0.5;
            xBeta += step * 0.5;
            
            // Bounds check
            xAlpha = Math.max(0.01, Math.min(0.5, xAlpha));
            xBeta = Math.max(0.5, Math.min(0.99, xBeta));
        }
        
        return null;  // Did not converge
    }
    
    private double evaluateDGDx(ModelSession session, double T, double x,
                                 EquilibriumPhaseDiagramRequest request) {
        final double DX = 1e-5;
        
        try {
            double[] comp_lo = {1 - (x - DX), x - DX};
            double[] comp_hi = {1 - (x + DX), x + DX};
            
            ThermodynamicResult res_lo = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, comp_lo));
            ThermodynamicResult res_hi = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, comp_hi));
            
            // dG/dx ≈ (G(x+dx) - G(x-dx)) / 2dx
            return (res_hi.gibbsEnergy - res_lo.gibbsEnergy) / (2 * DX);
        } catch (Exception e) {
            return Double.NaN;
        }
    }
    
    private double evaluateD2GDx2(ModelSession session, double T, double x,
                                   EquilibriumPhaseDiagramRequest request) {
        final double DX = 1e-5;
        
        try {
            double[] comp_lo = {1 - (x - DX), x - DX};
            double[] comp_mid = {1 - x, x};
            double[] comp_hi = {1 - (x + DX), x + DX};
            
            ThermodynamicResult res_lo = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, comp_lo));
            ThermodynamicResult res_mid = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, comp_mid));
            ThermodynamicResult res_hi = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, comp_hi));
            
            return (res_hi.gibbsEnergy - 2 * res_mid.gibbsEnergy + res_lo.gibbsEnergy)
                   / (DX * DX);
        } catch (Exception e) {
            return Double.NaN;
        }
    }
}
```

---

## Example 3: Heat Capacity Temperature Scan

**Physics:** C_v(T) from entropy: C_v = T · (∂²S/∂T²)

**Input:** Temperature range, composition

**Output:** Temperature vs. heat capacity

### Implementation Pattern

```java
public class HeatCapacityWorkflow {
    
    private final ThermodynamicWorkflow thermoWorkflow;
    
    public HeatCapacityWorkflow(ThermodynamicWorkflow thermoWorkflow) {
        this.thermoWorkflow = thermoWorkflow;
    }
    
    /**
     * Scans C_v(T) from entropy second derivative.
     * Model layer computes S(T); calculation layer computes dS/dT.
     */
    public List<HeatCapacityPoint> scanHeatCapacity(
            ModelSession session,
            double[] composition,
            double tStart,
            double tEnd,
            double tStep) throws Exception {
        
        List<HeatCapacityPoint> results = new ArrayList<>();
        
        for (double T = tStart; T <= tEnd; T += tStep) {
            final double DT = Math.max(1.0, T * 0.001);  // ~0.1% T
            
            // Three temperatures for numerical d²S/dT²
            ThermodynamicResult res_lo = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T - DT, composition));
            ThermodynamicResult res_mid = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T, composition));
            ThermodynamicResult res_hi = thermoWorkflow.runCalculation(
                session, new ThermodynamicRequest(T + DT, composition));
            
            // C_v = T · d²S/dT² = T · (S(T+dT) - 2*S(T) + S(T-dT)) / dT²
            double d2SdT2 = (res_hi.entropy - 2 * res_mid.entropy + res_lo.entropy) 
                          / (DT * DT);
            double cv = T * d2SdT2;
            
            results.add(new HeatCapacityPoint(T, cv));
        }
        
        return results;
    }
}

public class HeatCapacityPoint {
    public final double temperature;
    public final double heatCapacity;
    
    public HeatCapacityPoint(double T, double cv) {
        this.temperature = T;
        this.heatCapacity = cv;
    }
}
```

---

## Summary: Pattern for New Calculations

1. **Request DTO** — parameters specific to your calculation
2. **Result DTO** — output data structure
3. **Workflow class** — contains the algorithm
   - Loops/searches over parameters
   - Calls `thermoWorkflow.runCalculation()` repeatedly
   - Post-processes results (derivatives, searches, aggregates)
4. **Register in CalculationService** — expose as public entry point

The model layer is untouched — it always evaluates G/H/S given (T, x, ECI). The calculation layer decides **what to do with those evaluations**.
