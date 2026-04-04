# CVM Calculation Tracing Guide

This document describes the logging added to trace CVM thermodynamic calculations through the GUI.

## Quick Start: Enable Logging

To see detailed CVM traces, configure Java logging in your IDE or application:

### Option 1: System Property (at JVM startup)

```bash
java -Djava.util.logging.config.file=logging.properties ...
```

### Option 2: logging.properties Configuration

Create `logging.properties` in project root or classpath:

```properties
# Global level
.level = INFO

# CVM-specific loggers
org.ce.domain.engine.cvm.CVMEngine.level = INFO
org.ce.domain.engine.cvm.CVMPhaseModel.level = INFO
org.ce.domain.engine.cvm.NewtonRaphsonSolverSimple.level = FINE
org.ce.domain.cluster.ClusterVariableEvaluator.level = FINER
org.ce.domain.engine.cvm.CVMFreeEnergy.level = FINER

# Console handler
handlers = java.util.logging.ConsoleHandler
java.util.logging.ConsoleHandler.level = FINEST
java.util.logging.ConsoleHandler.formatter = java.util.logging.SimpleFormatter
java.util.logging.SimpleFormatter.format = [%4$-7s] %5$s%n
```

---

## Logging Flow: Entry to Exit

### 1. **CVMEngine.compute()** — Main entry point

**File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMEngine.java`

**Log Messages:**
```
=== CVMEngine.compute() START ===
Input parameters:
  systemId: <id>
  systemName: <name>
  temperature: <T> K
  composition: [<x1>, <x2>, ...]
  numComponents: <K>

Validating C-Matrix...
✓ C-Matrix found

Validating composition...
✓ Composition valid (sum=1.000000000)

Validating temperature...
✓ Temperature valid

Loading C-Matrix from AllClusterData...
✓ C-Matrix loaded

Creating CVMInput with topology...
✓ CVMInput created

Evaluating ECI at T=<T> K...
  Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
  CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
  Processing <N> CEC terms...
    eci[<idx>] (<cfName>) = <a> + <b>*<T> = <value>
  ECI array: [<values>]
✓ ECI evaluated (4 non-point terms)

Running CVM N-R minimization...
✓ CVM minimization converged in <N> iterations (||Gu||=<value>)
  G(eq) = <G> J/mol
  H(eq) = <H> J/mol
  S(eq) = <S> J/(mol·K)

Emitting N-R iteration trace to eventSink...
✓ Emitted <N> iteration snapshots

=== CVMEngine.compute() SUCCESS ===
```

**Key Points:**
- Validates all inputs (C-matrix, composition, temperature)
- Shows ECI transformation at temperature (equation: `eci = a + b*T`)
- Final thermodynamic values at equilibrium

---

### 2. **CVMPhaseModel.ensureMinimized()** — Lazy minimization trigger

**File:** `app/src/main/java/org/ce/domain/engine/cvm/CVMPhaseModel.java`

**Log Messages:**
```
CVMPhaseModel.ensureMinimized() — STARTING minimization
  systemId: <id>
  temperature: <T> K
  composition: [<x1>, <x2>, ...]
  ncf: 4, tcf: 6
  tolerance: 1.0e-6
  ECI: [[0]=<val>, [3]=<val>, ...]

  Calling NewtonRaphsonSolverSimple.solve()...
  Newton-Raphson solver returned
  Cached <N> iteration snapshots

  Solver converged!
  Equilibrium CFs (non-point): [<v1>, <v2>, <v3>, <v4>]
  Evaluating CVMFreeEnergy at equilibrium CFs...
  Thermodynamic values: G=<G>, H=<H>, S=<S>

CVMPhaseModel.minimize — SUCCESS: T=<T> K, x=[...], iterations=<N>, ||dG||=<norm>, G=<G> J/mol, elapsed=<ms> ms

CVMPhaseModel.ensureMinimized() — COMPLETE
```

**Key Points:**
- Shows parameters at start of minimization
- Reports equilibrium CF values and final thermodynamic values
- Elapsed time in milliseconds

---

### 3. **NewtonRaphsonSolverSimple.minimize()** — N-R iteration loop

**File:** `app/src/main/java/org/ce/domain/engine/cvm/NewtonRaphsonSolverSimple.java`

**Log Messages (FINE level):**
```
NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1000.0, xB=0.5, tolerance=1.0e-6
  Initial CF (random state): [v1, v2, v3, v4]

  [iteration loop continues, see FINEST below]

NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=42, ||Gu||=9.2e-7
```

**Log Messages (FINEST level — every 20 iterations or when close to convergence):**
```
NewtonRaphsonSolverSimple — iter   1: G=<G> H=<H> S=<S> ||Gu||=<norm>
NewtonRaphsonSolverSimple — iter  20: G=<G> H=<H> S=<S> ||Gu||=<norm>
NewtonRaphsonSolverSimple — iter  40: G=<G> H=<H> S=<S> ||Gu||=<norm>
```

**Error Cases:**
```
NewtonRaphsonSolverSimple — SINGULAR HESSIAN at iter=<N>: <reason>

NewtonRaphsonSolverSimple — NOT CONVERGED after <maxIter> iterations (||Gu||=<norm> > tolerance=<tol>)
  Final state: G=<G>, H=<H>, S=<S>
```

**Key Points:**
- Initial random-state CFs shown
- Convergence tracked via gradient norm ||Gu||
- Thermodynamic values updated every iteration
- Singular Hessian failures detected and reported

---

## Example: Full Trace for Binary Nb-Ti BCC_A2

```
INFO: === CVMEngine.compute() START ===
INFO: Input parameters:
INFO:   systemId: BCC_A2_bin
INFO:   systemName: BCC (A2)
INFO:   temperature: 1200.0 K
INFO:   composition: [0.500000, 0.500000]
INFO:   numComponents: 2
FINE: Validating C-Matrix...
FINE: ✓ C-Matrix found
FINE: Validating composition...
FINE: ✓ Composition valid (sum=1.000000000)
FINE: Validating temperature...
FINE: ✓ Temperature valid
FINE: Loading C-Matrix from AllClusterData...
FINE: ✓ C-Matrix loaded
FINE: Creating CVMInput with topology...
FINE: ✓ CVMInput created
FINE: Evaluating ECI at T=1200.0 K...
FINE:   Basis: BCC_A2, numComponents=2, ncf=4, totalCfs=6
FINE:   CF names: [v4AB, v3AB, v22AB, v21AB, xA, xB]
FINE:   Processing 4 CEC terms...
FINER:     eci[2] (v22AB) = -260.0 + 0.0*1200.0 = -2.6000e+02
FINER:     eci[3] (v21AB) = -390.0 + 0.0*1200.0 = -3.9000e+02
FINE:   ECI array: [[2]=-2.60e+02, [3]=-3.90e+02]
FINE: ✓ ECI evaluated (4 non-point terms)
INFO: Running CVM N-R minimization...
INFO: CVMPhaseModel.ensureMinimized() — STARTING minimization
INFO:   systemId: BCC_A2_bin
INFO:   temperature: 1200.0 K
INFO:   composition: [0.5, 0.5]
INFO:   ncf: 4, tcf: 6
FINE:   tolerance: 1.0e-6
FINE:   ECI: [[2]=-2.60e+02, [3]=-3.90e+02]
FINE:   Calling NewtonRaphsonSolverSimple.solve()...
FINE: NewtonRaphsonSolverSimple.minimize — ENTER: ncf=4, tcf=6, T=1200.0, xB=0.5, tolerance=1.0e-6
FINE:   Initial CF (random state): [0.0625, 0.0, 0.25, 0.25]
FINEST: NewtonRaphsonSolverSimple — iter   1: G=-3.05e+01 H=-6.50e+01 S=1.08e-01 ||Gu||=2.65e+01
FINEST: NewtonRaphsonSolverSimple — iter  20: G=-4.23e+01 H=-8.12e+01 S=1.08e-01 ||Gu||=1.82e-04
FINEST: NewtonRaphsonSolverSimple — iter  23: G=-4.25e+01 H=-8.15e+01 S=1.08e-01 ||Gu||=6.31e-07
FINE: NewtonRaphsonSolverSimple — CONVERGED (gradient norm < tolerance): iter=23, ||Gu||=6.31e-07
FINE:   Newton-Raphson solver returned
FINE:   Cached 23 iteration snapshots
FINE:   Solver converged!
FINE:   Equilibrium CFs (non-point): [0.126, -0.051, 0.302, 0.198]
FINE:   Evaluating CVMFreeEnergy at equilibrium CFs...
FINE:   Thermodynamic values: G=-4.25e+01, H=-8.15e+01, S=1.08e-01
INFO: CVMPhaseModel.minimize — SUCCESS: T=1200.0 K, x=[0.5, 0.5], iterations=23, ||dG||=6.31e-07, G=-4.25e+01 J/mol, elapsed=15 ms
INFO: CVMPhaseModel.ensureMinimized() — COMPLETE
INFO: ✓ CVM minimization converged in 23 iterations (||Gu||=6.31e-07)
INFO:   G(eq) = -4.25e+01 J/mol
INFO:   H(eq) = -8.15e+01 J/mol
INFO:   S(eq) = 1.08e-01 J/(mol·K)
```

---

## Debugging Tips

### 1. **Solver Not Converging?**
- Check `||Gu||` in iteration trace — should decrease monotonically
- If it plateaus, check for negative CVs (stpmx limiting step)
- Increase `maxIter` or check ECI values for physics inconsistency

### 2. **Singular Hessian?**
- Usually means the free-energy landscape is degenerate
- Check that composition is not at an extreme (0 or 1)
- Verify enthalpy term is correct: H = Σ eci[l] * v[l]

### 3. **Slow Convergence?**
- Check N-R iteration count in CVMPhaseModel logs
- For binary systems, 20-50 iterations is normal
- Very slow convergence near phase boundaries

### 4. **Comparing Old vs New Basis**
- Run same temperature/composition with old orthogonal code
- Both should give same G, H, S values (to floating-point precision)
- N-R iteration count may differ slightly due to different random init

---

## Log Level Reference

| Level | Usage |
|-------|-------|
| `SEVERE` | Validation failures, exceptions |
| `WARNING` | Solver non-convergence, singular Hessian |
| `INFO` | Major milestones (START, SUCCESS, iteration count) |
| `FINE` | Parameter summaries, step completion |
| `FINER` | ECI term details, intermediate CF values |
| `FINEST` | Every N-R iteration snapshot |

---

## Redirecting Logs to File

### Using logging.properties

```properties
handlers = java.util.logging.FileHandler
java.util.logging.FileHandler.pattern = cvm_trace_%u_%g.log
java.util.logging.FileHandler.limit = 10485760
java.util.logging.FileHandler.count = 3
java.util.logging.FileHandler.formatter = java.util.logging.SimpleFormatter
java.util.logging.SimpleFormatter.format = [%1$tb %1$td, %1$tY %1$tl:%1$tM:%1$tS %1$Tp] %4$s %3$s: %5$s%6$s%n
```

### Programmatically

```java
try {
    FileHandler fh = new FileHandler("cvm_trace.log");
    fh.setLevel(Level.FINEST);
    SimpleFormatter formatter = new SimpleFormatter();
    fh.setFormatter(formatter);

    Logger.getLogger("org.ce.domain.engine.cvm").addHandler(fh);
    Logger.getLogger("org.ce.domain.engine.cvm").setLevel(Level.FINE);
} catch (IOException e) {
    e.printStackTrace();
}
```

---

## Performance Notes

- **INFO level:** Minimal overhead (~1% slower)
- **FINE level:** Moderate overhead (~5% slower)
- **FINEST level:** Significant overhead (~20% slower) — use only for debugging

For production runs, keep logging at **INFO** level or higher.
