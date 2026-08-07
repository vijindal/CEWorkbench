# Composition API redesign — Map-based, pycalphad-style (NOT YET IMPLEMENTED)

> Saved plan, to be implemented in a future session. See conversation history
> for full context. Decisions below were confirmed with the user before this
> plan was saved.

## Context

CE Thermodynamics Workbench's CLI and GUI both specify composition, but through
three inconsistent, positional-`double[]`-based conventions that meet only inside
`CalculationService.deriveComposition()`, which hardcodes "element index 0 is
the dependent species." This already caused one real, silently-wrong-answer bug
(CLI dropped the first mole-fraction argument for 3+ component systems — now
patched as a minimal fix, committed as `60b51a1`).

The user wants a proper fix, not just the patch: adopt pycalphad's convention —
composition specified as named `element -> fraction` pairs (`Map<String,Double>`),
not positional arrays — because this class of bug becomes structurally
impossible once "which index is which element" no longer exists. This code is
meant to be called as a backend API by other tools, not just CLI/GUI, so the
API surface matters beyond these two frontends.

Decisions already made with the user:
1. **Missing composition → throw an error.** No default guess (today's silent
   zero-fill is exactly what caused the original bug).
2. **Keep the one-varying-composition-axis limit, but make it an explicit,
   clear error** instead of today's silent `break`/truncation. Multi-axis
   scanning is out of scope (would need `CalculationResult.Grid` to become
   N-D, which is a separate future task).
3. **Remove `Mode.FINITE_SIZE_SCALING`** — `runFiniteSizeScan` is unreachable
   from the GUI, hardcodes `L={12,16,24}` ignoring `MCS_L`, and discards all
   but the last result. Deleting it also simplifies the Parameter design
   (`COMPOSITION` becomes purely the single-point input for `ANALYSIS`).
4. **Clean break on CLI syntax** — no positional-double fallback. Bare numbers
   where a composition token is expected produce a clear migration-hint error.

## Design

### 1. Canonical element order — `SystemId.elementList()`

Add to `Workspace.SystemId` (`src/main/java/org/ce/model/storage/Workspace.java`),
since it's already the identity object carried into `ModelSession` and already
owns the one existing `elements.split("-")` consumer (`ncompSuffix()`):

```java
public record SystemId(String elements, String structure, String model) {
    /** Canonical element order, index 0 = dependent species by convention. Never empty. */
    public List<String> elementList() {
        String[] parts = elements.trim().split("-");
        List<String> out = new ArrayList<>(parts.length);
        for (String p : parts) {
            String t = p.trim();
            if (!t.isEmpty()) out.add(t);
        }
        if (out.isEmpty())
            throw new IllegalArgumentException("No elements parsed from: '" + elements + "'");
        return List.copyOf(out);
    }
    public int numComponents() { return elementList().size(); }
    public int indexOf(String element) {
        List<String> list = elementList();
        for (int i = 0; i < list.size(); i++)
            if (list.get(i).equalsIgnoreCase(element)) return i;
        return -1;
    }
}
```

Element-symbol matching is **case-insensitive** (CLI users will type `ti=0.3`).
`ModelSession.numComponents()` delegates to `systemId.numComponents()`; add
`ModelSession.elements()` → `systemId.elementList()` as the accessor the
calculation layer reaches for.

Only migrate the `elements.split("-")` call sites directly in the
composition-handling path (CalculationService, ThermodynamicWorkflow
validation, CLI parsing, GUI compSpinners keying) — not all 12+ scattered
occurrences across the codebase.

### 2. New composition value types — `org.ce.calculation.Composition`

New file, `src/main/java/org/ce/calculation/Composition.java`:

```java
package org.ce.calculation;

import java.util.*;

/** A single-point composition: element symbol -> mole fraction. May be partial (K-1 entries). */
public record Composition(Map<String, Double> fractions) {
    public Composition { fractions = Map.copyOf(fractions); }

    public static Composition of(Map<String, Double> m) { return new Composition(m); }

    /** Resolves against canonical order to a full-K array. Throws on any invalid input. */
    public double[] resolve(List<String> canonicalOrder) {
        int K = canonicalOrder.size();
        double[] full = new double[K];
        boolean[] given = new boolean[K];

        for (var e : fractions.entrySet()) {
            int idx = indexOfIgnoreCase(canonicalOrder, e.getKey());
            if (idx < 0) throw new IllegalArgumentException(
                "Unknown element '" + e.getKey() + "' for system " + canonicalOrder
                + ". Valid: " + canonicalOrder);
            if (given[idx]) throw new IllegalArgumentException(
                "Duplicate entry for element '" + e.getKey() + "'");
            double v = e.getValue();
            if (Double.isNaN(v) || v < 0.0 || v > 1.0) throw new IllegalArgumentException(
                "Mole fraction for '" + e.getKey() + "' must be in [0,1], got " + v);
            full[idx] = v; given[idx] = true;
        }

        int missing = 0, missingIdx = -1;
        for (int i = 0; i < K; i++) if (!given[i]) { missing++; missingIdx = i; }
        double sum = 0; for (int i = 0; i < K; i++) if (given[i]) sum += full[i];

        if (K == 0) throw new IllegalStateException("Empty canonical element order");
        if (fractions.isEmpty()) throw new IllegalArgumentException(
            "No composition specified. Provide at least " + (K - 1) + " of " + K
            + " element=fraction pairs for system " + canonicalOrder + ".");

        if (missing == 0) {
            if (Math.abs(sum - 1.0) > 1e-6) throw new IllegalArgumentException(
                "Fully-specified composition must sum to 1.0, got " + sum + ": " + fractions);
        } else if (missing == 1) {
            if (sum > 1.0 + 1e-6) throw new IllegalArgumentException(
                "Specified fractions sum to " + sum + " > 1; no room for derived element '"
                + canonicalOrder.get(missingIdx) + "'");
            full[missingIdx] = 1.0 - sum;
        } else {
            List<String> missingNames = new ArrayList<>();
            for (int i = 0; i < K; i++) if (!given[i]) missingNames.add(canonicalOrder.get(i));
            throw new IllegalArgumentException(
                "Underdetermined: " + missing + " elements unspecified " + missingNames
                + ", need at most 1 unspecified. Specify " + (K - 1) + " of " + K + " fractions.");
        }
        return full;
    }

    private static int indexOfIgnoreCase(List<String> list, String sym) {
        for (int i = 0; i < list.size(); i++)
            if (list.get(i).equalsIgnoreCase(sym)) return i;
        return -1;
    }
}

/** A scan range over one element's mole fraction. */
record ScanRange(double start, double end, double step) {
    static ScanRange fixed(double v) { return new ScanRange(v, v, 0.0); }
    boolean varies() { return Math.abs(start - end) > 1e-9; }
    int pointCount() { return varies() ? (int) Math.round((end - start) / step) + 1 : 1; }
    double valueAt(int i) { return varies() ? start + i * step : start; }
}

/** A composition scan: element symbol -> range. At most one range may vary. */
record CompositionScan(Map<String, ScanRange> ranges) {
    CompositionScan {
        ranges = new LinkedHashMap<>(ranges); // preserve insertion order for reproducible grids
        long varyingCount = ranges.values().stream().filter(ScanRange::varies).count();
        if (varyingCount > 1) {
            List<String> varying = ranges.entrySet().stream()
                .filter(e -> e.getValue().varies()).map(Map.Entry::getKey).toList();
            throw new IllegalArgumentException(
                "Only one composition axis may vary per scan; found " + varying
                + " both varying. Fix all but one to a single value.");
        }
        ranges = Collections.unmodifiableMap(ranges);
    }

    static CompositionScan fixedAt(Composition c) {
        Map<String, ScanRange> r = new LinkedHashMap<>();
        for (var e : c.fractions().entrySet()) r.put(e.getKey(), ScanRange.fixed(e.getValue()));
        return new CompositionScan(r);
    }

    int pointCount() {
        return ranges.values().stream().filter(ScanRange::varies)
            .findFirst().map(ScanRange::pointCount).orElse(1);
    }

    Composition compositionAt(int i) {
        Map<String, Double> m = new LinkedHashMap<>();
        for (var e : ranges.entrySet()) m.put(e.getKey(), e.getValue().valueAt(i));
        return new Composition(m);
    }

    void validateAgainst(List<String> canonicalOrder) {
        // cheap fail-fast: resolve the first and last grid point before running anything
        compositionAt(0).resolve(canonicalOrder);
        compositionAt(pointCount() - 1).resolve(canonicalOrder);
    }
}
```

`Composition.resolve()` is the single chokepoint replacing
`CalculationService.deriveComposition` — no hardcoded index 0 (the derived
element is whichever one is *omitted*), no possible length mismatch,
non-negativity checked, sum-overflow checked, unknown names rejected, empty
input throws per decision #1.

### 3. Parameter redesign

In `CalculationDescriptor.java`:
- Retype `COMPOSITION` to `Composition.class` (single point).
- Add `COMPOSITION_SCAN` of type `CompositionScan.class`.
- **Delete** `X_STARTS`, `X_ENDS`, `X_STEPS` (the K-1 positional triple).
- **Keep** `X_START`/`X_END`/`X_STEP` (singular) — GUI spinner-editor templates
  only, never read as job parameters; add a comment clarifying they're
  unrelated to the deleted plural trio.
- **Delete** `Mode.FINITE_SIZE_SCALING` and its `Registry.getRequirements` branch
  (decision #3).
- `Registry.getRequirements(ANALYSIS)` → `{T_START, T_END, T_STEP, COMPOSITION_SCAN}`.

In `JobSpecifications`, close the unchecked-cast hole (this is the mechanism
that let the original bug stay invisible):

```java
public void set(Parameter param, Object value) {
    if (value != null && !param.type.isInstance(value))
        throw new IllegalArgumentException(
            "Parameter '" + param.name + "' expects " + param.type.getSimpleName()
            + ", got " + value.getClass().getSimpleName());
    parameters.put(param, value);
}

public <T> T require(Parameter param, Class<T> type) {
    Object v = parameters.getOrDefault(param, param.defaultValue);
    if (v == null) throw new IllegalStateException("Required parameter '" + param.name + "' not set");
    return type.cast(v);
}
```

Keep `getOrDefault`/`get` for optional params; use `require` for mandatory ones.

### 4. `CalculationService.runAnalysis` rewrite

Delete `deriveComposition` and the `Varying` record. Replace the
float-accumulating `for (double val=start; val<=end+1e-9; val+=step)` loops
(a real latent bug — repeated `+=` can produce 10 or 11 points for a 0..1
step-0.1 scan depending on rounding) with index-based loops using
`ScanRange.pointCount()`/`valueAt()`.

```java
private CalculationResult.Grid runAnalysis(ModelSession session, JobSpecifications jobSpecs,
        Consumer<String> textSink, Consumer<ProgressEvent> eventSink) throws Exception {

    List<String> order = session.systemId.elementList();

    CompositionScan scan = jobSpecs.get(Parameter.COMPOSITION_SCAN)
        .map(CompositionScan.class::cast)
        .orElseGet(() -> {
            Composition c = jobSpecs.get(Parameter.COMPOSITION)
                .map(Composition.class::cast)
                .orElseThrow(() -> new IllegalStateException(
                    "No composition specified (neither COMPOSITION nor COMPOSITION_SCAN set)."));
            return CompositionScan.fixedAt(c);
        });

    scan.validateAgainst(order); // fail fast before running anything expensive

    double tStart = jobSpecs.require(Parameter.T_START, Double.class);
    double tEnd   = jobSpecs.require(Parameter.T_END, Double.class);
    double tStep  = jobSpecs.require(Parameter.T_STEP, Double.class);
    boolean tVaries = Math.abs(tStart - tEnd) > 1e-6;
    int tPoints = tVaries ? (int) Math.round((tEnd - tStart) / tStep) + 1 : 1;

    List<List<ThermodynamicResult>> grid = new ArrayList<>();
    for (int ti = 0; ti < tPoints; ti++) {
        double T = tVaries ? tStart + ti * tStep : tStart;
        List<ThermodynamicResult> row = new ArrayList<>();
        for (int xi = 0; xi < scan.pointCount(); xi++) {
            double[] x = scan.compositionAt(xi).resolve(order);
            row.add(thermoWorkflow.runCalculation(session, new ThermodynamicWorkflow.Request(
                T, x, jobSpecs.getProperty(), textSink, eventSink,
                jobSpecs.getOrDefault(Parameter.MCS_L), jobSpecs.getOrDefault(Parameter.MCS_NEQUIL),
                jobSpecs.getOrDefault(Parameter.MCS_NAVG), null)));
        }
        grid.add(row);
    }
    return new CalculationResult.Grid(grid);
}
```

Also delete `runFiniteSizeScan` entirely (decision #3) and its call in `execute()`'s
mode switch.

### 5. `ThermodynamicWorkflow.validateInputs` — narrow, don't extend

Since `Composition.resolve()` now owns full validation for the Map-based path,
`validateInputs` becomes a cheap defensive assertion for any hand-built
`Request` (e.g. `RunComparison.java`). Pass `session` in to add the one check
`resolve()` can't guarantee was run:

```java
private void validateInputs(ModelSession session, double T, double[] x) {
    if (T < 0) throw new IllegalArgumentException("Temperature cannot be negative: " + T);
    if (x == null || x.length == 0) throw new IllegalArgumentException("Composition array missing");
    if (x.length != session.numComponents()) throw new IllegalArgumentException(
        "Composition length " + x.length + " != numComponents " + session.numComponents());
    double sum = 0;
    for (double val : x) {
        if (val < 0 || val > 1) throw new IllegalArgumentException("Mole fraction out of [0,1]: " + val);
        sum += val;
    }
    if (Math.abs(sum - 1.0) > 1e-4)
        throw new IllegalArgumentException("Composition does not sum to 1.0: " + Arrays.toString(x));
}
```

### 6. CLI syntax (`Main.java`)

Named `Sym=frac` tokens, positionally free among the trailing args:

```
calc_min <elements> <structure> <model> <temp> <El>=<x> [<El>=<x> ...] [G|H|S] [--verbose]

  ./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"
  ./gradlew run --args="calc_min Nb-Ti-V-Zr BCC_A2 T 1000 Ti=0.3 V=0.2 Zr=0.25"
```

Scan form (`type2`), `=start:end:step`:

```
./gradlew run --args="type2 Nb-Ti BCC_A2 T --temp 800:1200:100 --x Ti=0.1:0.9:0.1"
```

Shared parser (used by both `calc_min` and `type2` so they can't diverge):

```java
/** Parses "Ti=0.3" or "Ti=0.1:0.9:0.1". Returns null if not a composition token. */
private static Map.Entry<String, ScanRange> parseCompToken(String tok) {
    int eq = tok.indexOf('=');
    if (eq <= 0) return null;
    String sym = tok.substring(0, eq).trim();
    String val = tok.substring(eq + 1).trim();
    String[] parts = val.split(":");
    try {
        return switch (parts.length) {
            case 1 -> Map.entry(sym, ScanRange.fixed(Double.parseDouble(parts[0])));
            case 3 -> Map.entry(sym, new ScanRange(Double.parseDouble(parts[0]),
                                                   Double.parseDouble(parts[1]),
                                                   Double.parseDouble(parts[2])));
            default -> throw new IllegalArgumentException(
                "Bad composition token '" + tok + "'. Use El=x or El=start:end:step");
        };
    } catch (NumberFormatException e) {
        throw new IllegalArgumentException("Non-numeric mole fraction in '" + tok + "'");
    }
}
```

Clean break (decision #4): if a trailing arg after `<temp>` doesn't match
`El=...` (i.e. `parseCompToken` returns null) and isn't `G`/`H`/`S`/`--verbose`,
emit a migration-hint error:

```
Error: positional mole fractions are no longer supported.
  Old: calc_min Nb-Ti-V BCC_A2 T 1000 0.3 0.4
  New: calc_min Nb-Ti-V BCC_A2 T 1000 Ti=0.3 V=0.4
```

Update `type2`'s equiatomic composition (currently the dead-COMPOSITION-write
bug) to build a real `Composition` via equal 1/K fractions for K-1 elements
and set `Parameter.COMPOSITION`.

Update `CLAUDE.md`'s "Build and run" example (`calc_min Nb-Ti BCC_A2 T 1000 0.5`
→ `... T 1000 Ti=0.5`).

### 7. GUI (`DynamicCalculationPanel.java`)

Fix the double-split root cause (`session.systemId.elements()` at form-build
time vs. `editorText(elementsCombo)` at request-build time can disagree) with
one shared accessor:

```java
private SystemId currentSystemId() {
    String el = editorText(elementsCombo);
    if (el.isBlank()) {
        ModelSession s = context.getActiveSession();
        if (s != null) el = s.systemId.elements();
    }
    return new SystemId(el, editorText(structureCombo), editorText(modelCombo));
}
```

Use it in both `rebuildParameterForm()` and `startExecution()`. In
`startExecution()`, build the `CompositionScan` map directly and **fail loudly**
(not silently leave 0.0) if a spinner is missing for the current element list:

```java
Map<String, ScanRange> ranges = new LinkedHashMap<>();
List<String> elems = currentSystemId().elementList();
for (String sym : elems.subList(1, elems.size())) {
    JSpinner[] sp = compSpinners.get(sym);
    if (sp == null) {
        logSink.accept("Error: composition input missing for '" + sym
            + "'. Re-select the property to rebuild the form.");
        return;
    }
    ranges.put(sym, new ScanRange((Double) sp[0].getValue(),
                                  (Double) sp[1].getValue(),
                                  (Double) sp[2].getValue()));
}
specs.set(Parameter.COMPOSITION_SCAN, new CompositionScan(ranges));
```

Add a listener on `elementsCombo` that calls `rebuildParameterForm()` on
change, so the form can't go stale relative to the combo in the first place
(root-cause fix; the guard above is the safety net).

`ParameterFieldFactory.createEditor` needs an explicit branch for
`Composition.class`/`CompositionScan.class` that returns null/throws, so a
future Parameter of these types can't silently fall through to a raw
`JTextField`.

## Files to change

| File | Change |
|---|---|
| `org/ce/calculation/Composition.java` | **New.** `Composition`, `ScanRange`, `CompositionScan`. |
| `org/ce/model/storage/Workspace.java` | Add `SystemId.elementList()`, `numComponents()`, `indexOf()`; route `ncompSuffix()` through them. |
| `org/ce/model/ModelSession.java` | `numComponents()` delegates to `systemId`; add `elements()` accessor. |
| `org/ce/calculation/CalculationDescriptor.java` | Retype `COMPOSITION`, add `COMPOSITION_SCAN`, delete `X_STARTS/X_ENDS/X_STEPS` and `Mode.FINITE_SIZE_SCALING`; update `Registry.getRequirements`; type-check `JobSpecifications.set`, add `require`. |
| `org/ce/calculation/workflow/CalculationService.java` | Rewrite `runAnalysis` (index-based grid, `CompositionScan`); delete `deriveComposition`, `Varying`, `runFiniteSizeScan`. |
| `org/ce/calculation/workflow/thermo/ThermodynamicWorkflow.java` | Pass `session` into `validateInputs`; add length/non-negativity checks. |
| `org/ce/ui/cli/Main.java` | `parseCompToken` helper shared by `calc_min`/`type2`; clean-break migration error for bare positional numbers; fix `type2`'s dead-COMPOSITION-write bug; update usage strings. |
| `org/ce/ui/gui/DynamicCalculationPanel.java` | `currentSystemId()` single source; build `CompositionScan` directly; fail-loud on missing spinner; rebuild form on elements-combo change. |
| `org/ce/ui/gui/ParameterFieldFactory.java` | Explicit branch for `Composition.class`/`CompositionScan.class`. |
| `org/ce/scratch/RunComparison.java` | Update or leave `double[]` Requests as-is (still valid at the Request boundary) — decide during implementation whether to migrate its examples to the new API for dogfooding. |
| `CLAUDE.md` | Update `calc_min` usage example. |

`MCSRunner`, `LatticeConfig`, `MetropolisMC`, `AlloyMC`, `CVMGibbsModel` —
**untouched**. The `double[]` boundary at `ThermodynamicWorkflow.Request` is
preserved exactly, so MCS quantization/re-derivation behavior (discretizing to
integer site counts, caching the *realized* composition) is unaffected.

## Verification plan

1. `./gradlew compileJava` — must build clean.
2. Regression-check the exact CLI invocations already validated this session:
   - `calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5` → G should still be **-3480.520906 J/mol** (binary, iteration-by-iteration validated earlier).
   - `calc_min Nb-Ti-V BCC_A2 T 1000 Ti=0.33 V=0.34` → G should still be **-7051.125730 J/mol** (ternary, x=(0.33,0.33,0.34)).
   - `calc_min Nb-Ti-V-Zr BCC_A2 T 1273 Ti=0.25 V=0.25 Zr=0.25 S` → S should still be **11.0812146249 J/(mol·K)** (quaternary, matches the published paper's 11.08).
3. Confirm the old positional syntax now errors with the migration-hint message rather than silently computing a wrong composition.
4. Confirm a composition Map with 2+ missing entries throws "Underdetermined"; a Map summing >1 with one omitted throws "no room for derived element"; an unknown element name throws "Unknown element".
5. Confirm `CompositionScan` with two varying ranges throws the one-axis-limit error.
6. GUI: launch `runGui`, change the elements combo after building the parameter form, confirm the form rebuilds (no stale spinners) and a `calc_min`-equivalent run through the panel produces the same G as the CLI for the same system/composition/temperature.
7. Confirm `Mode.FINITE_SIZE_SCALING` and `runFiniteSizeScan` are fully removed and nothing references them (`grep -r FINITE_SIZE_SCALING src/`).

## Other CLI/GUI architectural findings not yet fixed (from the same audit)

These were found during the "deep audit" that motivated this plan, but are
separate from the composition redesign above — worth fixing in a follow-up pass:

- **V4** (MED-HIGH): `CECManagementPanel.java:436-461` builds a binary→ternary
  species mapping (real model semantics) directly inside a Swing panel instead
  of delegating to `CECManagementWorkflow`.
- **V5** (MED): `CECManagementPanel.java:134-154` hand-rolls Hamiltonian/cluster
  ID derivation instead of using `SystemId`; diverges from the canonical
  version for CVCF models (no `_CVCF` special-case) and silently mislabels
  instead of throwing for ncomp≥5.
- **V2/V3** (MED): `Main.java:138` and `DataPreparationPanel.java:389` both
  call `ClusterCFIdentificationPipeline.runFullWorkflow` directly, bypassing
  the calculation layer, even though `CECManagementWorkflow` already wraps
  this call correctly.
- **V6** (LOW-MED): `MainWindow.java:91-99` reverse-engineers "was this a
  T-scan or x-scan" from result data using tolerance 1e-3, while the
  calculation layer used tolerance 1e-6 to decide it originally — potential
  mislabeling for fine-step scans. `CalculationResult.Grid` should carry the
  `Varying`/scan-axis descriptor instead of discarding it.
- **V7** (cosmetic): `OutputPanel.java:3` has an unused import from the model
  layer in a pure-rendering class.
