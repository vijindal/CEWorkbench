# API Reference

Calling the CE Thermodynamics Workbench from your own code.

Two interfaces, same engine underneath:

- **[Java API](#java-api)** — in-process calls from any JVM language.
- **[JSON API](#json-api)** — one JSON request in, one JSON response out, for
  Python, MATLAB, shell, or anything else.

The CLI and GUI both route through the same entry point, so results are
identical regardless of how you call it.

For running calculations interactively, see [USER_GUIDE.md](USER_GUIDE.md).

---

## Conventions

These apply to both interfaces.

**Composition is named, never positional.** You give `element → mole fraction`
pairs. Any one element may be omitted and is derived so the total is 1. Omitting
two or more is an error.

**One varying axis per scan.** Either temperature or one element's fraction may
sweep — not both. Multi-dimensional grids require repeated calls.

**ECIs are matched by name and validated strictly.** A name that doesn't resolve
against the basis is rejected outright rather than silently contributing zero.
Incomplete coverage is likewise rejected.

**ECIs must already be in the CVCF basis, in J/mol.** No basis transformation or
unit conversion is performed. Each term is evaluated as `a + b·T`.

**Check `converged`.** The CVM minimizer can hit its iteration limit and still
return plausible-looking numbers.

---

## Java API

`CalculationService.calculate` / `calculateScan` are the public entry points.

```java
CEWorkbenchContext ctx = new CEWorkbenchContext();
CalculationService service = ctx.getCalculationService();

// Build the session once per (elements, structure, model, engine) identity.
ModelSpecifications specs =
        new ModelSpecifications("Nb-Ti", "BCC_A2", "T", EngineConfig.CVM);
ModelSession session = service.getOrBuildSession(specs, null);

// Single point — Nb derived as 1 - 0.5.
ThermodynamicResult r = service.calculate(
        session,
        new Conditions(1000.0, Map.of("Ti", 0.5)),
        Property.GIBBS_ENERGY,
        null, null);

System.out.println(r.gibbsEnergy);   // -3480.5209063901
if (Boolean.FALSE.equals(r.converged)) { /* do not trust r */ }

// Temperature scan.
List<ThermodynamicResult> scan = service.calculateScan(
        session,
        new ConditionsScan(new Range(800, 1200, 100), Map.of("Ti", Range.fixed(0.5))),
        Property.GIBBS_ENERGY,
        null, null);
```

The two trailing arguments are optional progress sinks:
`Consumer<String>` for log lines and `Consumer<ProgressEvent>` for structured
events. Pass `null` to ignore them.

For Monte Carlo, supply algorithm parameters explicitly:

```java
service.calculate(session, conditions, Property.ENTHALPY,
                  new CalculationService.McsParams(8, 1000, 2000), null, null);
```

`McsParams` holds supercell size `L` and sweep counts — algorithm settings, not
physical conditions, which is why they are separate from `Conditions`.

### Result fields

| Field | Notes |
|---|---|
| `gibbsEnergy`, `enthalpy`, `entropy` | J/mol, J/mol, J/(mol·K). `NaN` when the engine doesn't produce it |
| `optimizedCFs` | Equilibrium correlation functions (CVM) |
| `avgCFs`, `stdCFs` | Sampled correlation functions (MCS) |
| `sro` | Cowley-Warren parameters by shell — see [below](#short-range-order) |
| `converged`, `iterations`, `finalGradientNorm` | Minimizer status; `null` for MCS |

### Using it as a dependency

During development the simplest route for another Gradle project is a composite
build — no publishing step, and changes are picked up immediately:

```groovy
// settings.gradle
includeBuild('../CEWorkbench')
```

---

## JSON API

For non-JVM callers. One JSON request on stdin, one JSON response on stdout.
Diagnostics go to **stderr**, so stdout is always machine-parseable.

```bash
./gradlew installDist                                    # once — builds the launcher
build/install/CEWorkbench/bin/CEWorkbench api < request.json
```

Exit code is 0 on success, 1 on error. A JSON error object is written either
way, so the reason is always readable.

> Set `CEWORKBENCH_DATA` to your `data/CEWorkbench` directory when invoking from
> another working directory, or the workbench won't find its inputs.

### Discovering what a system expects

```json
{"describe": {"elements":"Nb-Ti","structure":"BCC_A2","model":"T","engine":"CVM"}}
```
```json
{"ok":true, "supported":true, "ncf":4,
 "expectedEciNames":["v4AB","v3AB","v22AB","v21AB"],
 "calculations":["GIBBS_ENERGY","ENTHALPY","ENTROPY"],
 "notImplemented":["PHASE_EQUILIBRIUM"]}
```

Query this instead of hard-coding ECI names — it is the reliable guard against
supplying names that don't match the basis.

### Running a calculation

```json
{
  "system":      {"elements":"Nb-Ti","structure":"BCC_A2","model":"T","engine":"CVM"},
  "calculation": "GIBBS_ENERGY",
  "conditions":  {"temperature":1000, "composition":{"Ti":0.5}}
}
```

`conditions.temperature` and each composition entry accept either a number or a
`{"start":…,"end":…,"step":…}` range; at most one may vary per request.

`system.model` must be the bare model name (`"T"`), not `"T_CVCF"` — the suffix
is applied internally to the stored Hamiltonian ID.

### Supplying your own ECIs

Omit `"hamiltonian"` to use the stored ECI database, or include it to pass
values directly:

```json
"hamiltonian": {
  "basis": "CVCF", "units": "J/mol",
  "cecTerms": [
    {"name":"e4AB","a":0.0,"b":0.0},    {"name":"e3AB","a":0.0,"b":0.0},
    {"name":"e22AB","a":3120.0,"b":0.0},{"name":"e21AB","a":6240.0,"b":0.0}
  ]
}
```

`basis` and `units` are declared explicitly and anything other than
`CVCF` / `J/mol` is rejected, so a mismatch fails loudly instead of silently
producing wrong energies.

Accepted spellings per term: canonical `v…`, the `e…` alias, the published
`e2<species><shell>` pair notation, and `CF_<index>`.

### Response

```json
{"ok": true,
 "system": {"elements":["Nb","Ti"],"structure":"BCC_A2","model":"T","engine":"CVM"},
 "hamiltonianSource": "inline",
 "calculation": "GIBBS_ENERGY",
 "points": [{
   "temperature": 1000.0,
   "converged": true,
   "composition": {"Nb":0.5, "Ti":0.5},
   "gibbsEnergy": -3480.5209063901,
   "enthalpy": 2214.3770000563,
   "entropy": 5.6948979064,
   "correlationFunctions": {"v4AB":0.0536, "v3AB":0.0, "v22AB":0.2392, "v21AB":0.2353},
   "sro": { "1NN": {...}, "2NN": {...} }
 }]}
```

Quantities an engine doesn't produce are **omitted**, never emitted as `NaN`
(which isn't valid JSON). When `converged` is `false` the point also carries
`iterations` and `finalGradientNorm`.

### Short-range order

CVM points include Cowley-Warren parameters, computed from the converged cluster
probabilities:

```
α_PR = 1 − p_PR / (x_P x_R)      α < 0 ordering · α = 0 random · α > 0 clustering
```

```json
"sro": {
  "1NN": {"Mo-Nb": {"alpha": -0.05088, "probability": 0.11617, "random": 0.11111},
          "Mo-Ta": {"alpha": -0.05541, "probability": 0.11791, "random": 0.11111}},
  "2NN": { ... }
}
```

All species pairs are reported, like-pairs included.

### Errors

```json
{"ok":false,"error":"ECI_VALIDATION_FAILED",
 "message":"...",
 "unmatched":["e21XX"], "unmapped":["v21AB"],
 "expected":["v4AB","v3AB","v22AB","v21AB"]}
```

| `error` | Cause |
|---|---|
| `INVALID_SYSTEM` | Missing or malformed elements/structure/model/engine |
| `UNSUPPORTED_SYSTEM` | No CVCF basis registered for that structure + component count |
| `MISSING_CONDITIONS` / `INVALID_CONDITIONS` | Absent, malformed, or multi-axis conditions |
| `UNSUPPORTED_BASIS` / `UNSUPPORTED_UNITS` | `hamiltonian.basis`/`units` not `CVCF`/`J/mol` |
| `INVALID_HAMILTONIAN` | Malformed `cecTerms` |
| `ECI_VALIDATION_FAILED` | Names don't map cleanly onto the basis |
| `NOT_IMPLEMENTED` | Calculation not available for this engine |
| `INTERNAL_ERROR` | Anything else; `message` carries the detail |

### Python example

```python
import json, os, subprocess

CEW_HOME = os.path.abspath("build/install/CEWorkbench")
EXE = os.path.join(CEW_HOME, "bin",
                   "CEWorkbench.bat" if os.name == "nt" else "CEWorkbench")

def call(request):
    p = subprocess.run([EXE, "api"], input=json.dumps(request),
                       capture_output=True, text=True)
    resp = json.loads(p.stdout)          # stderr holds logs; stdout is pure JSON
    if not resp.get("ok"):
        raise RuntimeError(f"{resp['error']}: {resp.get('message','')}")
    return resp

# Ask which ECIs this system needs
names = call({"describe": {"elements":"Nb-Ti","structure":"BCC_A2",
                           "model":"T","engine":"CVM"}})["expectedEciNames"]

# Compute with your own ECIs (CVCF basis, J/mol)
resp = call({
    "system":      {"elements":"Nb-Ti","structure":"BCC_A2","model":"T","engine":"CVM"},
    "hamiltonian": {"basis":"CVCF","units":"J/mol","cecTerms":[
                        {"name":"e4AB","a":0.0},    {"name":"e3AB","a":0.0},
                        {"name":"e22AB","a":3120.0},{"name":"e21AB","a":6240.0}]},
    "calculation": "GIBBS_ENERGY",
    "conditions":  {"temperature":{"start":800,"end":1200,"step":100},
                    "composition":{"Ti":0.5}},
})

for pt in resp["points"]:
    if not pt.get("converged", True):
        print(f"  warning: T={pt['temperature']} did not converge")
        continue
    print(pt["temperature"], pt["gibbsEnergy"])
```

## Ternary isothermal sections (`ternary_grid`)

A separate JSON stdin/stdout subcommand for composition-grid scans over a
3-component system at fixed temperature — the data behind isothermal-section
plots (Figs. 15–20 of Jindal & Lele 2025). `api`'s `ConditionsScan` only
supports one varying axis at a time (temperature XOR one composition
element); `ternary_grid` sweeps the full 2-D composition triangle directly,
in-process, reusing the session cache — no per-point subprocess overhead.

```json
{"system":      {"elements":"Nb-Ti-V","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "GIBBS_ENERGY",
 "temperature": 1273,
 "n": 20}
```

`n` is the grid resolution (subdivisions per triangle edge; point count is
`(n+1)(n+2)/2`, default 20). `calculation` accepts `GIBBS_ENERGY`,
`ENTHALPY`, `ENTROPY`, or `SRO`.

### Short-range order (`"calculation": "SRO"`)

Plots the 1st-neighbour Cowley-Warren pair SRO parameter (α, Eq. 40 of
Jindal & Lele 2025) for one unlike species pair across the composition
triangle — the data behind Fig. 24 of that paper. Requires a `"pair"` field:
a 2-element array naming the two elements (both must be in `system.elements`):

```json
{"system":      {"elements":"Nb-Ti-V","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "SRO",
 "pair":        ["Nb","Ti"],
 "temperature": 1273,
 "n": 20}
```

Only 1st-neighbour pair SRO is exposed for now — 2nd-neighbour and
triangle/tetrahedron multi-site SRO are not yet available via this endpoint.
(Some CVCF correlation functions, like the ternary binary-triangle CFs
`v3AB`/`v3AC`/`v3BC`, are antisymmetric *differences* of two cluster
probabilities rather than a single probability, so they don't have a directly
meaningful Cowley-Warren-style alpha — see CLAUDE.md for the full discussion.)

`α` is mathematically undefined at a pure-element corner (division by a
mole fraction that's exactly 0) — no point is emitted there at all for an
SRO request, rather than a `NaN` value.

### Response

```json
{"ok": true,
 "elements": ["Nb","Ti","V"], "structure":"BCC_A2", "model":"T", "engine":"CVM",
 "temperature": 1273.0, "calculation": "GIBBS_ENERGY", "skipped": 0,
 "points": [
   {"Nb":1.0, "Ti":0.0, "V":0.0, "value":0.0, "region": "CORNER"},
   {"Nb":0.0, "Ti":0.1, "V":0.9, "value":-2340.2, "region": "EDGE"},
   {"Nb":0.4, "Ti":0.3, "V":0.3, "value":-8500.1, "region": "INTERIOR"}
 ]}
```

Every point falls into exactly one of three regions, marked by `"region"`:

- **`INTERIOR`** (all three mole fractions > 0) — a genuine ternary CVM
  solve. The ternary solver can fail to converge in a thin composition band
  adjacent to a binary edge (a known near-edge instability, not a bug); such
  points are simply omitted rather than interpolated or estimated.
- **`EDGE`** (exactly one mole fraction = 0) — a genuine **binary** CVM
  solve on the sub-Hamiltonian extracted for that pair, not the ternary
  Hamiltonian evaluated at zero composition (which was tried and found
  numerically fragile exactly there). See CLAUDE.md for why this is
  physically sound (binary CECs are inherited unchanged into the ternary
  Hamiltonian).
- **`CORNER`** (two mole fractions = 0, i.e. a pure element) — no
  calculation; G/H/S are analytically 0, and no point is emitted for SRO
  (undefined there).

`"skipped"` counts grid points where a calculation was attempted but failed
(solver non-convergence, or — for an edge — no matching binary CVCF terms
found for that pair). For an SRO request, only the requested pair's own edge
is attempted — the other two edges don't involve both elements of the pair
at all, so those points are never attempted and never counted in
`"skipped"`, exactly like corners. This means `"skipped"` is directly
comparable between a G/H/S request and an SRO request on the same grid: it
always means "a real computation failed," never "this quantity doesn't
apply here."

### Getting a rendered image, not just numbers

Add `"render": true` to have the server also produce a plotted PNG (mpltern
ternary contour, matching the GUI's ternary panel) and return it embedded as
base64 — no shared filesystem between caller and server required:

```json
{"system": {"elements":"Nb-Ti-V","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "GIBBS_ENERGY", "temperature": 1273, "n": 20, "render": true}
```

```json
{"ok": true, "...": "...",
 "image": {"format": "png", "base64": "iVBORw0KGgoAAAANS..."}}
```

Points are always included regardless of `render`, so callers who only want
raw numbers pay no rendering cost. Rendering shells out to
`scripts/isothermal_section.py` (mpltern) — requires a working `python` on
`PATH` with `matplotlib`/`mpltern`/`numpy` installed on the machine running
the server. If rendering fails, the response still carries `"ok": true` and
the computed `points`, with the failure reported as `"renderError"` instead
of `"image"` — a broken Python environment never blocks the physics.

```python
resp = call_ternary_grid({  # same call() pattern as above, subcommand "ternary_grid"
    "system": {"elements":"Nb-Ti-V","structure":"BCC_A2","model":"T","engine":"CVM"},
    "calculation": "GIBBS_ENERGY", "temperature": 1273, "n": 30, "render": True,
})
if "image" in resp:
    import base64
    with open("isothermal_section.png", "wb") as f:
        f.write(base64.b64decode(resp["image"]["base64"]))
elif "renderError" in resp:
    print("render failed, using raw points instead:", resp["renderError"])
```

---

## Quaternary square plots (`quaternary_square`)

A separate JSON stdin/stdout subcommand for composition-grid scans over a
4-component system at fixed temperature — the data behind the "square plot"
of Fig. 20, Jindal & Lele 2025. The unit `(X,Y)` square is mapped onto the
quaternary composition simplex via that figure's parametrization; see
CLAUDE.md for the formula and why it covers only a 2-D slice of the full 3-D
simplex.

**There is no `slotOrder` or `variant` request field.** Given
`system.elements` as `A-B-C-D`, the command always computes and returns
**exactly two** square parametrizations — `A-B-C-D` and `A-B-D-C` (the last
two elements swapped) — because together they reach all six binary edges of
the composition tetrahedron; a single square only reaches four. This is
deliberate: the caller doesn't have to know or choose which slot ordering to
ask for.

```json
{"system":      {"elements":"Nb-Ti-V-Zr","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "GIBBS_ENERGY",
 "temperature": 1273,
 "n": 50}
```

`n` is the grid resolution (subdivisions per square axis; point count is
`(n+1)^2`, default **50** — note this differs from `ternary_grid`'s default
of 20). `calculation` accepts `GIBBS_ENERGY`, `ENTHALPY`, `ENTROPY`, or `SRO`.

### Short-range order (`"calculation": "SRO"`)

Same idea as the ternary case: the 1st-neighbour Cowley-Warren pair SRO
parameter for one unlike species pair, requiring a `"pair"` field (both
elements must be in `system.elements`):

```json
{"system":      {"elements":"Nb-Ti-V-Zr","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "SRO",
 "pair":        ["Nb","Ti"],
 "temperature": 1273,
 "n": 50}
```

Only 1st-neighbour pair SRO is exposed, for the same reasons as the ternary
case — see that section above. `α` is undefined at a pure-element corner; no
point is emitted there for an SRO request.

### Response

The top-level response always carries a `"results"` array with **exactly two
entries**, one per slot ordering, each with its own `points`/`skipped`
(and, if requested, its own `image`):

```json
{"ok": true,
 "structure": "BCC_A2", "model": "T", "engine": "CVM",
 "temperature": 1273.0, "calculation": "GIBBS_ENERGY",
 "results": [
   {"slotOrder": "Nb-Ti-V-Zr",
    "elements": ["Nb","Ti","V","Zr"],
    "skipped": 0,
    "points": [
      {"x":0.0, "y":0.0, "Nb":0.0, "Ti":1.0, "V":0.0, "Zr":0.0, "value":0.0, "region":"CORNER"},
      {"x":0.5, "y":0.5, "Nb":0.25, "Ti":0.25, "V":0.25, "Zr":0.25, "value":-8123.4, "region":"INTERIOR"}
    ]},
   {"slotOrder": "Nb-Ti-Zr-V",
    "elements": ["Nb","Ti","Zr","V"],
    "skipped": 0,
    "points": [ "..." ]}
 ]}
```

Each result's points use that result's own `elements` order for the
per-element mole-fraction keys — `slotOrder: "Nb-Ti-Zr-V"` names its points'
fields `Nb`/`Ti`/`Zr`/`V`, not `Nb`/`Ti`/`V`/`Zr`.

Every point falls into exactly one of **three** regions — one fewer than the
ternary case, because the square's geometry is different from the triangle's:

- **`INTERIOR`** (all four mole fractions > 0) — a genuine quaternary CVM
  solve.
- **`SQUARE_EDGE_BINARY`** (exactly two mole fractions = 0) — a genuine
  **binary** CVM solve, extracted the same way `ternary_grid`'s `EDGE`
  points are (see CLAUDE.md), just for a 2-of-4 subset instead of 2-of-3.
  There is no `EDGE`-equivalent region here: a square boundary (one of X,Y
  at 0 or 1) always zeroes two mole fractions at once, not one, so what
  would be an "edge" case in the triangle is already a binary point here.
- **`CORNER`** (three mole fractions = 0, i.e. a pure element) — no
  calculation; G/H/S are analytically 0, and no point is emitted for SRO.

`"skipped"` means the same thing as in `ternary_grid`: a real computation
was attempted (an interior or square-edge-binary solve) and failed to
converge, never "this quantity doesn't apply here."

### Getting rendered images, not just numbers

Add `"render": true` to have the server also render each of the two results
and return each as embedded base64 PNG — **one image per result entry**,
not one top-level image:

```json
{"system": {"elements":"Nb-Ti-V-Zr","structure":"BCC_A2","model":"T","engine":"CVM"},
 "calculation": "GIBBS_ENERGY", "temperature": 1273, "n": 50, "render": true}
```

```json
{"ok": true, "...": "...",
 "results": [
   {"slotOrder": "Nb-Ti-V-Zr", "...": "...",
    "image": {"format": "png", "base64": "iVBORw0KGgoAAAANS..."}},
   {"slotOrder": "Nb-Ti-Zr-V", "...": "...",
    "image": {"format": "png", "base64": "iVBORw0KGgoAAAANS..."}}
 ]}
```

Points are always included regardless of `render`. Rendering shells out to
`scripts/square_section.py` — plain matplotlib, **not** mpltern; unlike the
ternary triangle, the `(X,Y)` square is already Cartesian, so no special
ternary-axis plotting library is needed. Requires a working `python` on
`PATH` with `matplotlib`/`numpy` installed. A rendering failure is reported
per-result as `"renderError"` instead of `"image"`, without failing the
request or the other result.

```python
resp = call_quaternary_square({  # same call() pattern as above, subcommand "quaternary_square"
    "system": {"elements":"Nb-Ti-V-Zr","structure":"BCC_A2","model":"T","engine":"CVM"},
    "calculation": "GIBBS_ENERGY", "temperature": 1273, "n": 50, "render": True,
})
import base64
for r in resp["results"]:
    if "image" in r:
        with open(f"square_{r['slotOrder']}.png", "wb") as f:
            f.write(base64.b64decode(r["image"]["base64"]))
    elif "renderError" in r:
        print(f"render failed for {r['slotOrder']}, using raw points instead:", r["renderError"])
```
