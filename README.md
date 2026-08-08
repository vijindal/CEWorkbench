# CE Thermodynamics Workbench

A scientific software framework for **Cluster Expansion (CE) based thermodynamic calculations** and **cluster identification**. The workbench provides tools to identify cluster basis functions, manage effective cluster interactions (ECI), and compute thermodynamic equilibrium states for alloy systems.

## Overview

Three classes of work are supported:

| Type | Name | Description |
|------|------|-------------|
| **1a** | Cluster Identification | Load ordered/disordered cluster files and symmetry groups. Unified pipeline handles geometric identification, C-matrix construction, and CVCF transformation in one sweep. |
| **1b** | Hamiltonian Scaffold | Auto-generate an empty ECI (Hamiltonian) JSON file from consolidated pipeline results, ready for editing. |
| **2** | Thermodynamic Equilibrium | Minimize free energy with CVM (Newton–Raphson, CVCF basis) or Monte Carlo (MCS) via a metadata-driven calculation interface. |

---

## Quick Start

### Prerequisites

- Java 21 or later
- Gradle 9.3+

### Launch the GUI

```bash
./gradlew runGui
```

Opens the VS Code-style dark workbench. Use the activity bar on the left to switch between the three panels.

### Run the CLI

```bash
# Full pipeline — explicit system
./gradlew run --args="all Nb-Ti BCC_A2 T"

# Individual stages
./gradlew run --args="type1a Nb-Ti BCC_A2 T"
./gradlew run --args="type1b Nb-Ti BCC_A2 T"
./gradlew run --args="type2  Nb-Ti BCC_A2 T"

# Single-point CVM with minimisation — composition as <El>=<x> pairs.
# Any one element may be omitted; its fraction is derived as 1 - sum(given).
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"
./gradlew run --args="calc_min Nb-Ti-V-Zr BCC_A2 T 1273 Ti=0.25 V=0.25 Zr=0.25 S"

# View loaded Hamiltonian
./gradlew run --args="view Nb-Ti BCC_A2 T"

# Verbose output for any mode
./gradlew run --args="type2 Nb-Ti BCC_A2 T --verbose"
```

**Mode reference:**

| Mode | Description |
|------|-------------|
| `type1a` | Cluster identification only |
| `type1b` | Hamiltonian scaffold only |
| `type2` | Thermodynamic calculation (temperature scan) |
| `all` | Runs type1a → type1b → type2 in sequence |
| `calc_min` | Single-point CVM with Newton–Raphson minimisation |
| `view` | Print Hamiltonian ECI table to stdout |

### Build

```bash
./gradlew build
```

---

## Using as a library

`CalculationService.calculate` / `calculateScan` are the public API — the same
entry point the CLI and GUI use internally. Composition is specified by element
name (never by array position), following pycalphad's `conditions` convention.
Any one element may be omitted; its mole fraction is derived so the composition
sums to 1.

```java
CEWorkbenchContext ctx = new CEWorkbenchContext();
CalculationService service = ctx.getCalculationService();

// Build the session once per (elements, structure, model, engine) identity.
ModelSpecifications specs =
        new ModelSpecifications("Nb-Ti", "BCC_A2", "T", EngineConfig.CVM);
ModelSession session = service.getOrBuildSession(specs, null);

// Single point — Nb is derived as 1 - 0.5.
ThermodynamicResult r = service.calculate(
        session,
        new Conditions(1000.0, Map.of("Ti", 0.5)),
        Property.GIBBS_ENERGY,
        null, null);
System.out.println(r.gibbsEnergy);

// Temperature scan (at most one axis — T or one element — may vary).
List<ThermodynamicResult> scan = service.calculateScan(
        session,
        new ConditionsScan(new Range(800, 1200, 100), Map.of("Ti", Range.fixed(0.5))),
        Property.GIBBS_ENERGY,
        null, null);
```

For MCS, pass a `CalculationService.McsParams(L, nEquil, nAvg)` as the fourth
argument; CVM callers can omit it.

---

## Calling from other languages — the JSON API

Non-JVM callers (Python, MATLAB, shell) use the `api` subcommand: one JSON
request on stdin, one JSON response on stdout. Diagnostics go to **stderr**, so
the payload is always machine-parseable.

```bash
./gradlew installDist                       # once — builds the launcher
build/install/CEWorkbench/bin/CEWorkbench api < request.json
```

Exit code is 0 on success, 1 on error (a JSON error object is still written, so
the reason is always machine-readable).

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

Query this rather than hard-coding ECI names — it is the reliable guard against
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
`{"start":…,"end":…,"step":…}` range; at most one axis may vary per request.
Composition follows the same rule as everywhere else — name your elements, omit
at most one, and it is derived.

**Supplying your own ECIs.** Omit `"hamiltonian"` to use the stored CEC
database, or include it to pass ECIs directly:

```json
"hamiltonian": {
  "basis": "CVCF", "units": "J/mol",
  "cecTerms": [
    {"name":"e4AB","a":0.0,"b":0.0},   {"name":"e3AB","a":0.0,"b":0.0},
    {"name":"e22AB","a":3120.0,"b":0.0},{"name":"e21AB","a":6240.0,"b":0.0}
  ]
}
```

ECIs must **already be in the CVCF basis, in J/mol** — no basis transformation
or unit conversion is applied. `basis` and `units` are declared explicitly and
anything else is rejected, so a mismatch fails loudly instead of silently
producing wrong energies. Each term is `a + b·T`.

> **Names are validated strictly.** Every supplied name must resolve against the
> basis, and every basis CF must be supplied. An unmatched name would otherwise
> leave that interaction at 0.0 and still return a plausible-looking number, so
> the API refuses to compute:
>
> ```json
> {"ok":false,"error":"ECI_VALIDATION_FAILED",
>  "unmatched":["e21XX"], "unmapped":["v21AB"],
>  "expected":["v4AB","v3AB","v22AB","v21AB"]}
> ```
>
> Accepted spellings per name: canonical `v…`, the `e…` alias, the published
> `e2<species><shell>` pair notation, and `CF_<index>`.

### Python example

```python
import json, os, subprocess

# Use an absolute path; on Windows the launcher is the .bat variant.
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
                        {"name":"e4AB","a":0.0},   {"name":"e3AB","a":0.0},
                        {"name":"e22AB","a":3120.0},{"name":"e21AB","a":6240.0}]},
    "calculation": "GIBBS_ENERGY",
    "conditions":  {"temperature":{"start":800,"end":1200,"step":100},
                    "composition":{"Ti":0.5}},
})

for pt in resp["points"]:
    print(pt["temperature"], pt["gibbsEnergy"])
```

Each response point carries `temperature`, `composition` (element → fraction),
the available quantities, `correlationFunctions` keyed by CF name, `sro` (below),
and `converged`. Quantities an engine doesn't produce are **omitted** rather than
emitted as `NaN`, which is not valid JSON.

> **Always check `converged`.** The CVM minimizer can hit its iteration limit and
> still return plausible-looking numbers. When `converged` is `false` the point
> also carries `iterations` and `finalGradientNorm`; treat the values as invalid.

### Short-range order

CVM points include Cowley-Warren SRO parameters (Jindal &amp; Lele, *Calphad* 89
(2025) 102825, Eq. 40), computed from the converged cluster probabilities:

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

`1NN` and `2NN` are the first- and second-nearest-neighbour shells. All species
pairs are reported, like-pairs included.

### Workspace data

CEWorkbench needs its `inputs/` data (cluster and symmetry files). When running
from another directory, point it at the data folder explicitly:

```bash
export CEWORKBENCH_DATA=/path/to/CEWorkbench/data/CEWorkbench
```

Otherwise it looks for `./data/CEWorkbench` relative to the working directory,
then `~/CEWorkbench`.

---

## Architecture

Three layers with strict one-way dependencies (`ui` → `calculation` → `model`),
plus a root-level application context that wires them together:

```
org.ce
 ├─ CEWorkbench.java          GUI entry point
 ├─ CEWorkbenchContext.java   Wires all layers; shared by GUI and CLI
 │
 ├─ model/       Physics evaluators, optimizers, cluster data, Hamiltonian, storage I/O
 ├─ calculation/ Public API, vocabulary, discovery, dispatch, statistics
 └─ ui/          GUI (Swing), CLI, and the JSON API
```

`ModelSession` is the central object: immutable, built once per
(elements, structure, model, engine) identity, and passed as the first argument to
every calculation.

**For the full architecture** — package layout, layer roles, the session and
calculation-API contracts, dataflow for cluster identification / CVM / MCS, and the
invariants to preserve when changing things — see [CLAUDE.md](CLAUDE.md). It is kept
current against the code; this README covers usage.

---

## Calculation Workflows

| Type | What it does | Entry point |
|---|---|---|
| **1a** | Cluster identification (geometry → CF basis → C-matrix → CVCF) | `type1a` CLI mode, Data Prep panel |
| **1b** | Scaffold an empty Hamiltonian JSON from the identification results | `type1b` CLI mode, Hamiltonian panel |
| **2** | Thermodynamic equilibrium (CVM or MCS) | `calc_min` / `type2` CLI, `api`, or `CalculationService.calculate` |

Type-1b writes `hamiltonians/<hamiltonianId>/hamiltonian.json` with all ECI terms at
`a = 0, b = 0`, ready to edit. Type-2 reads it back (or takes ECIs inline through the
JSON API).


## Input Data Files

Located in `~/CEWorkbench/inputs/` (runtime workspace). Defaults are bundled in `src/main/resources/` and copied on first run.

| File | Description |
|------|-------------|
| `clus/BCC_A2-T.txt` | BCC A2 (disordered) cluster coordinates |
| `clus/BCC_B2-T.txt` | BCC B2 (ordered) cluster coordinates |
| `clus/FCC_A1-TO.txt` | FCC A1 cluster coordinates |
| `sym/BCC_A2-SG.txt` | BCC A2 space group symmetry operations |
| `sym/BCC_B2-SG.txt` | BCC B2 space group symmetry operations |
| `sym/FCC_A1-SG.txt` | FCC A1 space group symmetry operations |
| `sym/HCP_A3-SG.txt` | HCP A3 space group symmetry operations |

Cluster files follow the naming convention `<structure>-<model>.txt`; symmetry files follow `<structure>-SG.txt`. The `ModelSession.Builder` derives file paths automatically from the system identity.

---

## Workspace Layout

All persistent data is stored under `~/CEWorkbench/`:

```
~/CEWorkbench/
 ├─ inputs/
 │   ├─ clus/      cluster coordinate files (*.txt)
 │   └─ sym/       symmetry group files (*.txt)
 └─ hamiltonians/
     └─ <hamiltonianId>/
         └─ hamiltonian.json    edit a and b ECI values here
```

`SystemId` derives IDs deterministically from elements / structure / model:
- `hamiltonianId` = `<elements>_<structure>_<model>` e.g. `Nb-Ti_BCC_A2_T`
- For CVM, the builder prefers a `_CVCF`-suffixed Hamiltonian (`Nb-Ti_BCC_A2_T_CVCF`) when available.

---

## Project Status

**Complete:**
- Three-layer architecture (Model / Calculation / UI) with one-way dependencies
- `ModelSession` — immutable session object; consolidated Stages 1–4 pipeline runs once per system identity, reused across all scan points
- Unified `ClusterCFIdentificationPipeline` — Stages 1–4 (cluster, CF, C-matrix, CVCF transformation) merged into a single orchestration flow
- Full elimination of disk-based cluster data stores in favor of on-the-fly identification
- CVM thermodynamic engine with Newton–Raphson core
- MCS thermodynamic engine with automatic ECI basis transformation
- Hamiltonian scaffold, load, edit, save workflow
- Modern VS Code-style GUI — Metadata-driven `DynamicCalculationPanel`, Activity Bar, Explorer, Output panel
- Bidirectional system identity sync across all GUI panels via `WorkbenchContext`
- CLI with all modes: `type1a`, `type1b`, `type2`, `all`, `calc_min`, `view`

**Planned:**
- Additional symmetry groups and crystal structures
- Phase diagram grid scan visualization
- Export results to CSV / JSON

---

## License

MIT — see [LICENSE](LICENSE).
