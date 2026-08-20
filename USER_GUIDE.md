# User Guide

How to run calculations with the CE Thermodynamics Workbench.

For calling the workbench from your own code, see [API.md](API.md).
For internal architecture, see [CLAUDE.md](CLAUDE.md).

**Contents**

- [Installation](#installation)
- [Core concepts](#core-concepts)
- [Supported systems](#supported-systems)
- [Using the GUI](#using-the-gui)
- [Using the command line](#using-the-command-line)
- [Specifying composition](#specifying-composition)
- [Reading the output](#reading-the-output)
- [Short-range order](#short-range-order)
- [Where your data lives](#where-your-data-lives)
- [Adding a new system](#adding-a-new-system)
- [Troubleshooting](#troubleshooting)

---

## Installation

Requires **Java 21 or later**. Nothing else — Gradle ships with the repository.

```bash
git clone https://github.com/vijindal/CEWorkbench.git
cd CEWorkbench
./gradlew build
```

Verify the install by reproducing a known value:

```bash
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"
```

You should get `G = -3480.5209063901 J/mol`. If you do, everything works.

---

## Core concepts

Four things identify a calculation.

**Elements** — the alloy system, hyphenated in a fixed order, e.g. `Nb-Ti` or
`Mo-Nb-Ta`. The order matters: it determines which element is `A`, `B`, `C` in
ECI names, and which one is derived when you omit a mole fraction.

**Structure** — the crystal lattice, e.g. `BCC_A2` (disordered body-centred
cubic), `FCC_A1`, `HCP_A3`, `BCC_B2` (ordered).

**Model** — the cluster approximation. `T` means the tetrahedron approximation:
the largest cluster retained is a 4-site tetrahedron. This sets how much
short-range order the model can represent.

**Engine** — how equilibrium is found:

- **CVM** minimizes the Gibbs energy analytically (Newton–Raphson). Fast,
  deterministic, gives G/H/S and SRO. This is the mature path.
- **MCS** samples configurations by Metropolis Monte Carlo on a supercell.
  Slower and statistical, but not limited by the cluster approximation's
  entropy expression.

Together these pick a **Hamiltonian**: the set of effective cluster interactions
(ECIs) describing how much energy each cluster configuration contributes. Without
ECIs for your system, no calculation is possible — see
[Adding a new system](#adding-a-new-system).

---

## Supported systems

**Cluster basis** (the mathematics is implemented):

| Structure | Model | Components |
|---|---|---|
| BCC_A2 | T | 2, 3, 4 |
| FCC_A1 | T | 2, 3, 4 |
| HCP_A3 | T | 2, 3, 4 |
| BCC_B2 | T | 2 |

**Shipped ECI data** (calculations work out of the box):

| System | Structure |
|---|---|
| Nb-Ti | BCC_A2 |
| Nb-Zr | BCC_A2 |
| Nb-Ti-V | BCC_A2 |
| Nb-Ti-V-Zr | BCC_A2 |

These two tables are different, and the difference matters. A structure being in
the first table means the workbench *can* handle it; only systems in the second
have interaction data included. For anything else you supply ECIs yourself.

---

## Using the GUI

```bash
./gradlew runGui
```

The activity bar on the left has five panels:

| Icon | Panel | Use |
|---|---|---|
| **1a** | Data Prep | Run cluster identification for a structure |
| **1b** | Hamilt. | Create, view, and edit ECI tables |
| **TH** | Thermo | Run calculations — the main panel |
| **TP** | Ternary | Ternary isothermal-section composition scans and plots |
| **QS** | Quat. Sq. | Quaternary square-plot composition scans (Fig. 20 reproduction) |

### Running a calculation

1. In the **Thermo** panel, set **Elements**, **Structure**, **Model**, and
   **Engine**. The fields are editable comboboxes — type a value that isn't in
   the dropdown if you need to.
2. Pick a **Property** (Gibbs Energy, Enthalpy, Entropy for CVM).
3. Fill in the parameter form. It rebuilds itself when you change the property
   or the element list, so composition rows always match your system.
4. Click **Run Calculation**.

There is no separate "build session" step — the session is constructed on demand
and reused for subsequent runs on the same system.

Each range row takes **Start**, **End**, and **Step**. Leave Start = End for a
single point. To scan, set them apart — but **only one axis may vary per run**
(either temperature or one element, not both).

Results appear in the output panel on the right, with a log below.

### Running a quaternary square-plot scan

1. In the **Quat. Sq.** panel, set **Elements** (exactly 4, e.g.
   `Nb-Ti-V-Zr`), **Structure**, and **Model**.
2. Set **Temperature** and pick a **Quantity** (Gibbs Energy, Enthalpy,
   Entropy, or SRO — SRO reveals a **Pair** dropdown).
3. Optionally adjust **Grid resolution (n)** (default 50).
4. Click **Compute & Plot**.

The panel always computes and shows **two** plots side by side — one for
each of the two square parametrizations needed to cover all six binary
edges of the quaternary composition tetrahedron. There is no manual variant
picker; this is deliberate.

---

## Using the command line

```bash
./gradlew run --args="<mode> <elements> <structure> <model> [options]"
```

| Mode | Purpose |
|---|---|
| `calc_min` | Single-point calculation at a given T and composition |
| `type2` | Temperature scan at equiatomic composition |
| `type1a` | Cluster identification only |
| `type1b` | Scaffold an empty Hamiltonian |
| `all` | `type1a` → `type1b` → `type2` |
| `view` | Print a Hamiltonian's ECI table |
| `api` | JSON in/out — see [API.md](API.md) |

Add `--verbose` to any mode for iteration-level detail.

> `ternary_grid` and `quaternary_square` are separate JSON stdin/stdout
> subcommands (not part of the `<mode> <elements> <structure> <model>`
> pattern above) for composition-grid scans — see [API.md](API.md).

### Examples

```bash
# Binary at 1000 K, equiatomic
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"

# Quaternary, asking for entropy specifically
./gradlew run --args="calc_min Nb-Ti-V-Zr BCC_A2 T 1273 Ti=0.25 V=0.25 Zr=0.25 S"

# Inspect the ECIs being used
./gradlew run --args="view Nb-Ti BCC_A2 T"
```

Append `G`, `H`, or `S` to `calc_min` to request a specific property. `G` is the
default and also reports H and S.

---

## Specifying composition

Composition is given as **`Element=fraction`** pairs. Never as bare numbers —
positional mole fractions are not accepted.

```bash
Ti=0.5                    # binary Nb-Ti: Nb derived as 0.5
Ti=0.33 V=0.34            # ternary Nb-Ti-V: Nb derived as 0.33
Ti=0.25 V=0.25 Zr=0.25    # quaternary: Nb derived as 0.25
```

**You may omit one element** — its fraction is derived so the total is 1. Omit
none and the values must already sum to 1. Omit two or more and you get an
`Underdetermined` error, because the composition isn't uniquely defined.

This is deliberate. An earlier positional convention silently produced wrong
answers for systems with three or more components, because there was no way to
tell which number belonged to which element. Naming them removes the ambiguity.

---

## Reading the output

```
  T           :         1000.0 K
  x           : [0.5000, 0.5000]
  G (J/mol)       :     -3480.5209063901
  H (J/mol)       :      2214.3770000563
  S (J/mol·K)     :         5.6948979064
```

| Quantity | Meaning |
|---|---|
| **G** | Gibbs energy of mixing, J/mol. Negative favours the solid solution over separated pure elements. |
| **H** | Enthalpy of mixing, J/mol. Negative = ordering tendency, positive = clustering. |
| **S** | Configurational entropy, J/(mol·K). Compare against the ideal value R·ln(K) for K components; the CVM result is lower because short-range order removes configurational freedom. |

Useful reference points: R·ln 2 = 5.763, R·ln 3 = 9.134, R·ln 4 = 11.526
J/(mol·K). An entropy close to these means a nearly random solution; well below
means significant SRO.

Correlation functions describe the equilibrium cluster populations. They are
mainly diagnostic — for a physical picture of ordering, use SRO instead.

---

## Short-range order

CVM calculations produce **Cowley-Warren SRO parameters**, which quantify how
far the atomic arrangement departs from random:

```
α = 1 − p(unlike pair) / (x_P · x_R)
```

| α | Meaning |
|---|---|
| **< 0** | Unlike neighbours *enriched* — ordering |
| **= 0** | Random solution |
| **> 0** | Unlike neighbours *depleted* — clustering / phase-separation tendency |

Reported per neighbour shell (`1NN`, `2NN`) and per element pair. Magnitudes
typically fall as temperature rises, since thermal disorder washes order out —
though in a multicomponent system an individual pair need not decrease
monotonically, because pairs compete for the same atoms.

SRO is currently exposed through the [JSON API](API.md#short-range-order) and
the Java API.

---

## Where your data lives

```
data/CEWorkbench/            (or ~/CEWorkbench/)
 ├─ inputs/
 │   ├─ clus/                cluster geometry, e.g. BCC_A2-T.txt
 │   └─ sym/                 symmetry operations, e.g. BCC_A2-SG.txt
 └─ hamiltonians/
     └─ Nb-Ti_BCC_A2_T_CVCF/
         └─ hamiltonian.json      ← ECI values
```

The workbench looks for this directory in order:

1. `$CEWORKBENCH_DATA` (or `-Dceworkbench.data=...`)
2. `./data/CEWorkbench` relative to where you launched it
3. `~/CEWorkbench`

**Set `CEWORKBENCH_DATA` if you run from anywhere other than the project
directory** — otherwise the inputs won't be found.

### The Hamiltonian file

```json
{
  "elements": "Nb-Ti",
  "structurePhase": "BCC_A2",
  "cecTerms": [
    { "name": "e21AB", "a": 6240.0, "b": 0.0 },
    { "name": "e22AB", "a": 3120.0, "b": 0.0 }
  ]
}
```

Each term is an ECI evaluated as **`a + b·T`** — put a temperature-independent
value in `a` and leave `b` at 0 unless you have a temperature-dependent fit.
Units are J/mol.

Term names identify clusters: `e21`/`e22` are 1st- and 2nd-nearest-neighbour
pairs, `e3` triangles, `e4` tetrahedra. The trailing letters are element pairs
by position — for `Nb-Ti-V`, `AB` is Nb-Ti, `AC` is Nb-V, `BC` is Ti-V.

**Names must match the basis exactly.** A name that doesn't match is not a
warning — the calculation refuses to run, because a silently-ignored ECI would
produce plausible but wrong energies. To see the exact names your system
expects, use the `describe` call in [API.md](API.md#discovering-what-a-system-expects).

---

## Adding a new system

To calculate a system with no shipped ECIs:

1. Confirm the structure and component count are in
   [Supported systems](#supported-systems).
2. Obtain ECIs — typically by fitting to first-principles energies with a
   cluster-expansion tool, then converting to the CVCF basis in J/mol.
3. Create `hamiltonians/<elements>_<structure>_<model>_CVCF/hamiltonian.json`
   with those values.

`type1b` will scaffold an empty ECI file for you:

```bash
./gradlew run --args="type1b Nb-Ti FCC_A1 T"
```

> **Known gap.** `type1b` writes the directory *without* the `_CVCF` suffix,
> but the calculator only loads `..._CVCF`. Rename the directory (and keep the
> term names) before calculating, or write the file directly at the `_CVCF`
> path.

Passing ECIs **inline** through the [JSON API](API.md) avoids the file layout
entirely, and is the better route when another program generates them.

---

## Troubleshooting

### "Required CVCF Hamiltonian file not found"

There are no ECIs for that system. Check the exact directory name — it must end
in `_CVCF`, e.g. `Nb-Ti_BCC_A2_T_CVCF`. See
[Adding a new system](#adding-a-new-system).

### "Failed to load cluster file"

The workbench can't find its `inputs/` directory. Set `CEWORKBENCH_DATA` to your
`data/CEWorkbench` path, or run from the project root.

### "Unknown element 'Xx' for system [Nb, Ti]"

A composition name doesn't match the element list. Spelling is flexible on case
but the symbol must be one of the system's elements.

### "Underdetermined: N elements unspecified"

You omitted more than one mole fraction. Specify all but one.

### "Only one condition axis may vary per scan"

Both temperature and a composition have differing Start and End. Fix one of them
to a single value. Multi-dimensional scans require repeated calls.

### "ECI_VALIDATION_FAILED"

Supplied ECI names don't line up with the basis. The error lists `unmatched`
(names not recognized), `unmapped` (basis entries you didn't supply), and
`expected` (the full correct list). Both problems are fatal by design — see
[the Hamiltonian file](#the-hamiltonian-file).

### When a calculation does not converge

The CVM minimizer is iterative and can fail, most often at **low temperature**
or **near-dilute compositions**. A failed run still returns numbers, and they
can look reasonable while being wrong — a common signature is entropy *rising*
as temperature falls, which is backwards.

The JSON API reports `converged` per point, with `iterations` and
`finalGradientNorm` when it fails. **Check it.** In the CLI, run with
`--verbose` and look for a non-convergence warning.

If a point won't converge:

- Move to higher temperature and confirm the trend is sensible there.
- Step in from a converged neighbouring point rather than jumping straight to
  the difficult condition.
- Treat compositions near a pure element with particular suspicion.

### Results look physically wrong

Sanity checks worth running before trusting a number:

- Entropy should be **below** R·ln(K) and should **increase** with temperature.
- Enthalpy of mixing should approach zero as temperature rises.
- Gibbs energy should **decrease** monotonically with temperature.
- At a pure element, G, H, and S should all be ≈ 0.

If any of these fail, suspect non-convergence first, then the ECI values.
