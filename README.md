# CE Thermodynamics Workbench

Compute the thermodynamics of alloy solid solutions from first-principles
interaction parameters, using the **Cluster Variation Method** (CVM) or
**Monte Carlo** simulation.

Given a set of effective cluster interactions (ECIs) for an alloy system, the
workbench returns Gibbs energy, enthalpy, entropy, equilibrium correlation
functions, and Cowley-Warren short-range order parameters — for binary through
quaternary systems.

```bash
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"
```
```
  T           :         1000.0 K
  x           : [0.5000, 0.5000]
  G (J/mol)       :     -3480.5209063901
  H (J/mol)       :      2214.3770000563
  S (J/mol·K)     :         5.6948979064
```

---

## What it does

The workbench covers three kinds of task:

| | Task | Purpose |
|---|---|---|
| **Type-1a** | Cluster identification | Derive the cluster basis for a crystal structure — symmetry orbits, correlation functions, C-matrix, and the CVCF transformation |
| **Type-1b** | Hamiltonian scaffold | Generate an empty ECI file for a system, ready to fill in with fitted values |
| **Type-2** | Equilibrium calculation | Minimize free energy at given temperature and composition (CVM), or sample it (Monte Carlo) |

Most users only need Type-2: supply ECIs, get thermodynamics. Types 1a/1b exist
for setting up a system the workbench hasn't seen before.

## Why the CVM

Unlike a regular-solution or Redlich-Kister model, the CVM treats **short-range
order** explicitly — atoms are not randomly mixed. This matters for concentrated
solid solutions and high-entropy alloys, where SRO measurably changes the
entropy and stability. The workbench reports SRO directly, as Cowley-Warren α
parameters per neighbour shell.

The method follows Jindal & Lele, *Calphad* **89** (2025) 102825.

---

## Status and scope

**What works today**

- CVM equilibrium (tetrahedron approximation) for **BCC_A2**, binary through
  quaternary — the path that is regression-tested.
- Cluster identification for BCC_A2, BCC_B2, FCC_A1, HCP_A3 (see
  [USER_GUIDE](USER_GUIDE.md#supported-systems) for exact K coverage).
- Gibbs energy, enthalpy, entropy, correlation functions, Cowley-Warren SRO.
- Temperature and composition scans.
- Ternary isothermal-section and quaternary square-plot composition scans
  and rendered contour plots (GUI, CLI, JSON API).
- GUI, CLI, a Java API, and a JSON API for calling from other languages.

**Known limitations — read before relying on results**

- **Shipped ECI data is BCC_A2 only.** Four Hamiltonians are included (Nb-Ti,
  Nb-Zr, Nb-Ti-V, Nb-Ti-V-Zr). Other systems need ECIs supplied by you.
- **CVM minimization does not always converge**, particularly at low
  temperature and near-dilute compositions. Results carry a `converged` flag —
  check it. See [USER_GUIDE](USER_GUIDE.md#when-a-calculation-does-not-converge).
- **Monte Carlo is less mature than CVM**: it returns enthalpy and correlation
  functions, but heat capacity and error bars are not yet wired up.
- **No phase equilibria.** The workbench computes properties of a *single*
  phase at a given composition. Common-tangent construction, phase diagrams,
  and multi-phase equilibria are not implemented.
- Type-1b scaffolding does not yet produce a file the calculator can load
  directly — see [USER_GUIDE](USER_GUIDE.md#adding-a-new-system).

---

## Getting started

Requires **Java 21+**. Gradle is included via the wrapper.

```bash
./gradlew runGui        # graphical interface
./gradlew run --args="calc_min Nb-Ti BCC_A2 T 1000 Ti=0.5"    # command line
```

The **[User Guide](USER_GUIDE.md)** walks through both, explains the inputs, and
covers what to do when something fails.

---

## Documentation

| Document | For |
|---|---|
| **[USER_GUIDE.md](USER_GUIDE.md)** | Running calculations — GUI, CLI, inputs, interpreting output, troubleshooting |
| **[API.md](API.md)** | Calling the workbench from your own code (Java, Python, any language) |
| **[CLAUDE.md](CLAUDE.md)** | Internal architecture — for contributors and code-modifying agents |
| **[CHANGELOG.md](CHANGELOG.md)** | Development history and in-progress work |

---

## Citing

The CVM formulation, CVCF basis, and SRO definitions implemented here are from:

> V. Jindal and S. Lele, "Multicomponent cluster variation method: Application
> to high entropy alloys", *Calphad* **89** (2025) 102825.

---

## License

MIT — see [LICENSE](LICENSE).
