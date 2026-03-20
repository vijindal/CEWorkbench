# CVM Thermodynamic Properties: Enthalpy, Entropy, and Gibbs Energy

## The Complete Thermodynamic Calculation

After entropy (Sc), the CVM engine computes enthalpy (Hc) and Gibbs free energy (Gc):

```mathematica
(* Configurational Enthalpy *)
Hc = Sum[
  Simplify[
    mhdis[[i]] * eList[[i]] * (uList[[i]] - (1-xB)*uA[[i]] - xB*uB[[i]])
  ],
  {i, tcf}
];

(* Apply substitutions for point correlations *)
Hc = Hc /. {u[tcdis][1][1] -> (2*xB - 1)} /. {u[tcdis+1][1][1] -> 1} // Simplify;

(* Gibbs Free Energy *)
Gc = Hc - T*Sc;
```

This completes the thermodynamic model for the given phase at any composition and temperature.

---

## Component 1: Configurational Enthalpy (Hc)

### Structure of the Enthalpy Expression

```mathematica
Hc = Sum[
  mhdis[[i]] * eList[[i]] * (uList[[i]] - (1-xB)*uA[[i]] - xB*uB[[i]])
],
  {i, 1, tcf}
]
```

This sums **over all correlation functions** (i = 1 to tcf).

---

## Data Sources: Where Each Variable Comes From

### **mhdis[[i]]** — HSP Multiplicity

**From:** Stage 1 (HSP Clusters)
**Source class:** `ClusterIdentificationResult.disorderedClusterResult.multiplicities`
**Role:** Weighting factor for each ECI contribution

- One multiplicity per correlation function
- Derived from cluster topology and symmetry
- Accounts for orbital averaging

Example: `mhdis = {6, 12, 24, 8, ...}` (different symmetric environment counts)

---

### **eList[[i]]** — Effective Cluster Interaction (ECI)

**From:** External input (fitted CE model)
**Source class:** `HamiltonianStore.load()` → `CEHamiltonian.eciParameters`
**Role:** Energy coefficient for each CF

- One ECI per correlation function
- Fitted from DFT or other first-principles calculations
- **Critical input from CE fitting procedure**
- Units: energy (typically eV or meV)

Example: `eList = {-0.5, 0.02, -0.01, 0.003, ...}` (from CE training)

**Relationship:** The Hamiltonian stores these parameters:
```java
CEHamiltonian {
  eciParameters: [e1, e2, e3, ...]  // One per CF
  cluster: ClusterData (defines tcf)
  ...
}
```

---

### **uList[[i]]** — Correlation Function Value (Runtime Variable)

**From:** Runtime computation (CVM solver)
**Computed from:** Occupation variables `{σᵢ, σᵢσⱼ, ...}`
**Role:** Actual CF value at this thermodynamic state

- One value per CF
- Computed via C-matrix projection: `uList = cmat · occupationProbs`
- Changes at each CVM iteration
- Range: typically [-1, 1] in conventional basis

**Computation pathway:**
```
Occupation variables {P₁, P₁₁, P₁₂, ...}
  ↓ [via C-matrix]
  ↓ cmat · P = uList
Correlation values {u[1,1,1], u[1,1,2], ..., u[i,j,k], ...}
```

---

### **uA[[i]] and uB[[i]]** — Pure-Phase Reference Correlations

**From:** Pre-computed reference states
**Role:** Baseline for enthalpy relative to pure phases

- One pair (uA, uB) per CF
- Computed once at composition boundaries
- uA = CF value in pure A phase
- uB = CF value in pure B phase

**Composition dependence:**
```
Reference CF = (1-xB) * uA[[i]] + xB * uB[[i]]
```

This is a linear interpolation between pure-phase CF values.

**Why needed:** Enthalpy is computed relative to a reference state:
```
Hc = Sum[mhdis[[i]] * eList[[i]] * (u_actual - u_reference)]
```

Without subtracting the reference, you'd include the enthalpy of the pure phases, which isn't desired in a phase diagram.

---

### **xB** — Mole Fraction of Component B

**From:** External input (thermodynamic calculation parameter)
**Role:** Composition variable

- For binary system: xB ∈ [0, 1]
- xA = 1 - xB
- Variable for which the phase diagram is computed

---

### **T** — Temperature

**From:** External input (thermodynamic calculation parameter)
**Role:** Temperature variable

- Kelvin (or whatever unit was used in ECI fitting)
- Variable for which the phase diagram is computed

---

## Component 2: Point Correlation Substitutions

After computing Hc symbolically, specific correlations are substituted:

```mathematica
Hc = Hc /. {u[tcdis][1][1] -> (2*xB - 1)}
        /. {u[tcdis+1][1][1] -> 1}
        // Simplify;
```

### **u[tcdis][1][1]** — Point Correlation of Given Phase

**Indexing:** `[tcdis][1][1]`
- `tcdis` = total number of HSP cluster types (the last type is the point cluster)
- `[1]` = first (only) point cluster in that class
- `[1]` = first (only) CF for the point cluster

**Physical meaning:** Site occupation (concentration variable)

**Substitution rule:**
```
u[tcdis][1][1] = 2*xB - 1
```

**Why this formula?** In conventional basis:
- Pure A phase (xB = 0): u = -1 (all A atoms)
- Equimolar (xB = 0.5): u = 0 (equal A and B)
- Pure B phase (xB = 1): u = +1 (all B atoms)

This is the **constraint** that fixes occupation to match composition.

---

### **u[tcdis+1][1][1]** — Reference Point Correlation

**Indexing:** `[tcdis+1]` — fictitious cluster type beyond the highest
**Substitution rule:**
```
u[tcdis+1][1][1] = 1
```

**Purpose:** Numerical reference point for free energy calculations

---

## Component 3: Gibbs Free Energy (Gc)

### The Fundamental Equation

```mathematica
Gc = Hc - T*Sc;
```

Where:
- **Hc** — configurational enthalpy (from CFs and ECI)
- **T** — temperature
- **Sc** — configurational entropy (from CVM cluster variation)

### Thermodynamic Interpretation

**Gibbs free energy determines phase stability:**
- Lower Gc → more stable phase
- At equilibrium between phases: Gc is equal
- Phase diagram points: where phases of different structures have equal Gc

---

## Data Integration Diagram: Full Thermodynamic Calculation

```
Type-1a Output (cluster_data.json)
├─ mhdis[i]                    ← Stage 1 (HSP multiplicities)
├─ kbdis[i], nijTable          ← Stage 1 (KB coefficients)
├─ mh[itc][inc]                ← Stage 2 (given-phase multiplicities)
├─ lc[itc]                     ← Stage 2 (cluster counts)
├─ lcf[itc][inc]               ← Stage 4 (CF counts per cluster)
├─ tcf = total CFs             ← Stage 4 (total CF count)
├─ cmat, lcv, wcv              ← C-Matrix (projection matrix)
└─ cv                          ← C-Matrix (grouped correlations)

Type-2 Runtime Input
├─ eList[i]                    ← External CE model (fitted ECI)
├─ xB                          ← Composition input
└─ T                           ← Temperature input

Type-2 CVM Solver
├─ Solve Lagrange multipliers
├─ Find occupation variables
└─ Compute:
   ├─ uList[i] = cmat · occupations          (CF values)
   ├─ uA[i], uB[i] = boundary correlations   (reference states)
   │
   ├─ Sc = -R * Sum[kbdis × mhdis × mh × wcv × uList × Log[uList]]
   │        (entropy from cluster data)
   │
   ├─ Hc = Sum[mhdis × eList × (uList - reference)]
   │        (enthalpy from CE model)
   │
   └─ Gc = Hc - T*Sc
            (Gibbs energy)

Output: Thermodynamic Properties
├─ Gc (Gibbs free energy)
├─ Hc (enthalpy)
├─ Sc (entropy)
└─ Other derived properties (heat capacity, etc.)
```

---

## Runtime Flow in Java

The thermodynamic calculation maps to:

```java
// CVM Engine (Type-2)
public EquilibriumState calculateEquilibrium(
    ClusterData clusterData,        // From cluster_data.json
    CEHamiltonian hamiltonian,      // External CE model
    double composition,              // xB
    double temperature              // T
) {
    // 1. Load cluster data
    double[] mhdis = clusterData.getMultiplicities();           // Stage 1
    double[] kbdis = clusterData.getKBCoefficients();           // Stage 1
    double[][] mh = clusterData.getMultiplicityRatios();        // Stage 2
    int[][] lc = clusterData.getClusterCounts();                // Stage 2
    Matrix cmat = clusterData.getCMatrix();                     // C-Matrix

    // 2. Load CE model
    double[] eList = hamiltonian.getECIParameters();            // External input

    // 3. Solve for occupation variables via Lagrange multipliers
    double[] occupations = solveCVM(clusterData, composition);

    // 4. Compute correlations
    double[] uList = cmat.multiply(occupations);                // u[i] = cmat·occ
    double[] uA = computeReferenceState(clusterData, false);    // Pure A
    double[] uB = computeReferenceState(clusterData, true);     // Pure B

    // 5. Compute entropy
    double entropy = computeEntropy(
        kbdis, mhdis, mh, clusterData.getWcv(),
        uList, clusterData.getLcv()
    );

    // 6. Compute enthalpy
    double enthalpy = 0.0;
    for (int i = 0; i < eList.length; i++) {
        double reference = (1 - composition) * uA[i] + composition * uB[i];
        enthalpy += mhdis[i] * eList[i] * (uList[i] - reference);
    }

    // 7. Compute Gibbs energy
    double gibbsEnergy = enthalpy - temperature * entropy;

    // 8. Return thermodynamic state
    return new EquilibriumState(
        gibbsEnergy, enthalpy, entropy, uList, composition, temperature
    );
}
```

---

## Why All Cluster Data Elements Are Essential

| Element | Hc | Sc | Gc | Why |
|---------|----|----|----|----|
| `mhdis[i]` | ✓ | ✓ | ✓ | Weighting factor in both H and S |
| `kbdis[i]` | - | ✓ | ✓ | KB coefficient in entropy (always × mhdis) |
| `mh[itc][inc]` | - | ✓ | ✓ | Cluster multiplicity scaling in S |
| `eList[i]` | ✓ | - | ✓ | ECI energy coefficients (external) |
| `cmat` | ✓ | ✓ | ✓ | Converts occupations to CFs |
| `uList[i]` | ✓ | ✓ | ✓ | CF values (computed at runtime) |
| `uA[i]`, `uB[i]` | ✓ | - | ✓ | Reference correlations for H |
| `T` | - | - | ✓ | Temperature for G = H - TS |

---

## Numerical Example

**System:** Binary B2 compound at xB = 0.5, T = 500 K

**Input:**
- `mhdis = {6, 12, 8}`
- `eList = {-0.5, 0.02, -0.01}` (eV)
- Solved: `uList = {1.0, 0.2, -0.1}` (from CVM)
- Reference: `uA = {-1, -0.8, -0.5}`, `uB = {1, 0.4, 0.3}`

**Enthalpy calculation:**

```
For i=1: mhdis[1] * eList[1] * (uList[1] - reference[1])
       = 6 * (-0.5) * (1.0 - 0) = -3.0 eV
         (reference at xB=0.5: (1-0.5)*(-1) + 0.5*(1) = 0)

For i=2: 12 * 0.02 * (0.2 - (-0.2)) = 0.096 eV
         (reference: (1-0.5)*(-0.8) + 0.5*(0.4) = -0.2)

For i=3: 8 * (-0.01) * (-0.1 - (-0.1)) = 0 eV
         (reference: (1-0.5)*(-0.5) + 0.5*(0.3) = -0.1)

Hc = -3.0 + 0.096 + 0 = -2.904 eV (per formula unit)
```

**Entropy calculation:** (from previous entropy doc)
```
Sc = 0.082 eV/K (example value from CVM)
```

**Gibbs energy:**
```
Gc = Hc - T*Sc = -2.904 - 500*0.082 = -2.904 - 41 = -43.904 eV
```

This Gc value determines the phase stability at T=500K, xB=0.5.

---

## Verification Checklist

After CVM calculation:

- [ ] `mhdis.length == tcf` — one multiplicity per CF
- [ ] `eList.length == tcf` — one ECI per CF (from CE model)
- [ ] `uList.length == tcf` — CF values computed correctly
- [ ] `uA.length == tcf`, `uB.length == tcf` — reference states complete
- [ ] Hc computation includes all ECI terms
- [ ] Point correlation substitution applied: `u[tcdis][1][1]` set to `2*xB - 1`
- [ ] Sc computed from entropy formula (from previous doc)
- [ ] Gc = Hc - T*Sc (no sign errors)
- [ ] Gc values sensible: typically large negative values in eV

---

## Complete Thermodynamic Properties Calculation Chain

```
Input: Cluster data + CE model + (xB, T)
  ↓
Compute occupations via CVM solver (Lagrange multipliers)
  ↓
uList = cmat · occupations  (CF values from C-matrix)
  ↓
Sc = entropy calculation    (using mhdis, kbdis, mh, lcv, wcv, uList)
  ↓
Hc = enthalpy calculation   (using mhdis, eList, uList, reference CFs)
  ↓
Gc = Hc - T·Sc              (Gibbs free energy)
  ↓
Phase stability determined by Gc
  → Lower Gc = more stable at this (xB, T)
  → Equal Gc across phases = equilibrium coexistence
```

---

## References

- [CLUSTER_DATA_STRUCTURE.md](CLUSTER_DATA_STRUCTURE.md) — source of mhdis, kbdis
- [CLUSTER_DATA_GENERATION_FLOW.md](CLUSTER_DATA_GENERATION_FLOW.md) — how cluster data generated
- [CMATRIX_CALCULATION_FLOW.md](CMATRIX_CALCULATION_FLOW.md) — source of cmat, lcv, wcv
- [CVM_ENTROPY_CALCULATION.md](CVM_ENTROPY_CALCULATION.md) — Sc calculation details
- [COMPLETE_DATA_FLOW.md](COMPLETE_DATA_FLOW.md) — end-to-end flow diagram
- [SIMPLIFIED_ARCHITECTURE.md](SIMPLIFIED_ARCHITECTURE.md) — CVM engine architecture
