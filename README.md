# CE Thermodynamics Workbench

A scientific software framework for **Cluster Expansion (CE) based thermodynamic calculations** and **cluster identification**. This workbench provides tools to generate reusable scientific data and compute thermodynamic equilibrium states for alloy systems.

## Overview

The CE Thermodynamics Workbench implements two fundamental classes of work:

1. **Type-1: Generate reusable scientific data** → Cluster identification and correlation function analysis
2. **Type-2: Compute thermodynamic equilibrium** → Free energy minimization and Monte Carlo simulations

## Quick Start

### Prerequisites

- Java 25 or later
- Gradle 9.3+

### Build

```bash
./gradlew clean build
```

### Run

```bash
./gradlew run
```

This runs the complete CVM identification pipeline with the A2 system example.

### Run Tests

```bash
./gradlew test
```

## Architecture

The system is organized into **four clean layers** with strict dependency rules:

```
org.ce
 ├─ domain      → physics models and algorithms
 ├─ workflow    → calculation orchestration
 ├─ storage     → disk IO and resource management
 └─ ui          → GUI and CLI interfaces
```

### Dependency Rules (Golden Rule)

```
domain → (no dependencies on workflow, storage, or ui)
workflow → domain, storage
storage → domain
ui → workflow
```

### Domain Layer (Physics)

Scientific concepts and algorithms. **Never depends on other layers.**

**Subpackages:**

- **cluster** - Structural cluster information
  - `Cluster`, `ClusterData`, `AllClusterData`
  - `ClusterIdentifier`, `CFIdentifier`, `CMatrixBuilder`
  - `Vector3D`, `Position`, `Site`, `Sublattice`

- **engine** - Thermodynamic solvers
  - `ThermodynamicEngine`, `CVMEngine`, `MCSEngine`
  - Compute equilibrium and minimize free energy

- **result** - Result objects
  - `EquilibriumState`, `EngineMetrics`
  - Represent immutable thermodynamic results

### Workflow Layer (Orchestration)

Coordinates domain operations without implementing physics.

**Key classes:**

- `ClusterIdentificationWorkflow` - Type-1: generate cluster data
- `ClusterIdentificationRequest` - Configuration for cluster identification
- `CalculationService` - General calculation orchestration (planned)

### Storage Layer (IO)

Handles all disk IO and data persistence.

**Key classes:**

- `InputLoader` - Parse cluster files and symmetry groups
- `ClusterDataStore` - Persist cluster data (planned)
- `HamiltonianStore` - Persist Hamiltonian models (planned)
- `Workspace` - Manage workspace and file paths (planned)

### UI Layer (User Interface)

User-facing interfaces.

**Current:**

- **cli** - Command-line interface
  - `Main.java` - Entry point for complete pipeline

**Planned:**

- **gui** - Graphical interface with calculation panels

## Package Structure

```
app/src/main/java/org/ce
├─ domain
│  ├─ cluster
│  │  ├─ ClusterIdentifier.java
│  │  ├─ CFIdentifier.java
│  │  ├─ CMatrixBuilder.java
│  │  ├─ AllClusterData.java
│  │  ├─ Vector3D.java
│  │  └─ (40+ supporting classes)
│  ├─ engine
│  │  ├─ CVMEngine.java
│  │  ├─ MCSEngine.java
│  │  └─ ThermodynamicEngine.java
│  └─ result
│     └─ EquilibriumState.java
│
├─ workflow
│  ├─ ClusterIdentificationWorkflow.java
│  ├─ ClusterIdentificationRequest.java
│  └─ CalculationService.java
│
├─ storage
│  ├─ input
│  │  ├─ InputLoader.java
│  │  ├─ ClusterParser.java
│  │  └─ SpaceGroupParser.java
│  ├─ ClusterDataStore.java
│  ├─ HamiltonianStore.java
│  └─ Workspace.java
│
└─ ui
   └─ cli
      └─ Main.java
```

## Calculation Workflows

### Type-1: Cluster Identification (Stage 1-3)

Generate cluster and correlation function data:

```
ClusterIdentificationWorkflow
 ↓
1. Load disordered and ordered clusters
2. Load symmetry operations
 ↓
ClusterIdentifier → Stage 1: Identify clusters
 ↓
CFIdentifier → Stage 2: Identify correlation functions
 ↓
CMatrixBuilder → Stage 3: Build C-matrix
 ↓
AllClusterData (results bundled)
 ↓
ClusterDataStore.save()
```

**Example usage:**

```java
ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
    .disorderedClusterFile("clus/A2-T.txt")
    .orderedClusterFile("clus/A2-T.txt")
    .disorderedSymmetryGroup("A2-SG")
    .orderedSymmetryGroup("A2-SG")
    .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})
    .translationVector(new Vector3D(0, 0, 0))
    .numComponents(2)
    .build();

AllClusterData result = ClusterIdentificationWorkflow.identify(config);
```

### Type-2: Thermodynamic Equilibrium (Planned)

Compute equilibrium states given cluster data and Hamiltonian:

```
CalculationService
 ↓
ClusterDataStore.load()
HamiltonianStore.load()
 ↓
ThermodynamicEngine (CVMEngine or MCSEngine)
 ↓
EquilibriumState
 ↓
UI → Display results
```

## Resources

Input data files are located in `app/src/main/resources/`:

- **clus/** - Cluster coordinate files (Mathematica format)
  - A1-TO.txt, A2-T.txt, B2-T.txt

- **sym/** - Space group symmetry operations
  - A1-SG.txt, A2-SG.txt, B2-SG.txt

## Key Concepts

### Vector3D
Immutable 3D vector for coordinates and transformations.
- Fields: x, y, z (doubles)
- Methods: getters, distance(), equals(), hashCode()

### AllClusterData
Bundles cluster and CF identification results for both disordered and ordered phases.
- Maintains semantic distinction between phases
- Internally: each result object contains both disordered and ordered data
- Methods: getters, printResults(), getSummary()

### ClusterIdentificationResult
Contains complete cluster analysis for both phases.
- Includes Nij table, Kikuchi-Baker coefficients
- Stores cluster types, multiplicities, orbit information
- Separates HSP (disordered) and ordered phase metrics

## Testing

Run the test suite:

```bash
./gradlew test
```

Example test: `ClusterIdentificationWorkflowTest.testA2B2SystemWithBinaryBasis()`

Tests verify:
- Cluster identification produces correct results
- CF identification completes successfully
- C-matrix construction works properly
- Results are non-null and properly structured

## Design Principles

This architecture follows these principles:

1. **Physics-first** - Domain layer contains pure physics/algorithms
2. **Minimal abstractions** - No unnecessary factories, repositories, adapters
3. **Clean separation** - Each layer has clear, single responsibility
4. **Easy extension** - Add new workflows, engines, or UI components independently
5. **Scientific clarity** - Structure mirrors the actual thermodynamic workflow

## Building from Source

### Prerequisites

- Java 25 (JDK)
- Gradle 9.3.1 (included in repo)

### Development

Edit source files in `app/src/main/java/org/ce/`

Build and test:

```bash
./gradlew build
```

### Distribution

Create executable JAR:

```bash
./gradlew assemble
```

Output: `build/libs/CEWorkbench-1.0.jar`

## Project Status

✅ **Complete:**
- Layered architecture with clean dependencies
- Cluster identification (Stage 1)
- Correlation function identification (Stage 2)
- C-matrix construction (Stage 3)
- CLI interface

🔄 **In Progress:**
- Storage layer for persistence
- Thermodynamic engines (CVM, MCS)

📋 **Planned:**
- GUI interface
- Type-2 equilibrium calculations
- Additional symmetry groups and systems

## Contributing

The codebase is organized for easy extension:

- **New domain algorithms** → Add to `domain/` packages
- **New workflows** → Create in `workflow/` package
- **New calculations** → Add `ThermodynamicWorkflow` class
- **New UI** → Extend `ui/` with gui components

## References

- **Cluster Expansion Theory** - See comments in domain layer
- **Symmetry Operations** - Load from resource files in sym/
- **Mathematica Format** - ClusterParser handles .txt input format

## License

[Add your license here]

## Authors

- Project: CE Thermodynamics Workbench
- Refactored architecture: Claude
- Original algorithms: [Your names here]
