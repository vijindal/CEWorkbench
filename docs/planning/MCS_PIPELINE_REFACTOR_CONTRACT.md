# Architecture-Driven MCS Pipeline Refactoring

This plan aligns the CEWorkbench calculation pipeline with the strict three-layer contract: UI harvests specifications, Calculation Layer owns persistent models and queues jobs, and the Model Layer validates and executes physics.

## Context & Objectives

The goal is to simplify the MCS calculation pipeline by:
1.  **Strictly classifying parameters**: Moving physical parameters (like supercell size `L`) that require a model rebuild into `ModelSpecifications`.
2.  **Clear Ownership**: `CalculationService` becomes the authoritative owner of the persistent model instance.
3.  **Queueing**: Implementing a synchronized execution queue in the Calculation Layer.
4.  **Streaming Results**: Providing point-by-point thermodynamic feedback to the UI while deferred plotting occurs at the end of a scan.

## 1. Specification Redesign (`org.ce.calculation`)

### [MODIFY] [CalculationDescriptor.java](file:///c:/Users/admin/codes/CEWorkbench/src/main/java/org/ce/calculation/CalculationDescriptor.java)
- Relocate `Parameter.MCS_L` from `JobSpecifications` to `ModelSpecifications`.
- Update `Registry.getRequirements` logic to ensure `L` is treated as a construction-time requirement for MCS models.

## 2. Model Layer Autonomy (`org.ce.model`)

### [NEW] [ThermodynamicModel.java](file:///c:/Users/admin/codes/CEWorkbench/src/main/java/org/ce/model/ThermodynamicModel.java)
A unified domain interface for physics engines.
- `void validate(ModelSpecifications specs)`: Internal validation for physical consistency (e.g. $L \ge 1$).
- `void calculate(T, x, JobSpecifications specs, Consumer<ThermodynamicResult> resultSink)`: Performs the actual physics and streams results.

### [NEW] `McsModel.java` & `CvmModel.java`
- Concrete implementations of `ThermodynamicModel`.
- `McsModel` encapsulates `MCSGeometry` and ensures it persists for the life of the engine instance.

## 3. Authorized Ownership & Queueing (`org.ce.calculation.workflow`)

### [MODIFY] [CalculationService.java](file:///c:/Users/admin/codes/CEWorkbench/src/main/java/org/ce/calculation/workflow/CalculationService.java)
- **Active Model Management**: Hold `private ThermodynamicModel activeModel`.
- **Job Queueing**: Implement a serialized runner (e.g. `singleThreadExecutor`) to prevent overlapping calculations on a single persistent engine.

## 4. Result Streaming & Plotting (`org.ce.calculation.workflow.thermo`)

### [MODIFY] [ThermodynamicWorkflow.java](file:///c:/Users/admin/codes/CEWorkbench/src/main/java/org/ce/calculation/workflow/thermo/ThermodynamicWorkflow.java)
- Update implementation to stream results point-by-point via the `resultSink`.
- Append a completion signal to notify the UI to finalize plots and statistical summaries.

## Verification Plan

### Automated
- **Queueing Verification**: Concurrent request testing to ensure serial execution.
- **Contract Verification**: Unit tests for `ModelSpecification` validation within the engines.

### Manual
- **GUI Feedback**: Monitoring the "Building..." vs "Running..." states in the terminal.
- **Scan Visibility**: Verifying that chart points appear incrementally during a temperature scan.
