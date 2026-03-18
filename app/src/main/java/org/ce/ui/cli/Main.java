package org.ce.ui.cli;

import org.ce.domain.cluster.*;
import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.Workspace;
import org.ce.storage.HamiltonianStore;
import org.ce.storage.InputLoader;
import org.ce.workflow.CalculationService;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import java.util.List;

/**
 * Main entry point for the CE Workbench CLI.
 *
 * Usage:
 *   ./gradlew run --args="type1a"   -- cluster identification + save cluster_data.json
 *   ./gradlew run --args="type1b"   -- scaffold empty Hamiltonian from saved cluster data
 *   ./gradlew run --args="type2"    -- thermodynamic calculation (requires type1a + type1b)
 *   ./gradlew run                   -- full pipeline: type1a -> type1b -> type2
 */
public class Main {

    // cluster data ID:   {structure}_{model}_{ncomp}  — element-agnostic
    private static final String CLUSTER_ID     = "BCC_B2_T_bin";
    // hamiltonian ID:    {elements}_{structure}_{model}  — element-specific
    private static final String HAMILTONIAN_ID = "A-B_BCC_B2_T";

    public static void main(String[] args) {

        String mode = args.length > 0 ? args[0].toLowerCase() : "all";

        if (!mode.equals("type1a") && !mode.equals("type1b")
                && !mode.equals("type2") && !mode.equals("all")) {
            System.err.println("Unknown mode: " + mode);
            System.err.println("Valid modes: type1a  type1b  type2  all");
            System.exit(1);
        }

        System.out.println("================================================================================");
        System.out.println("                    CE THERMODYNAMICS WORKBENCH");
        System.out.println("================================================================================");
        System.out.println("Mode: " + mode);

        try {
            Workspace workspace         = new Workspace();
            ClusterDataStore clusterStore   = new ClusterDataStore(workspace);
            HamiltonianStore hamiltonianStore = new HamiltonianStore(workspace);
            CECManagementWorkflow cecWorkflow = new CECManagementWorkflow(hamiltonianStore, clusterStore);

            java.nio.file.Path inputsDir = workspace.inputsDir();

            // ------------------------------------------------------------------
            // TYPE-1a: Cluster Identification + Save
            // ------------------------------------------------------------------
            if (mode.equals("type1a") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1a: Cluster Identification ===\n");

                // Transformation matrix and translation vector are automatically extracted from symmetry group files
                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .disorderedClusterFile("clus/A2-T.txt")
                        .orderedClusterFile("clus/B2-T.txt")
                        .disorderedSymmetryGroup("A2-SG")
                        .orderedSymmetryGroup("B2-SG")
                        .numComponents(2)
                        .build();

                AllClusterData clusterData = ClusterIdentificationWorkflow.identify(config);
                System.out.println("Identification complete: " + clusterData.getSummary());

                List<Cluster> maxClusters = InputLoader.parseClusterFileFromPath(inputsDir, "clus/A2-T.txt");
                CMatrixResult cmatrix = CMatrixBuilder.build(
                        clusterData.getDisorderedClusterResult(),
                        clusterData.getDisorderedCFResult(),
                        maxClusters,
                        config.getNumComponents()
                );
                System.out.println("C-matrix built: " + cmatrix.getLcv().length + " cluster types");

                AllClusterData fullData = new AllClusterData(
                        clusterData.getDisorderedClusterResult(),
                        clusterData.getOrderedClusterResult(),
                        clusterData.getDisorderedCFResult(),
                        clusterData.getOrderedCFResult(),
                        cmatrix
                );
                clusterStore.save(CLUSTER_ID, fullData);
                System.out.println("Saved: ~/CEWorkbench/cluster-data/" + CLUSTER_ID + "/cluster_data.json");
            }

            // ------------------------------------------------------------------
            // TYPE-1b: Scaffold empty Hamiltonian from saved cluster data
            // ------------------------------------------------------------------
            if (mode.equals("type1b") || mode.equals("all")) {

                System.out.println("\n=== TYPE-1b: Hamiltonian Scaffold ===\n");

                if (hamiltonianStore.exists(HAMILTONIAN_ID)) {
                    System.out.println("Hamiltonian already exists: " + HAMILTONIAN_ID);
                    System.out.println("  Delete ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/ to re-scaffold.");
                } else {
                    cecWorkflow.scaffoldFromClusterData(HAMILTONIAN_ID, "A-B", "A2", "T");
                    System.out.println("Saved: ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/hamiltonian.json");
                    System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
                }
            }

            // ------------------------------------------------------------------
            // TYPE-2: Thermodynamic Calculation
            // ------------------------------------------------------------------
            if (mode.equals("type2") || mode.equals("all")) {

                System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");

                CVMEngine cvmEngine = new CVMEngine();
                ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
                        clusterStore, cecWorkflow, cvmEngine, null
                );
                CalculationService service = new CalculationService(thermoWorkflow);

                double temperature   = 1000.0;
                double[] composition = {0.5, 0.5};

                System.out.println("Running CVM at T=" + temperature + " K, x_B=" + composition[1]);
                ThermodynamicResult result = service.runSinglePoint(
                        CLUSTER_ID, HAMILTONIAN_ID, temperature, composition, "CVM"
                );

                System.out.println("\nResult:");
                System.out.println("  Temperature : " + result.temperature + " K");
                System.out.println("  Composition : x_B = " + result.composition[1]);
                System.out.println("  Free energy : " + String.format("%.6f", result.gibbsEnergy) + " J/mol");
                System.out.println("  Enthalpy    : " + String.format("%.6f", result.enthalpy) + " J/mol");
            }

        } catch (Exception e) {
            System.err.println("\nError: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("\n================================================================================");
    }
}
