package org.ce.ui.cli;

import org.ce.domain.cluster.*;
import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.SystemId;
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
 *   ./gradlew run                                  -- full pipeline, built-in defaults
 *   ./gradlew run --args="type1a"                  -- cluster identification only
 *   ./gradlew run --args="type1b"                  -- Hamiltonian scaffold only
 *   ./gradlew run --args="type2"                   -- thermodynamic calculation only
 *   ./gradlew run --args="all  Nb-Ti BCC_A2 T"     -- explicit system, all modes
 *   ./gradlew run --args="type2 Nb-Ti BCC_A2 T"    -- explicit system, type-2 only
 *
 * Arguments:  mode  [elements  structure  model]
 *   mode       : type1a | type1b | type2 | all  (default: all)
 *   elements   : element pair, e.g. Nb-Ti  (default: A-B)
 *   structure  : structure ID,  e.g. BCC_A2 (default: BCC_B2)
 *   model      : model ID,      e.g. T      (default: T)
 */
public class Main {

    public static void main(String[] args) {

        // Parse positional arguments — same three values as the GUI's system context
        String mode      = args.length > 0 ? args[0].toLowerCase() : "all";
        String elements  = args.length > 1 ? args[1] : "A-B";
        String structure = args.length > 2 ? args[2] : "BCC_B2";
        String model     = args.length > 3 ? args[3] : "T";

        if (!mode.equals("type1a") && !mode.equals("type1b")
                && !mode.equals("type2") && !mode.equals("all")) {
            System.err.println("Unknown mode: " + mode);
            System.err.println("Usage: [mode] [elements] [structure] [model]");
            System.err.println("  mode: type1a | type1b | type2 | all");
            System.exit(1);
        }

        SystemId system       = new SystemId(elements, structure, model);
        String   CLUSTER_ID   = system.clusterId();
        String   HAMILTONIAN_ID = system.hamiltonianId();

        System.out.println("================================================================================");
        System.out.println("                    CE THERMODYNAMICS WORKBENCH");
        System.out.println("================================================================================");
        System.out.println("Mode      : " + mode);
        System.out.println("System    : " + elements + "  /  " + structure + "  /  " + model);
        System.out.println("Cluster ID: " + CLUSTER_ID + "   Hamiltonian ID: " + HAMILTONIAN_ID);

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
                        .disorderedClusterFile("clus/BCC_A2-T.txt")   // A2: disordered (HSP)
                        .orderedClusterFile("clus/BCC_B2-T.txt")       // B2: ordered phase
                        .disorderedSymmetryGroup("BCC_A2-SG")
                        .orderedSymmetryGroup("BCC_B2-SG")
                        .numComponents(2)
                        .build();

                AllClusterData clusterData = ClusterIdentificationWorkflow.identify(config);
                System.out.println("Identification complete: " + clusterData.getSummary());

                List<Cluster> maxClusters = InputLoader.parseClusterFile(inputsDir, "clus/BCC_A2-T.txt");
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
                    cecWorkflow.scaffoldFromClusterData(HAMILTONIAN_ID, elements, structure, model);
                    System.out.println("Saved: ~/CEWorkbench/hamiltonians/" + HAMILTONIAN_ID + "/hamiltonian.json");
                    System.out.println("  -> Edit hamiltonian.json to add real ECI values, then run type2.");
                }
            }

            // ------------------------------------------------------------------
            // TYPE-2: Thermodynamic Calculation (CVM temperature scan)
            // ------------------------------------------------------------------
            if (mode.equals("type2") || mode.equals("all")) {

                System.out.println("\n=== TYPE-2: Thermodynamic Calculation (CVM) ===\n");

                CVMEngine cvmEngine = new CVMEngine();
                ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
                        clusterStore, cecWorkflow, cvmEngine, null
                );
                CalculationService service = new CalculationService(thermoWorkflow);

                double[] composition = {0.5, 0.5};
                double tStart = 300.0;
                double tEnd   = 2000.0;
                double tStep  = 100.0;

                System.out.println("System      : " + CLUSTER_ID + " / " + HAMILTONIAN_ID);
                System.out.println("Composition : x_B = " + composition[1]);
                System.out.println("T range     : " + tStart + " K to " + tEnd + " K, step " + tStep + " K\n");

                List<ThermodynamicResult> results = service.runLineScanTemperature(
                        CLUSTER_ID, HAMILTONIAN_ID, composition, tStart, tEnd, tStep, "CVM"
                );

                System.out.println(String.format("  %-8s  %-16s  %-16s  %-16s",
                        "T (K)", "G (J/mol)", "H (J/mol)", "S (J/mol/K)"));
                System.out.println("  " + "-".repeat(62));
                for (ThermodynamicResult r : results) {
                    System.out.println(String.format("  %-8.1f  %-16.4f  %-16.4f  %-16.6f",
                            r.temperature, r.gibbsEnergy, r.enthalpy, r.entropy));
                }
            }

        } catch (Exception e) {
            System.err.println("\nError: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("\n================================================================================");
    }
}
