package org.ce;

import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.engine.MCSEngine;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.Workspace;
import org.ce.storage.HamiltonianStore;
import org.ce.ui.gui.MainWindow;
import org.ce.workflow.CalculationService;
import org.ce.workflow.cec.CECManagementWorkflow;
import org.ce.workflow.thermo.ThermodynamicWorkflow;

import javax.swing.*;
import java.nio.file.Path;

/**
 * GUI entry point for the CE Thermodynamics Workbench.
 *
 * <p>Builds the application infrastructure (workspace, stores, workflows, service)
 * and launches the Swing main window on the Event Dispatch Thread.</p>
 *
 * <p>For the CLI pipeline demo see {@link org.ce.ui.cli.Main}.</p>
 */
public class CEWorkbench {

    public static void main(String[] args) {

        // Infrastructure
        Workspace workspace = new Workspace();
        ClusterDataStore  clusterStore      = new ClusterDataStore(workspace);
        HamiltonianStore  hamiltonianStore  = new HamiltonianStore(workspace);
        CECManagementWorkflow cecWorkflow   = new CECManagementWorkflow(hamiltonianStore, clusterStore);
        CVMEngine         cvmEngine         = new CVMEngine();
        MCSEngine         mcsEngine         = new MCSEngine();
        ThermodynamicWorkflow thermoWorkflow = new ThermodynamicWorkflow(
                clusterStore, cecWorkflow, cvmEngine, mcsEngine
        );
        CalculationService service = new CalculationService(thermoWorkflow);

        // Launch GUI on EDT
        SwingUtilities.invokeLater(() -> {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception ignored) {}

            MainWindow window = new MainWindow(service, clusterStore, cecWorkflow);
            window.setVisible(true);
        });
    }
}
