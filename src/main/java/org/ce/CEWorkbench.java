package org.ce;

import org.ce.domain.engine.cvm.CVMEngine;
import org.ce.domain.engine.mcs.MCSEngine;
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
        CEWorkbenchContext appCtx = new CEWorkbenchContext();

        // Launch GUI on EDT
        SwingUtilities.invokeLater(() -> {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception ignored) {}

            MainWindow window = new MainWindow(appCtx);
            window.setVisible(true);
        });
    }
}
