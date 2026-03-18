package org.ce.ui.gui;

import org.ce.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;

/**
 * Main application window.
 *
 * <p>A tabbed JFrame containing four panels that mirror the two fundamental
 * task types of the CE Workbench:
 * <ol>
 *   <li>Data Preparation — run cluster identification (Type-1)</li>
 *   <li>CEC Management  — scaffold / inspect CEC database</li>
 *   <li>Calculate       — single-point thermodynamic calculation (Type-2)</li>
 *   <li>Results         — display the last EquilibriumState</li>
 * </ol>
 * </p>
 */
public class MainWindow extends JFrame {

    private final ResultsPanel resultsPanel;

    public MainWindow(CalculationService calculationService,
                      org.ce.storage.ClusterDataStore clusterStore,
                      org.ce.workflow.cec.CECManagementWorkflow cecWorkflow) {

        super("CE Thermodynamics Workbench");

        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(900, 680);
        setMinimumSize(new Dimension(700, 500));
        setLocationRelativeTo(null);

        resultsPanel = new ResultsPanel();

        JTabbedPane tabs = new JTabbedPane();
        tabs.addTab("Data Prep",       new DataPreparationPanel(clusterStore));
        tabs.addTab("CEC Management",  new CECManagementPanel(cecWorkflow));
        tabs.addTab("Calculate",       new CalculationPanel(calculationService, resultsPanel, tabs));
        tabs.addTab("Results",         resultsPanel);

        add(tabs, BorderLayout.CENTER);

        // Status bar
        JLabel statusBar = new JLabel(" Ready");
        statusBar.setBorder(BorderFactory.createEtchedBorder());
        add(statusBar, BorderLayout.SOUTH);
    }
}
