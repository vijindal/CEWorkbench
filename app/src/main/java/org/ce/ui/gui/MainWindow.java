package org.ce.ui.gui;

import org.ce.domain.engine.ProgressEvent;
import org.ce.storage.ClusterDataStore;
import org.ce.workflow.CalculationService;
import org.ce.workflow.cec.CECManagementWorkflow;

import javax.swing.*;
import java.awt.*;

/**
 * Main application window — VS Code-style three-column dark layout.
 *
 * <pre>
 *  ┌──────────────────────────────────────────────────────────────┐
 *  │  HeaderBar  (dark title bar, app name + system identity)     │
 *  ├────┬───────────────────┬─────────────────────────────────────┤
 *  │    │                   │                                     │
 *  │ A  │  ExplorerPanel    │  OutputPanel                        │
 *  │ c  │  (parameters for  │  ┌───────────────────────────────┐  │
 *  │ t  │  selected section)│  │  RESULTS  (top)               │  │
 *  │ i  │                   │  ├───────────────────────────────┤  │
 *  │ v  │                   │  │  OUTPUT / log  (bottom)       │  │
 *  │ i  │                   │  └───────────────────────────────┘  │
 *  │ t  │                   │                                     │
 *  │ y  │                   │                                     │
 *  ├────┴───────────────────┴─────────────────────────────────────┤
 *  │  StatusBar  (VS Code blue footer)                            │
 *  └──────────────────────────────────────────────────────────────┘
 * </pre>
 *
 * <p>All panels share a {@link WorkbenchContext} for system-identity propagation.
 * Log output from every panel flows to {@link OutputPanel#appendLog(String)}.
 * Completed results flow to {@link OutputPanel#showResult}.</p>
 */
public class MainWindow extends JFrame {

    private final StatusBar     statusBar;
    private final ActivityBar   activityBar;
    private final ExplorerPanel explorerPanel;

    public MainWindow(CalculationService calculationService,
                      ClusterDataStore clusterStore,
                      CECManagementWorkflow cecWorkflow) {

        super("CE Thermodynamics Workbench");
        applyDarkTheme();

        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1200, 760);
        setMinimumSize(new Dimension(900, 580));
        setLocationRelativeTo(null);

        // ── shared session state ──────────────────────────────────────────────
        WorkbenchContext context = new WorkbenchContext();

        // ── chrome ───────────────────────────────────────────────────────────
        statusBar = new StatusBar();
        HeaderBar header = new HeaderBar(context);

        // ── output panel (always visible on the right) ────────────────────────
        OutputPanel outputPanel = new OutputPanel(context);

        // ── sinks wired to the output panel ───────────────────────────────────
        java.util.function.Consumer<String> logSink    = outputPanel::appendLog;
        java.util.function.Consumer<String> statusSink = this::postStatus;
        java.util.function.Consumer<org.ce.domain.result.ThermodynamicResult> resultSink =
                outputPanel::showResult;
        java.util.function.Consumer<ProgressEvent> chartSink = outputPanel::onChartEvent;

        // ── parameter panels (go into the explorer) ───────────────────────────
        DataPreparationPanel dataPrepPanel = new DataPreparationPanel(
                clusterStore, context, statusSink, logSink);

        CECManagementPanel cecPanel = new CECManagementPanel(
                cecWorkflow, context, statusSink, logSink);

        CalculationPanel calcPanel = new CalculationPanel(
                calculationService, context, statusSink, logSink, resultSink, chartSink);

        // ── explorer panel ────────────────────────────────────────────────────
        explorerPanel = new ExplorerPanel();
        explorerPanel.addCard(dataPrepPanel, 0);
        explorerPanel.addCard(cecPanel,      1);
        explorerPanel.addCard(calcPanel,     2);

        // ── activity bar ──────────────────────────────────────────────────────
        Runnable[] navCallbacks = {
            () -> navigate(0),
            () -> navigate(1),
            () -> navigate(2),
        };
        activityBar = new ActivityBar(navCallbacks);

        // ── centre zone: ActivityBar | splittable (Explorer + Output) ─────────
        JSplitPane contentSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                explorerPanel, outputPanel);
        contentSplit.setDividerLocation(290);
        contentSplit.setDividerSize(1);
        contentSplit.setContinuousLayout(true);
        contentSplit.setBorder(null);
        contentSplit.setBackground(new Color(0x1A1A1A));

        JPanel centre = new JPanel(new BorderLayout());
        centre.setBackground(new Color(0x1E1E1E));
        centre.add(activityBar,   BorderLayout.WEST);
        centre.add(contentSplit,  BorderLayout.CENTER);

        // ── assemble ──────────────────────────────────────────────────────────
        setLayout(new BorderLayout());
        add(header,   BorderLayout.NORTH);
        add(centre,   BorderLayout.CENTER);
        add(statusBar, BorderLayout.SOUTH);

        navigate(2);    // start on Calculation panel
    }

    // =========================================================================
    // Navigation
    // =========================================================================

    private void navigate(int index) {
        activityBar.setActive(index);
        explorerPanel.showCard(index);
    }

    // =========================================================================
    // Status
    // =========================================================================

    /** Posts a status message (safe to call from any thread). */
    public void postStatus(String message) {
        if (SwingUtilities.isEventDispatchThread()) {
            statusBar.post(message);
        } else {
            SwingUtilities.invokeLater(() -> statusBar.post(message));
        }
    }

    // =========================================================================
    // Look & Feel — VS Code Dark+
    // =========================================================================

    private static void applyDarkTheme() {
        try {
            for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (Exception ignored) {}

        // ── global dark overrides ─────────────────────────────────────────────
        // Panel / button background
        UIManager.put("control",                new Color(0x3C3C3C));
        // Text field / text area background
        UIManager.put("nimbusLightBackground",  new Color(0x3C3C3C));
        // Base colour used to derive most Nimbus chrome
        UIManager.put("nimbusBase",             new Color(0x1A1A2E));
        UIManager.put("nimbusBlueGrey",         new Color(0x444444));
        // Default foreground text
        UIManager.put("text",                   new Color(0xD4D4D4));
        // Selection
        UIManager.put("nimbusSelectedText",     Color.WHITE);
        UIManager.put("nimbusSelectionBackground", new Color(0x094771));
        // Focus ring
        UIManager.put("nimbusFocus",            new Color(0x007ACC));
        // Tooltip background
        UIManager.put("info",                   new Color(0x2D2D2D));
        UIManager.put("nimbusInfoBackground",   new Color(0x2D2D2D));
        // Titled border text
        UIManager.put("TitledBorder.titleColor", new Color(0x9CDCFE));
        // Table
        UIManager.put("Table.background",       new Color(0x252526));
        UIManager.put("Table.alternateRowColor", new Color(0x2A2D2E));
        UIManager.put("Table.gridColor",        new Color(0x3C3C3C));
        UIManager.put("TableHeader.background", new Color(0x2D2D2D));
        UIManager.put("TableHeader.foreground", new Color(0xBBBBBB));
        // Scroll bar
        UIManager.put("ScrollBar.thumb",        new Color(0x424242));
        UIManager.put("ScrollBar.track",        new Color(0x1E1E1E));
        // Split pane
        UIManager.put("SplitPane.background",           new Color(0x1A1A1A));
        UIManager.put("SplitPaneDivider.draggingColor", new Color(0x007ACC));
    }
}
