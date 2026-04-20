package org.ce.ui.gui;

import org.ce.calculation.QuantityDescriptor;
import org.ce.model.ProgressEvent;
import org.ce.model.hamiltonian.CECEntry;

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

    public MainWindow(org.ce.CEWorkbenchContext appCtx) {

        super("CE Thermodynamics Workbench");
        applyDarkTheme();

        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        // Size window responsively based on available screen space
        Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        int width = Math.min(1200, (int)(screenSize.width * 0.9));
        int height = Math.min(760, (int)(screenSize.height * 0.9));
        setSize(width, height);
        setMinimumSize(new Dimension(900, 580));
        setLocationRelativeTo(null);

        // ── shared session state ──────────────────────────────────────────────
        WorkbenchContext context = new WorkbenchContext();

        // ── shared quantity selection model ───────────────────────────────────
        QuantityDescriptor.SelectionModel quantityModel = new QuantityDescriptor.SelectionModel();

        // ── chrome ───────────────────────────────────────────────────────────
        statusBar = new StatusBar();
        HeaderBar header = new HeaderBar(context);

        // ── output panel (always visible on the right) ────────────────────────
        OutputPanel outputPanel = new OutputPanel(context, quantityModel);

        // ── sinks wired to the output panel ───────────────────────────────────
        // Unified log sink: GUI output box + CLI terminal
        java.util.function.Consumer<String> logSink = msg -> {
            outputPanel.appendLog(msg);
            System.out.println(msg);
        };
        appCtx.setLogSink(logSink);

        java.util.function.Consumer<String> statusSink = this::postStatus;
        java.util.function.BiConsumer<CECEntry, CECEntry> cecResultSink = outputPanel::showCECResult;
        java.util.function.Function<CECEntry, Boolean> cecEditApplySink = outputPanel::applyCECEdits;
        java.util.function.Consumer<Object> resultSink = result -> {
            switch (result) {
                case org.ce.calculation.CalculationResult.Single s -> {
                    outputPanel.showResult(s.value());
                }
                case org.ce.calculation.CalculationResult.Grid g -> {
                    java.util.List<java.util.List<org.ce.model.ThermodynamicResult>> grid = g.values();
                    if (grid.size() == 1 && grid.get(0).size() == 1) {
                        outputPanel.showResult(grid.get(0).get(0));
                    } else if (grid.size() == 1) {
                        java.util.List<org.ce.model.ThermodynamicResult> row = grid.get(0);
                        boolean tVaries = row.size() > 1 &&
                                Math.abs(row.get(0).temperature - row.get(row.size()-1).temperature) > 1e-3;
                        outputPanel.showLineScanResult(row, tVaries ? "T" : "X");
                        logSink.accept("\n" + org.ce.calculation.ResultFormatter.table(row));
                    } else if (!grid.isEmpty() && grid.get(0).size() == 1) {
                        java.util.List<org.ce.model.ThermodynamicResult> col = grid.stream()
                                .map(r -> r.get(0)).toList();
                        boolean tVaries = col.size() > 1 &&
                                Math.abs(col.get(0).temperature - col.get(col.size()-1).temperature) > 1e-3;
                        outputPanel.showLineScanResult(col, tVaries ? "T" : "X");
                        logSink.accept("\n" + org.ce.calculation.ResultFormatter.table(col));
                    } else {
                        outputPanel.showMapResult(grid);
                    }
                }
                default -> { /* unknown result type — ignore */ }
            }
        };
        java.util.function.Consumer<ProgressEvent> chartSink = outputPanel::onChartEvent;

        // ── parameter panels (go into the explorer) ───────────────────────────
        DataPreparationPanel dataPrepPanel = new DataPreparationPanel(
                appCtx, context, statusSink);

        CECManagementPanel cecPanel = new CECManagementPanel(
                appCtx, context, statusSink, cecResultSink, cecEditApplySink);

        DynamicCalculationPanel calcPanel = new DynamicCalculationPanel(
                appCtx, context, statusSink, logSink, resultSink, chartSink, quantityModel);

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
        contentSplit.setDividerLocation(290);  // pixel fallback until first resize
        contentSplit.addComponentListener(new java.awt.event.ComponentAdapter() {
            private boolean initialized = false;
            @Override
            public void componentResized(java.awt.event.ComponentEvent e) {
                if (!initialized && contentSplit.getWidth() > 0) {
                    initialized = true;
                    contentSplit.setDividerLocation(0.27);
                    contentSplit.removeComponentListener(this);
                }
            }
        });
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
        add(header, BorderLayout.NORTH);
        add(centre,    BorderLayout.CENTER);
        add(statusBar, BorderLayout.SOUTH);

        navigate(2);    // start on Thermodynamics panel
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
