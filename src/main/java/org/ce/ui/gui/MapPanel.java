package org.ce.ui.gui;

import org.ce.calculation.QuantityDescriptor;
import org.ce.model.ModelSession;
import org.ce.model.result.ThermodynamicResult;
import org.ce.calculation.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;
import java.util.List;
import java.util.function.Consumer;

/**
 * Parameter panel for 2D grid scans (T × x composition).
 *
 * <p>Session identity and engine selection are managed by {@link SessionBar}.
 * Results are forwarded to {@link OutputPanel#showMapResult}.</p>
 */
public class MapPanel extends JPanel {

    // ── VS Code dark colours ──────────────────────────────────────────────────
    private static final Color BG       = new Color(0x252526);
    private static final Color LABEL_FG = new Color(0xCCCCCC);
    private static final Color DIM_FG   = new Color(0x858585);

    private final CalculationService        service;
    private final WorkbenchContext          context;
    private final Consumer<String>          statusSink;
    private final Consumer<String>          logSink;
    private final OutputPanel               outputPanel;

    // Temperature range
    private final JSpinner tStartField = new JSpinner(new SpinnerNumberModel(300.0,  1.0, 10000.0, 100.0));
    private final JSpinner tEndField   = new JSpinner(new SpinnerNumberModel(1500.0, 1.0, 10000.0, 100.0));
    private final JSpinner tStepField  = new JSpinner(new SpinnerNumberModel(100.0,  1.0,  1000.0,  50.0));

    // Composition range
    private final JSpinner xStartField = new JSpinner(new SpinnerNumberModel(0.05, 0.0, 1.0, 0.05));
    private final JSpinner xEndField   = new JSpinner(new SpinnerNumberModel(0.95, 0.0, 1.0, 0.05));
    private final JSpinner xStepField  = new JSpinner(new SpinnerNumberModel(0.05, 0.001, 0.5, 0.01));

    // MCS parameters (shown only when engine is MCS)
    private final JSpinner lField      = new JSpinner(new SpinnerNumberModel(4, 1, 32, 1));
    private final JSpinner nEquilField = new JSpinner(new SpinnerNumberModel(1000, 100, 100000, 100));
    private final JSpinner nAvgField   = new JSpinner(new SpinnerNumberModel(2000, 100, 100000, 100));
    private JPanel mcsParamsPanel;

    // Buttons
    private final JButton runButton   = new JButton("Run Map");
    private final JButton abortButton = new JButton("Abort");

    private SwingWorker<?, ?> activeWorker = null;

    public MapPanel(org.ce.CEWorkbenchContext appCtx,
                    WorkbenchContext context,
                    Consumer<String> statusSink,
                    Consumer<String> logSink,
                    OutputPanel outputPanel,
                    QuantityDescriptor.SelectionModel quantityModel) {
        this.service     = appCtx.getCalculationService();
        this.context     = context;
        this.statusSink  = statusSink;
        this.logSink     = logSink;
        this.outputPanel = outputPanel;

        setBackground(BG);
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(8, 12, 8, 12));

        JPanel content = new JPanel();
        content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
        content.setBackground(BG);

        // QuantitiesPanel at top
        QuantitiesPanel qp = new QuantitiesPanel(quantityModel, context);
        qp.setAlignmentX(LEFT_ALIGNMENT);
        qp.setBorder(BorderFactory.createEmptyBorder(0, 0, 8, 0));
        content.add(qp);

        buildForm(content);
        add(content, BorderLayout.NORTH);

        context.addSessionListener(session -> SwingUtilities.invokeLater(() -> updateFromSession(session)));
        updateFromSession(context.getActiveSession());
    }

    // =========================================================================
    // Form construction
    // =========================================================================

    private void buildForm(JPanel container) {
        // ── Temperature range ─────────────────────────────────────────────────
        container.add(makeTitle("TEMPERATURE RANGE (K)"));
        container.add(makeSpinnerRow("T start:", tStartField));
        container.add(makeSpinnerRow("T end:",   tEndField));
        container.add(makeSpinnerRow("T step:",  tStepField));

        container.add(Box.createVerticalStrut(6));

        // ── Composition range ─────────────────────────────────────────────────
        container.add(makeTitle("COMPOSITION RANGE  (x₂)"));
        container.add(makeSpinnerRow("x start:", xStartField));
        container.add(makeSpinnerRow("x end:",   xEndField));
        container.add(makeSpinnerRow("x step:",  xStepField));

        container.add(Box.createVerticalStrut(8));

        // ── MCS parameters (hidden when CVM) ──────────────────────────────────
        mcsParamsPanel = new JPanel();
        mcsParamsPanel.setLayout(new BoxLayout(mcsParamsPanel, BoxLayout.Y_AXIS));
        mcsParamsPanel.setBackground(BG);
        mcsParamsPanel.setAlignmentX(LEFT_ALIGNMENT);

        mcsParamsPanel.add(makeTitle("MCS PARAMETERS"));
        mcsParamsPanel.add(makeSpinnerRow("Lattice L:",    lField));
        mcsParamsPanel.add(makeSpinnerRow("Equil. sweeps:", nEquilField));
        mcsParamsPanel.add(makeSpinnerRow("Avg. sweeps:",  nAvgField));
        mcsParamsPanel.setVisible(false);
        container.add(mcsParamsPanel);

        container.add(Box.createVerticalStrut(10));

        // ── Buttons ───────────────────────────────────────────────────────────
        runButton.setEnabled(false);
        abortButton.setEnabled(false);
        styleButton(runButton,   new Color(0x0E639C));
        styleButton(abortButton, new Color(0x5A1D1D));

        runButton.addActionListener(e -> startScan());
        abortButton.addActionListener(e -> { if (activeWorker != null) activeWorker.cancel(true); });

        JPanel btnRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        btnRow.setBackground(BG);
        btnRow.setAlignmentX(LEFT_ALIGNMENT);
        btnRow.add(runButton);
        btnRow.add(abortButton);
        container.add(btnRow);
    }

    // =========================================================================
    // Session listener
    // =========================================================================

    private void updateFromSession(ModelSession session) {
        boolean hasSession = (session != null);
        boolean isMcs = hasSession && session.engineConfig.isMcs();
        runButton.setEnabled(hasSession);
        if (mcsParamsPanel != null) {
            mcsParamsPanel.setVisible(isMcs);
            revalidate();
            repaint();
        }
    }

    // =========================================================================
    // Scan execution
    // =========================================================================

    private void startScan() {
        ModelSession session = context.getActiveSession();
        if (session == null) {
            statusSink.accept("No active session — rebuild session using the Session Bar first.");
            return;
        }

        double tStart = ((Number) tStartField.getValue()).doubleValue();
        double tEnd   = ((Number) tEndField.getValue()).doubleValue();
        double tStep  = ((Number) tStepField.getValue()).doubleValue();
        double xStart = ((Number) xStartField.getValue()).doubleValue();
        double xEnd   = ((Number) xEndField.getValue()).doubleValue();
        double xStep  = ((Number) xStepField.getValue()).doubleValue();

        int nT = (int) Math.round((tEnd - tStart) / tStep) + 1;
        int nX = (int) Math.round((xEnd - xStart) / xStep) + 1;

        runButton.setEnabled(false);
        abortButton.setEnabled(true);
        statusSink.accept(String.format("Map scan: %d × %d = %d points…", nT, nX, nT * nX));
        logSink.accept(String.format("━━━ Map Scan (T × x) ━━━"));
        logSink.accept(String.format("  T: %.0f – %.0f K, step %.0f K  (%d points)", tStart, tEnd, tStep, nT));
        logSink.accept(String.format("  x: %.3f – %.3f, step %.3f  (%d points)", xStart, xEnd, xStep, nX));
        logSink.accept(String.format("  Total: %d × %d = %d calculations", nT, nX, nT * nX));
        logSink.accept("  Running...");

        activeWorker = new SwingWorker<List<List<ThermodynamicResult>>, String>() {
            @Override
            protected List<List<ThermodynamicResult>> doInBackground() throws Exception {
                return service.runGridScan(session, tStart, tEnd, tStep, xStart, xEnd, xStep);
            }
            @Override
            protected void done() {
                try {
                    List<List<ThermodynamicResult>> grid = get();
                    int total = grid.stream().mapToInt(List::size).sum();
                    logSink.accept("  Done — " + total + " points computed.");
                    outputPanel.showMapResult(grid);
                    statusSink.accept(String.format("Map scan complete: %d points.", total));
                } catch (java.util.concurrent.CancellationException ignored) {
                    logSink.accept("  Map scan cancelled.");
                    statusSink.accept("Map scan cancelled.");
                } catch (Exception ex) {
                    logSink.accept("  Map scan error: " + ex.getMessage());
                    statusSink.accept("Map scan error: " + ex.getMessage());
                } finally {
                    runButton.setEnabled(context.getActiveSession() != null);
                    abortButton.setEnabled(false);
                    activeWorker = null;
                }
            }
        };
        activeWorker.execute();
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private JLabel makeTitle(String text) {
        JLabel lbl = new JLabel(text);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 9));
        lbl.setForeground(DIM_FG);
        lbl.setAlignmentX(LEFT_ALIGNMENT);
        lbl.setBorder(BorderFactory.createEmptyBorder(6, 0, 2, 0));
        return lbl;
    }

    private JPanel makeSpinnerRow(String label, JSpinner spinner) {
        JPanel row = new JPanel(new BorderLayout(8, 0));
        row.setBackground(BG);
        row.setAlignmentX(LEFT_ALIGNMENT);
        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, 28));

        JLabel lbl = new JLabel(label);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        lbl.setForeground(LABEL_FG);
        lbl.setPreferredSize(new Dimension(100, 24));

        row.add(lbl, BorderLayout.WEST);
        row.add(spinner, BorderLayout.CENTER);
        row.setBorder(BorderFactory.createEmptyBorder(1, 0, 1, 0));
        return row;
    }

    private void styleButton(JButton btn, Color bg) {
        btn.setBackground(bg);
        btn.setForeground(Color.WHITE);
        btn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 11));
        btn.setFocusPainted(false);
        btn.setBorderPainted(false);
        btn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
    }
}
