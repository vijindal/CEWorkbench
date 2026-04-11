package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.ResultFormatter;
import org.ce.model.ModelSession;
import org.ce.model.ThermodynamicResult;
import org.ce.calculation.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;
import java.util.List;
import java.util.function.Consumer;

/**
 * Parameter panel for 1-D line scans — either a temperature sweep at fixed
 * composition, or a composition sweep at fixed temperature.
 *
 * <p>Session identity and engine selection are managed by {@link SessionBar}.
 * This panel reads the active {@link ModelSession} from {@link WorkbenchContext}
 * and enables its buttons only when a session is ready.</p>
 *
 * <p>Results are forwarded to {@link OutputPanel#showLineScanResult}.</p>
 */
public class LineScanPanel extends JPanel {

    // ── VS Code dark colours ──────────────────────────────────────────────────
    private static final Color BG        = new Color(0x252526);
    private static final Color LABEL_FG  = new Color(0xCCCCCC);
    private static final Color DIM_FG    = new Color(0x858585);

    private final org.ce.CEWorkbenchContext appCtx;
    private final CalculationService service;
    private final WorkbenchContext context;
    private final Consumer<String> statusSink;
    private final Consumer<String> logSink;
    private final OutputPanel outputPanel;

    // Scan type
    private final JRadioButton tempScanBtn = new JRadioButton("Temperature scan");
    private final JRadioButton compScanBtn = new JRadioButton("Composition scan");

    // Temperature scan parameters
    private final JSpinner tStartField = new JSpinner(new SpinnerNumberModel(500.0, 1.0, 10000.0, 100.0));
    private final JSpinner tEndField   = new JSpinner(new SpinnerNumberModel(1500.0, 1.0, 10000.0, 100.0));
    private final JSpinner tStepField  = new JSpinner(new SpinnerNumberModel(50.0, 1.0, 1000.0, 10.0));
    private final JTextField fixedCompField = new JTextField("0.5, 0.5", 12);

    // Composition scan parameters
    private final JSpinner xStartField = new JSpinner(new SpinnerNumberModel(0.05, 0.0, 1.0, 0.05));
    private final JSpinner xEndField   = new JSpinner(new SpinnerNumberModel(0.95, 0.0, 1.0, 0.05));
    private final JSpinner xStepField  = new JSpinner(new SpinnerNumberModel(0.05, 0.001, 0.5, 0.01));
    private final JTextField fixedTempField = new JTextField("1000.0", 10);

    // MCS parameters (shown only when engine is MCS)
    private final JSpinner lField      = new JSpinner(new SpinnerNumberModel(4, 1, 32, 1));
    private final JSpinner nEquilField = new JSpinner(new SpinnerNumberModel(1000, 100, 100000, 100));
    private final JSpinner nAvgField   = new JSpinner(new SpinnerNumberModel(2000, 100, 100000, 100));
    private JPanel mcsParamsPanel;

    // Parameter group panels (for show/hide)
    private JPanel tempParamsPanel;
    private JPanel compParamsPanel;

    // Buttons
    private final JButton scanButton  = new JButton("Scan");
    private final JButton abortButton = new JButton("Abort");

    private SwingWorker<?, ?> activeWorker = null;

    public LineScanPanel(org.ce.CEWorkbenchContext appCtx,
                         WorkbenchContext context,
                         Consumer<String> statusSink,
                         Consumer<String> logSink,
                         OutputPanel outputPanel,
                         QuantityDescriptor.SelectionModel quantityModel) {
        this.appCtx = appCtx;
        this.service = appCtx.getCalculationService();
        this.context = context;
        this.statusSink = statusSink;
        this.logSink = logSink;
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

        hookListeners();
        updateFromSession(context.getActiveSession());
        updateParamVisibility();
    }

    // =========================================================================
    // Form construction
    // =========================================================================

    private void buildForm(JPanel container) {
        // ── Scan type ─────────────────────────────────────────────────────────
        container.add(makeTitle("SCAN TYPE"));
        ButtonGroup group = new ButtonGroup();
        group.add(tempScanBtn);
        group.add(compScanBtn);
        tempScanBtn.setSelected(true);
        styleRadio(tempScanBtn);
        styleRadio(compScanBtn);

        JPanel radioRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        radioRow.setBackground(BG);
        radioRow.add(tempScanBtn);
        radioRow.add(Box.createHorizontalStrut(12));
        radioRow.add(compScanBtn);
        radioRow.setAlignmentX(LEFT_ALIGNMENT);
        container.add(radioRow);
        container.add(Box.createVerticalStrut(8));

        // ── Temperature scan params ───────────────────────────────────────────
        tempParamsPanel = new JPanel();
        tempParamsPanel.setLayout(new BoxLayout(tempParamsPanel, BoxLayout.Y_AXIS));
        tempParamsPanel.setBackground(BG);
        tempParamsPanel.setAlignmentX(LEFT_ALIGNMENT);

        tempParamsPanel.add(makeTitle("TEMPERATURE RANGE (K)"));
        tempParamsPanel.add(makeSpinnerRow("T start:", tStartField));
        tempParamsPanel.add(makeSpinnerRow("T end:", tEndField));
        tempParamsPanel.add(makeSpinnerRow("T step:", tStepField));
        tempParamsPanel.add(Box.createVerticalStrut(6));
        tempParamsPanel.add(makeTitle("FIXED COMPOSITION"));
        tempParamsPanel.add(makeFieldRow("x (x₁,x₂):", fixedCompField));
        container.add(tempParamsPanel);

        // ── Composition scan params ───────────────────────────────────────────
        compParamsPanel = new JPanel();
        compParamsPanel.setLayout(new BoxLayout(compParamsPanel, BoxLayout.Y_AXIS));
        compParamsPanel.setBackground(BG);
        compParamsPanel.setAlignmentX(LEFT_ALIGNMENT);

        compParamsPanel.add(makeTitle("COMPOSITION RANGE  (x₂)"));
        compParamsPanel.add(makeSpinnerRow("x start:", xStartField));
        compParamsPanel.add(makeSpinnerRow("x end:", xEndField));
        compParamsPanel.add(makeSpinnerRow("x step:", xStepField));
        compParamsPanel.add(Box.createVerticalStrut(6));
        compParamsPanel.add(makeTitle("FIXED TEMPERATURE (K)"));
        compParamsPanel.add(makeFieldRow("T (K):", fixedTempField));
        container.add(compParamsPanel);

        container.add(Box.createVerticalStrut(8));

        // ── MCS parameters (hidden when CVM) ──────────────────────────────────
        mcsParamsPanel = new JPanel();
        mcsParamsPanel.setLayout(new BoxLayout(mcsParamsPanel, BoxLayout.Y_AXIS));
        mcsParamsPanel.setBackground(BG);
        mcsParamsPanel.setAlignmentX(LEFT_ALIGNMENT);

        mcsParamsPanel.add(makeTitle("MCS PARAMETERS"));
        mcsParamsPanel.add(makeSpinnerRow("Lattice L:", lField));
        mcsParamsPanel.add(makeSpinnerRow("Equil. sweeps:", nEquilField));
        mcsParamsPanel.add(makeSpinnerRow("Avg. sweeps:", nAvgField));
        mcsParamsPanel.setVisible(false);
        container.add(mcsParamsPanel);

        container.add(Box.createVerticalStrut(10));

        // ── Buttons ───────────────────────────────────────────────────────────
        scanButton.setEnabled(false);
        abortButton.setEnabled(false);

        styleButton(scanButton,  new Color(0x0E639C));
        styleButton(abortButton, new Color(0x5A1D1D));

        JPanel btnRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        btnRow.setBackground(BG);
        btnRow.setAlignmentX(LEFT_ALIGNMENT);
        btnRow.add(scanButton);
        btnRow.add(abortButton);
        container.add(btnRow);
    }

    // =========================================================================
    // Listeners
    // =========================================================================

    private void hookListeners() {
        tempScanBtn.addActionListener(e -> updateParamVisibility());
        compScanBtn.addActionListener(e -> updateParamVisibility());

        context.addSessionListener(session -> SwingUtilities.invokeLater(() -> updateFromSession(session)));

        scanButton.addActionListener(e -> startScan());
        abortButton.addActionListener(e -> cancelScan());
    }

    private void updateParamVisibility() {
        boolean isTemp = tempScanBtn.isSelected();
        tempParamsPanel.setVisible(isTemp);
        compParamsPanel.setVisible(!isTemp);
        revalidate();
        repaint();
    }

    private void updateFromSession(ModelSession session) {
        boolean hasSession = (session != null);
        boolean isMcs = hasSession && session.engineConfig.isMcs();
        scanButton.setEnabled(hasSession);
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

        scanButton.setEnabled(false);
        abortButton.setEnabled(true);

        boolean isTempScan = tempScanBtn.isSelected();
        final int     mcsL      = ((Number) lField.getValue()).intValue();
        final int     mcsNEquil = ((Number) nEquilField.getValue()).intValue();
        final int     mcsNAvg   = ((Number) nAvgField.getValue()).intValue();

        ModelSpecifications modelSpecs = new ModelSpecifications(
                session.systemId.elements, session.systemId.structure, session.systemId.model, session.engineConfig);
        CalculationSpecifications calcSpecs = new CalculationSpecifications(Property.GIBBS_ENERGY, Mode.LINE_SCAN);
        calcSpecs.set(Parameter.MCS_L, mcsL);
        calcSpecs.set(Parameter.MCS_NEQUIL, mcsNEquil);
        calcSpecs.set(Parameter.MCS_NAVG, mcsNAvg);

        String scanLabel;
        if (isTempScan) {
            double tStart = ((Number) tStartField.getValue()).doubleValue();
            double tEnd   = ((Number) tEndField.getValue()).doubleValue();
            double tStep  = ((Number) tStepField.getValue()).doubleValue();
            double[] comp = parseComposition(fixedCompField.getText());
            if (comp == null) {
                statusSink.accept("Invalid composition format. Use comma-separated values, e.g. 0.5, 0.5");
                scanButton.setEnabled(true);
                abortButton.setEnabled(false);
                return;
            }
            calcSpecs.set(Parameter.T_START, tStart);
            calcSpecs.set(Parameter.T_END, tEnd);
            calcSpecs.set(Parameter.T_STEP, tStep);
            calcSpecs.set(Parameter.COMPOSITION, comp);
            scanLabel = String.format("Temperature scan: %.0f–%.0f K, step %.0f K, x=%s …",
                    tStart, tEnd, tStep, fixedCompField.getText().trim());
        } else {
            double xStart = ((Number) xStartField.getValue()).doubleValue();
            double xEnd   = ((Number) xEndField.getValue()).doubleValue();
            double xStep  = ((Number) xStepField.getValue()).doubleValue();
            double temp;
            try {
                temp = Double.parseDouble(fixedTempField.getText().trim());
            } catch (NumberFormatException ex) {
                statusSink.accept("Invalid temperature value.");
                scanButton.setEnabled(true);
                abortButton.setEnabled(false);
                return;
            }
            calcSpecs.set(Parameter.X_START, xStart);
            calcSpecs.set(Parameter.X_END, xEnd);
            calcSpecs.set(Parameter.X_STEP, xStep);
            calcSpecs.set(Parameter.TEMPERATURE, temp);
            scanLabel = String.format("Composition scan: x=%.3f–%.3f, step %.3f, T=%.0f K …",
                    xStart, xEnd, xStep, temp);
        }

        statusSink.accept(scanLabel);

        activeWorker = new SwingWorker<List<ThermodynamicResult>, String>() {
            @Override
            protected List<ThermodynamicResult> doInBackground() throws Exception {
                logSink.accept("━━━ " + (isTempScan ? "Temperature" : "Composition") + " Scan ━━━");
                logSink.accept("  Using Discovery-based execution dispatcher");
                logSink.accept("  Running...");

                // Calculation Layer Role: Unified Execution
                List<ThermodynamicResult> results = service.executeScan(modelSpecs, calcSpecs, this::publish);

                logSink.accept(ResultFormatter.tableHeader());
                logSink.accept(ResultFormatter.tableSeparator());
                for (ThermodynamicResult r : results) {
                    logSink.accept(ResultFormatter.tableRow(r));
                }
                logSink.accept("  Done — " + results.size() + " points computed.");
                return results;
            }
            @Override
            protected void done() {
                try {
                    List<ThermodynamicResult> results = get();
                    outputPanel.showLineScanResult(results, isTempScan ? "T" : "X");
                    statusSink.accept((isTempScan ? "Temperature" : "Composition") + " scan complete: " + results.size() + " points.");
                } catch (java.util.concurrent.CancellationException ignored) {
                    logSink.accept("  Scan cancelled.");
                    statusSink.accept("Scan cancelled.");
                } catch (Exception ex) {
                    String msg = ex.getCause() != null ? ex.getCause().getMessage() : ex.getMessage();
                    logSink.accept("  Scan error: " + msg);
                    statusSink.accept("Scan error: " + msg);
                } finally {
                    scanButton.setEnabled(context.getActiveSession() != null);
                    abortButton.setEnabled(false);
                    activeWorker = null;
                }
            }
            @Override
            protected void process(List<String> chunks) {
                chunks.forEach(logSink);
            }
        };

        activeWorker.execute();
    }

    private void cancelScan() {
        if (activeWorker != null) activeWorker.cancel(true);
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private double[] parseComposition(String text) {
        try {
            String[] parts = text.split(",");
            double[] comp = new double[parts.length];
            for (int i = 0; i < parts.length; i++) comp[i] = Double.parseDouble(parts[i].trim());
            return comp;
        } catch (NumberFormatException e) {
            return null;
        }
    }

    private JLabel makeTitle(String text) {
        JLabel lbl = new JLabel(text);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 9));
        lbl.setForeground(DIM_FG);
        lbl.setAlignmentX(LEFT_ALIGNMENT);
        lbl.setBorder(BorderFactory.createEmptyBorder(6, 0, 2, 0));
        return lbl;
    }

    private JPanel makeSpinnerRow(String label, JSpinner spinner) {
        return makeRowPanel(label, spinner);
    }

    private JPanel makeFieldRow(String label, JTextField field) {
        field.setBackground(new Color(0x3C3C3C));
        field.setForeground(LABEL_FG);
        field.setCaretColor(LABEL_FG);
        field.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createLineBorder(new Color(0x3C3C3C)),
                BorderFactory.createEmptyBorder(2, 4, 2, 4)));
        return makeRowPanel(label, field);
    }

    private JPanel makeRowPanel(String labelText, JComponent field) {
        JPanel row = new JPanel(new BorderLayout(8, 0));
        row.setBackground(BG);
        row.setAlignmentX(LEFT_ALIGNMENT);
        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, 28));

        JLabel lbl = new JLabel(labelText);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        lbl.setForeground(LABEL_FG);
        lbl.setPreferredSize(new Dimension(90, 24));

        row.add(lbl, BorderLayout.WEST);
        row.add(field, BorderLayout.CENTER);
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

    private void styleRadio(JRadioButton btn) {
        btn.setBackground(BG);
        btn.setForeground(LABEL_FG);
        btn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        btn.setFocusPainted(false);
    }
}
