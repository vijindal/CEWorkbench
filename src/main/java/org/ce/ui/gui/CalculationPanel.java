package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.calculation.QuantityDescriptor;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;
import org.ce.calculation.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;
import java.util.List;
import java.util.function.Consumer;

/**
 * Parameter panel for single-point thermodynamic calculations (shown in the Explorer column).
 *
 * <p>Session identity and engine selection are managed by {@link SessionBar}.
 * This panel reads the active {@link ModelSession} from {@link WorkbenchContext}
 * and enables its buttons only when a session is ready.</p>
 */
public class CalculationPanel extends JPanel {

    // ── VS Code dark colours ──────────────────────────────────────────────────
    private static final Color BG       = new Color(0x252526);
    private static final Color LABEL_FG = new Color(0xCCCCCC);

    private final org.ce.CEWorkbenchContext      appCtx;
    private final CalculationService             service;
    private final WorkbenchContext               context;
    private final Consumer<String>               statusSink;
    private final Consumer<ThermodynamicResult>  resultSink;
    private final Consumer<ProgressEvent>        chartSink;

    // Calculation parameters
    private final JSpinner   temperatureField = new JSpinner(new SpinnerNumberModel(1000.0, 1.0, 10000.0, 100.0));
    private final JTextField compositionField = new JTextField("0.5, 0.5", 12);

    // MCS-specific parameters (shown only when session engine is MCS)
    private final JSpinner lField      = new JSpinner(new SpinnerNumberModel(4, 1, 32, 1));
    private final JSpinner nEquilField = new JSpinner(new SpinnerNumberModel(1000, 100, 100000, 100));
    private final JSpinner nAvgField   = new JSpinner(new SpinnerNumberModel(2000, 100, 100000, 100));
    private JLabel lLabel, nEquilLabel, nAvgLabel;

    private final JButton calcButton  = new JButton("Calculate");
    private final JButton abortButton = new JButton("Abort");
    private final JButton prodButton  = new JButton("Production Run");
    private SwingWorker<?, ?> activeWorker = null;

    public CalculationPanel(org.ce.CEWorkbenchContext appCtx,
                            WorkbenchContext context,
                            Consumer<String> statusSink,
                            Consumer<ThermodynamicResult> resultSink,
                            Consumer<ProgressEvent> chartSink,
                            QuantityDescriptor.SelectionModel quantityModel) {
        this.appCtx     = appCtx;
        this.service    = appCtx.getCalculationService();
        this.context    = context;
        this.statusSink = statusSink;
        this.resultSink = resultSink;
        this.chartSink  = chartSink;

        setBackground(BG);

        // Observe session lifecycle
        context.addSessionListener(this::onSessionChanged);

        abortButton.setEnabled(false);
        abortButton.setForeground(new Color(0xF44747));
        abortButton.addActionListener(e -> { if (activeWorker != null) activeWorker.cancel(true); });

        prodButton.setToolTipText("Run MCS at L=12, 16, 24 and extrapolate to thermodynamic limit (1/N → 0)");
        prodButton.addActionListener(e -> runProductionScan());

        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));

        JPanel content = new JPanel(new BorderLayout(0, 6));
        content.setBackground(BG);
        content.add(new QuantitiesPanel(quantityModel, context), BorderLayout.NORTH);
        content.add(buildForm(), BorderLayout.CENTER);

        add(content, BorderLayout.NORTH);

        updateButtonState();
    }

    // =========================================================================
    // Session observation
    // =========================================================================

    private void onSessionChanged(ModelSession session) {
        updateButtonState();
    }

    private void updateButtonState() {
        boolean hasSession = context.hasActiveSession();
        boolean isMcs = hasSession && context.getActiveSession().engineConfig.isMcs();

        calcButton.setEnabled(hasSession);
        prodButton.setEnabled(hasSession && isMcs);

        if (lLabel != null) {
            lLabel.setVisible(isMcs);
            lField.setVisible(isMcs);
            nEquilLabel.setVisible(isMcs);
            nEquilField.setVisible(isMcs);
            nAvgLabel.setVisible(isMcs);
            nAvgField.setVisible(isMcs);
            revalidate();
            repaint();
        }
    }

    // =========================================================================
    // Form layout
    // =========================================================================

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBackground(BG);
        form.setBorder(BorderFactory.createTitledBorder("Calculation Parameters"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor  = GridBagConstraints.WEST;
        lc.gridx   = 0;
        lc.fill    = GridBagConstraints.HORIZONTAL;
        lc.weightx = 1.0;
        lc.insets  = new Insets(4, 6, 1, 6);

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill    = GridBagConstraints.HORIZONTAL;
        fc.weightx = 1.0;
        fc.gridx   = 0;
        fc.insets  = new Insets(0, 6, 4, 8);

        int row = 0;

        // ── Temperature ───────────────────────────────────────────────────────
        lc.gridy = row++;  form.add(makeLabel("Temperature (K):", LABEL_FG), lc);
        fc.gridy = row++;  form.add(temperatureField, fc);

        // ── Composition ───────────────────────────────────────────────────────
        lc.gridy = row++;  form.add(makeLabel("Composition:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(compositionField, fc);

        // ── MCS parameters (hidden when CVM) ──────────────────────────────────
        lLabel = makeLabel("Lattice size L:", LABEL_FG);
        lc.gridy = row++;  form.add(lLabel, lc);
        fc.gridy = row++;  form.add(lField, fc);

        nEquilLabel = makeLabel("Equil. sweeps:", LABEL_FG);
        lc.gridy = row++;  form.add(nEquilLabel, lc);
        fc.gridy = row++;  form.add(nEquilField, fc);

        nAvgLabel = makeLabel("Avg. sweeps:", LABEL_FG);
        lc.gridy = row++;  form.add(nAvgLabel, lc);
        fc.gridy = row++;  form.add(nAvgField, fc);

        // ── Buttons ───────────────────────────────────────────────────────────
        JPanel btnRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        btnRow.setBackground(BG);
        calcButton.addActionListener(e -> runCalculation());
        btnRow.add(calcButton);
        btnRow.add(abortButton);
        btnRow.add(prodButton);

        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx   = 0;
        bc.gridy   = row;
        bc.fill    = GridBagConstraints.HORIZONTAL;
        bc.weightx = 1.0;
        bc.anchor  = GridBagConstraints.WEST;
        bc.insets  = new Insets(8, 0, 4, 6);
        form.add(btnRow, bc);

        return form;
    }

    private static JLabel makeLabel(String text, Color fg) {
        JLabel l = new JLabel(text);
        l.setForeground(fg);
        return l;
    }

    // =========================================================================
    // Calculation
    // =========================================================================

    private void runProductionScan() {
        ModelSession session = context.getActiveSession();
        if (session == null) { appCtx.log("No active session."); return; }

        double temperature = ((Number) temperatureField.getValue()).doubleValue();
        final double[] composition;
        try {
            composition = parseComposition(compositionField.getText().trim(), session.numComponents());
        } catch (IllegalArgumentException ex) {
            appCtx.log("Invalid composition: " + ex.getMessage());
            statusSink.accept("Error: Invalid composition");
            return;
        }
        int mcsNEquil = ((Number) nEquilField.getValue()).intValue();
        int mcsNAvg   = ((Number) nAvgField.getValue()).intValue();

        ModelSpecifications modelSpecs = new ModelSpecifications(
                session.systemId.elements, session.systemId.structure, session.systemId.model, session.engineConfig);
        CalculationSpecifications calcSpecs = new CalculationSpecifications(Property.ENTHALPY, Mode.FINITE_SIZE_SCALING);
        calcSpecs.set(Parameter.TEMPERATURE, temperature);
        calcSpecs.set(Parameter.COMPOSITION, composition);
        calcSpecs.set(Parameter.MCS_NEQUIL, mcsNEquil);
        calcSpecs.set(Parameter.MCS_NAVG, mcsNAvg);

        calcButton.setEnabled(false);
        prodButton.setEnabled(false);
        abortButton.setEnabled(true);
        appCtx.clearLog();
        appCtx.log("Production Run (FSS): L=12, 16, 24 — MCS with 1/N extrapolation");
        appCtx.log("  Session: " + session.label());
        appCtx.log("  T=" + temperature + " K  composition=" + java.util.Arrays.toString(composition));
        appCtx.log("  nEquil=" + mcsNEquil + "  nAvg=" + mcsNAvg + " per L-value");
        statusSink.accept("Production Run at T=" + temperature + " K…");

        SwingWorker<ThermodynamicResult, Object> worker = new SwingWorker<>() {
            @Override
            protected ThermodynamicResult doInBackground() throws Exception {
                // Calculation Layer Role: Unified Execution
                return service.execute(modelSpecs, calcSpecs, msg -> publish((Object) msg));
            }
            @Override
            protected void process(List<Object> chunks) {
                for (Object o : chunks) {
                    if (o instanceof String)        appCtx.log((String) o);
                    else if (o instanceof ProgressEvent) chartSink.accept((ProgressEvent) o);
                }
            }
            @Override
            protected void done() {
                calcButton.setEnabled(true);
                prodButton.setEnabled(true);
                abortButton.setEnabled(false);
                activeWorker = null;
                if (isCancelled()) { appCtx.log("Aborted."); statusSink.accept("Aborted."); return; }
                try {
                    ThermodynamicResult r = get();
                    appCtx.log("\nProduction run complete.");
                    appCtx.log("  H(∞) = " + String.format("%+.6f ± %.6f", r.enthalpy, r.stdEnthalpy) + " J/mol");
                    statusSink.accept("Done — H(∞) = " + String.format("%+.4f", r.enthalpy) + " J/mol");
                    resultSink.accept(r);
                } catch (Exception ex) {
                    String msg = ex.getCause() != null ? ex.getCause().getMessage() : ex.getMessage();
                    appCtx.log("Error: " + msg);
                    statusSink.accept("Error: " + msg);
                }
            }
        };
        activeWorker = worker;
        worker.execute();
    }

    private void runCalculation() {
        ModelSession session = context.getActiveSession();
        if (session == null) { appCtx.log("No active session."); return; }

        double temperature = ((Number) temperatureField.getValue()).doubleValue();
        double[] composition;
        try {
            composition = parseComposition(compositionField.getText().trim(), session.numComponents());
        } catch (IllegalArgumentException ex) {
            appCtx.log("Invalid composition: " + ex.getMessage());
            statusSink.accept("Error: Invalid composition");
            return;
        }

        int mcsL      = ((Number) lField.getValue()).intValue();
        int mcsNEquil = ((Number) nEquilField.getValue()).intValue();
        int mcsNAvg   = ((Number) nAvgField.getValue()).intValue();

        ModelSpecifications modelSpecs = new ModelSpecifications(
                session.systemId.elements, session.systemId.structure, session.systemId.model, session.engineConfig);
        CalculationSpecifications calcSpecs = new CalculationSpecifications(Property.GIBBS_ENERGY, Mode.SINGLE_POINT);
        calcSpecs.set(Parameter.TEMPERATURE, temperature);
        calcSpecs.set(Parameter.COMPOSITION, composition);
        calcSpecs.set(Parameter.MCS_L, mcsL);
        calcSpecs.set(Parameter.MCS_NEQUIL, mcsNEquil);
        calcSpecs.set(Parameter.MCS_NAVG, mcsNAvg);

        calcButton.setEnabled(false);
        abortButton.setEnabled(true);
        appCtx.clearLog();
        appCtx.log("Running " + session.engineConfig + " calculation...");
        appCtx.log("  Session: " + session.label());
        appCtx.log("  T=" + temperature + " K  composition=" + java.util.Arrays.toString(composition));
        statusSink.accept("Running " + session.engineConfig + " at T=" + temperature + " K...");

        SwingWorker<ThermodynamicResult, Object> worker = new SwingWorker<>() {
            @Override
            protected ThermodynamicResult doInBackground() throws Exception {
                // Calculation Layer Role: Unified Execution
                return service.execute(modelSpecs, calcSpecs, msg -> publish((Object) msg));
            }
            @Override
            protected void process(List<Object> chunks) {
                for (Object o : chunks) {
                    if (o instanceof String)        appCtx.log((String) o);
                    else if (o instanceof ProgressEvent) chartSink.accept((ProgressEvent) o);
                }
            }
            @Override
            protected void done() {
                calcButton.setEnabled(true);
                abortButton.setEnabled(false);
                activeWorker = null;
                if (isCancelled()) { appCtx.log("Aborted."); statusSink.accept("Aborted."); return; }
                try {
                    ThermodynamicResult r = get();
                    appCtx.log("Done. G=" + String.format("%.4f", r.gibbsEnergy) + " J/mol");
                    appCtx.log("     H=" + String.format("%.4f", r.enthalpy) + " J/mol");
                    statusSink.accept("Done — G = " + String.format("%.4f", r.gibbsEnergy) + " J/mol");
                    resultSink.accept(r);
                } catch (Exception ex) {
                    String msg = ex.getCause() != null ? ex.getCause().getMessage() : ex.getMessage();
                    appCtx.log("Error: " + msg);
                    statusSink.accept("Error: " + msg);
                }
            }
        };
        activeWorker = worker;
        worker.execute();
    }

    private static double[] parseComposition(String text, int numComponents) {
        if (text.isBlank())
            throw new IllegalArgumentException("Composition cannot be empty.");
        String[] tokens = text.split(",");
        double[] comp = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
            try { comp[i] = Double.parseDouble(tokens[i].trim()); }
            catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid value: '" + tokens[i] + "'");
            }
            if (comp[i] < 0 || comp[i] > 1)
                throw new IllegalArgumentException("Values must be 0–1: " + comp[i]);
        }
        if (comp.length == 1 && numComponents == 2) {
            double xB = comp[0];
            return new double[]{1.0 - xB, xB};
        }
        if (comp.length != numComponents)
            throw new IllegalArgumentException(
                    "Expected " + numComponents + " values, got " + comp.length);
        double sum = 0;
        for (double x : comp) sum += x;
        if (Math.abs(sum - 1.0) > 1e-8)
            throw new IllegalArgumentException("Must sum to 1.0; got " + sum);
        return comp;
    }
}
