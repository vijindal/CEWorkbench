package org.ce.ui.gui;

import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;
import org.ce.calculation.CalculationRegistry;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationSpecifications;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Metadata-driven calculation panel that replaces legacy specific panels.
 * Dynamically builds input forms based on the selected Property and Mode.
 * Simplified to support multicomponent systems and standard thermodynamic scans.
 */
public class DynamicCalculationPanel extends JPanel {

    private static final Color BG        = new Color(0x1E1E1E);
    private static final Color LABEL_FG  = new Color(0xBBBBBB);

    private final org.ce.CEWorkbenchContext appCtx;
    private final WorkbenchContext          context;
    private final CalculationService        service;
    private final Consumer<String>          statusSink;
    private final Consumer<String>          logSink;
    private final Consumer<Object>          resultSink;
    private final Consumer<ProgressEvent>   chartSink;

    private final JComboBox<Property> propertyCombo = new JComboBox<>();
    private final JPanel              parameterForm = new JPanel(new GridBagLayout());
    private final Map<Parameter, JComponent> parameterFields = new HashMap<>();
    private final Map<String, JSpinner[]>    compSpinners    = new HashMap<>();

    private final JButton runButton   = new JButton("Run Calculation");
    private final JButton abortButton = new JButton("Abort");
    private SwingWorker<?, ?> activeWorker = null;

    public DynamicCalculationPanel(org.ce.CEWorkbenchContext appCtx,
                                   WorkbenchContext context,
                                   Consumer<String> statusSink,
                                   Consumer<String> logSink,
                                   Consumer<Object> resultSink,
                                   Consumer<ProgressEvent> chartSink,
                                   QuantityDescriptor.SelectionModel quantityModel) {
        this.appCtx     = appCtx;
        this.context    = context;
        this.service    = appCtx.getCalculationService();
        this.statusSink = statusSink;
        this.logSink    = logSink;
        this.resultSink = resultSink;
        this.chartSink  = chartSink;

        setBackground(BG);
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        // ── Observers ────────────────────────────────────────────────────────
        context.addSessionListener(this::onSessionChanged);

        propertyCombo.addActionListener(e -> {
            rebuildParameterForm();
            updateSelectionModel(quantityModel);
        });

        // ── Parameter Form (Scrollable) ──────────────────────────────────────
        parameterForm.setBackground(BG);
        JScrollPane scroll = new JScrollPane(parameterForm);
        scroll.setBorder(null);
        scroll.setOpaque(false);
        scroll.getViewport().setOpaque(false);

        // ── Buttons ──────────────────────────────────────────────────────────
        runButton.addActionListener(e -> startExecution());
        abortButton.addActionListener(e -> { if (activeWorker != null) activeWorker.cancel(true); });
        abortButton.setEnabled(false);
        abortButton.setForeground(new Color(0xF44747));

        JPanel buttonPanel = new JPanel(new GridLayout(1, 2, 6, 0));
        buttonPanel.setBackground(BG);
        buttonPanel.add(runButton);
        buttonPanel.add(abortButton);
        buttonPanel.setBorder(BorderFactory.createEmptyBorder(10, 0, 0, 0));

        // ── Assemble ─────────────────────────────────────────────────────────
        add(scroll, BorderLayout.CENTER);
        add(buttonPanel, BorderLayout.SOUTH);

        // Initial populate
        onSessionChanged(context.getActiveSession());
    }

    private void onSessionChanged(ModelSession session) {
        ModelSession.EngineConfig engine = (session != null) 
                ? session.engineConfig 
                : ModelSession.EngineConfig.cvm();

        rebuildPropertyOptions(engine);
        updateButtonState();
    }

    private void rebuildPropertyOptions(ModelSession.EngineConfig engine) {
        propertyCombo.removeAllItems();
        List<Property> props = CalculationRegistry.INSTANCE.getAvailableProperties(engine);
        for (Property p : props) propertyCombo.addItem(p);
        rebuildParameterForm();
    }

    private void rebuildParameterForm() {
        parameterForm.removeAll();
        parameterFields.clear();
        compSpinners.clear();

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 1.0;
        gbc.gridx = 0;
        int row = 0;

        // ── Header: Property Selection ──
        gbc.insets = new Insets(0, 0, 4, 0);
        gbc.gridy = row++;
        parameterForm.add(createLabel("Property"), gbc);
        gbc.gridy = row++;
        parameterForm.add(propertyCombo, gbc);

        // Spacer before parameters
        gbc.insets = new Insets(12, 0, 0, 0);

        Property prop = (Property) propertyCombo.getSelectedItem();
        if (prop == null) {
            updateUIState();
            return;
        }

        // Default to ANALYSIS mode
        Mode mode = Mode.ANALYSIS;
        ModelSession session = context.getActiveSession();
        ModelSession.EngineConfig engine = (session != null) ? session.engineConfig : ModelSession.EngineConfig.cvm();
        List<Parameter> requirements = CalculationRegistry.INSTANCE.getRequirements(prop, mode, engine);

        // ── Parameter Rows ──
        int i = 0;
        while (i < requirements.size()) {
            Parameter p = requirements.get(i);
            gbc.gridy = row++;
            
            if (p == Parameter.T_START && i + 2 < requirements.size() &&
                requirements.get(i+1) == Parameter.T_END && requirements.get(i+2) == Parameter.T_STEP) {
                
                parameterForm.add(createRangeRow("T (Temperature)", Parameter.T_START, Parameter.T_END, Parameter.T_STEP), gbc);
                i += 3;
            } else if (p == Parameter.X_STARTS && i + 2 < requirements.size() &&
                       requirements.get(i+1) == Parameter.X_ENDS && requirements.get(i+2) == Parameter.X_STEPS) {
                
                if (session != null) {
                    String[] elements = session.systemId.elements.split("-");
                    for (int n = 1; n < elements.length; n++) {
                        parameterForm.add(createCompRangeRow(elements[n]), gbc);
                        if (n < elements.length - 1) gbc.gridy = row++;
                    }
                } else {
                    parameterForm.add(createCompRangeRow("Composition (X)"), gbc);
                }
                i += 3;
            } else {
                JPanel pnl = new JPanel(new BorderLayout(0, 2));
                pnl.setOpaque(false);
                pnl.add(createLabel(p.name), BorderLayout.NORTH);
                JComponent editor = ParameterFieldFactory.createEditor(p, null);
                parameterFields.put(p, editor);
                pnl.add(editor, BorderLayout.CENTER);
                parameterForm.add(pnl, gbc);
                i++;
            }
            gbc.insets = new Insets(10, 0, 0, 0); // Spacing for subsequent rows
        }

        // Vertical spacer at bottom
        gbc.gridy = row++;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        parameterForm.add(new JPanel() {{ setOpaque(false); }}, gbc);

        updateUIState();
    }

    private JPanel createRangeRow(String label, Parameter start, Parameter end, Parameter step) {
        JPanel row = new JPanel(new BorderLayout(0, 2));
        row.setOpaque(false);
        row.add(createLabel(label), BorderLayout.NORTH);

        JPanel grid = new JPanel(new GridLayout(1, 3, 4, 0));
        grid.setOpaque(false);
        grid.add(createMiniField(start));
        grid.add(createMiniField(end));
        grid.add(createMiniField(step));
        
        row.add(grid, BorderLayout.CENTER);
        return row;
    }

    private JPanel createCompRangeRow(String label) {
        JPanel row = new JPanel(new BorderLayout(0, 2));
        row.setOpaque(false);
        row.add(createLabel(label), BorderLayout.NORTH);

        JSpinner sStart = (JSpinner) ParameterFieldFactory.createEditor(Parameter.X_START, null);
        JSpinner sEnd   = (JSpinner) ParameterFieldFactory.createEditor(Parameter.X_END, null);
        JSpinner sStep  = (JSpinner) ParameterFieldFactory.createEditor(Parameter.X_STEP, null);
        compSpinners.put(label, new JSpinner[]{ sStart, sEnd, sStep });

        JPanel grid = new JPanel(new GridLayout(1, 3, 4, 0));
        grid.setOpaque(false);
        grid.add(createCustomMiniField("Start", sStart));
        grid.add(createCustomMiniField("End", sEnd));
        grid.add(createCustomMiniField("Step", sStep));
        
        row.add(grid, BorderLayout.CENTER);
        return row;
    }

    private JPanel createCustomMiniField(String label, JSpinner spinner) {
        JPanel pnl = new JPanel(new BorderLayout(0, 1));
        pnl.setOpaque(false);
        JLabel l = new JLabel(label);
        l.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 8));
        l.setForeground(new Color(0x888888));
        l.setHorizontalAlignment(SwingConstants.CENTER);
        pnl.add(l, BorderLayout.NORTH);
        spinner.setMinimumSize(new Dimension(30, 22));
        pnl.add(spinner, BorderLayout.CENTER);
        return pnl;
    }

    private JPanel createMiniField(Parameter p) {
        JSpinner spinner = (JSpinner) ParameterFieldFactory.createEditor(p, null);
        parameterFields.put(p, spinner);
        return createCustomMiniField(p.name.split(" ")[1], spinner);
    }

    private JLabel createLabel(String text) {
        JLabel label = new JLabel(text);
        label.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        label.setForeground(LABEL_FG);
        return label;
    }

    private void updateUIState() {
        runButton.setEnabled(context.hasActiveSession());
        parameterForm.revalidate();
        parameterForm.repaint();
    }

    private void updateButtonState() {
        runButton.setEnabled(context.hasActiveSession());
    }

    private void startExecution() {
        ModelSession session = context.getActiveSession();
        if (session == null) return;

        Property prop = (Property) propertyCombo.getSelectedItem();
        Mode mode = Mode.ANALYSIS;

        ModelSpecifications modelSpecs = new ModelSpecifications(
                session.systemId.elements,
                session.systemId.structure,
                session.systemId.model,
                session.engineConfig);

        CalculationSpecifications specs = new CalculationSpecifications(prop, mode);
        
        // Map individual fields
        parameterFields.forEach((p, editor) -> {
            if (editor instanceof JSpinner) {
                specs.set(p, ((JSpinner) editor).getValue());
            }
        });

        // Map component arrays
        if (compSpinners.size() > 0) {
            double[] starts = new double[compSpinners.size()];
            double[] ends   = new double[compSpinners.size()];
            double[] steps  = new double[compSpinners.size()];
            
            String[] elements = session.systemId.elements.split("-");
            for (int n = 1; n < elements.length; n++) {
                JSpinner[] spinners = compSpinners.get(elements[n]);
                if (spinners != null) {
                    starts[n-1] = (Double) spinners[0].getValue();
                    ends[n-1]   = (Double) spinners[1].getValue();
                    steps[n-1]  = (Double) spinners[2].getValue();
                }
            }
            specs.set(Parameter.X_STARTS, starts);
            specs.set(Parameter.X_ENDS, ends);
            specs.set(Parameter.X_STEPS, steps);
        }

        activeWorker = new SwingWorker<Object, String>() {
            @Override
            protected Object doInBackground() throws Exception {
                return service.execute(modelSpecs, specs, this::publish, chartSink);
            }

            @Override
            protected void process(List<String> chunks) {
                chunks.forEach(logSink::accept);
            }

            @Override
            protected void done() {
                try {
                    Object result = get();
                    if (result != null) resultSink.accept(result);
                } catch (Exception e) {
                    statusSink.accept("Error: " + e.getLocalizedMessage());
                } finally {
                    runButton.setEnabled(true);
                    abortButton.setEnabled(false);
                    activeWorker = null;
                }
            }
        };

        runButton.setEnabled(false);
        abortButton.setEnabled(true);
        activeWorker.execute();
    }

    private void updateSelectionModel(QuantityDescriptor.SelectionModel model) {
        Property prop = (Property) propertyCombo.getSelectedItem();
        if (prop != null) {
            try {
                org.ce.calculation.QuantityDescriptor q = org.ce.calculation.QuantityDescriptor.valueOf(prop.name());
                model.setExclusive(q);
            } catch (Exception ignored) {}
        }
    }
}
