package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.CalculationRegistry;
import org.ce.calculation.CalculationSpecifications;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.ProgressEvent;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Metadata-driven calculation explorer.
 * Dynamically builds its input form based on CalculationRegistry requirements.
 */
public class DynamicCalculationPanel extends JPanel {

    private static final Color BG       = new Color(0x252526);
    private static final Color LABEL_FG = new Color(0xCCCCCC);

    private final org.ce.CEWorkbenchContext appCtx;
    private final WorkbenchContext          context;
    private final CalculationService        service;
    private final Consumer<String>          statusSink;
    private final Consumer<ThermodynamicResult> resultSink;
    private final Consumer<ProgressEvent>   chartSink;

    private final JComboBox<Property> propertyCombo = new JComboBox<>();
    private final JComboBox<Mode>     modeCombo     = new JComboBox<>();
    private final JPanel              parameterForm = new JPanel(new GridBagLayout());
    private final Map<Parameter, JComponent> parameterFields = new HashMap<>();

    private final JButton runButton   = new JButton("Run Calculation");
    private final JButton abortButton = new JButton("Abort");
    private SwingWorker<?, ?> activeWorker = null;

    public DynamicCalculationPanel(org.ce.CEWorkbenchContext appCtx,
                                   WorkbenchContext context,
                                   Consumer<String> statusSink,
                                   Consumer<ThermodynamicResult> resultSink,
                                   Consumer<ProgressEvent> chartSink,
                                   QuantityDescriptor.SelectionModel quantityModel) {
        this.appCtx     = appCtx;
        this.context    = context;
        this.service    = appCtx.getCalculationService();
        this.statusSink = statusSink;
        this.resultSink = resultSink;
        this.chartSink  = chartSink;

        setBackground(BG);
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));

        // ── Observers ────────────────────────────────────────────────────────
        context.addSessionListener(this::onSessionChanged);

        propertyCombo.addActionListener(e -> rebuildModeOptions());
        modeCombo.addActionListener(e -> rebuildParameterForm());

        // ── Header (Selection) ──────────────────────────────────────────────
        JPanel header = new JPanel(new GridBagLayout());
        header.setBackground(BG);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 1.0;
        gbc.insets = new Insets(0, 0, 4, 0);

        gbc.gridy = 0;
        header.add(createLabel("Property"), gbc);
        gbc.gridy = 1;
        header.add(propertyCombo, gbc);

        gbc.gridy = 2;
        header.add(createLabel("Mode"), gbc);
        gbc.gridy = 3;
        header.add(modeCombo, gbc);

        // ── Quantities Header ────────────────────────────────────────────────
        JPanel topZone = new JPanel(new BorderLayout());
        topZone.setBackground(BG);
        topZone.add(new QuantitiesPanel(quantityModel, context), BorderLayout.NORTH);
        topZone.add(header, BorderLayout.CENTER);

        // ── Parameter Form ──────────────────────────────────────────────────
        parameterForm.setBackground(BG);

        // ── Buttons ──────────────────────────────────────────────────────────
        runButton.addActionListener(e -> startExecution());
        abortButton.addActionListener(e -> { if (activeWorker != null) activeWorker.cancel(true); });
        abortButton.setEnabled(false);
        abortButton.setForeground(new Color(0xF44747));

        JPanel buttonPanel = new JPanel(new GridLayout(1, 2, 4, 0));
        buttonPanel.setBackground(BG);
        buttonPanel.add(runButton);
        buttonPanel.add(abortButton);

        // ── Assemble ─────────────────────────────────────────────────────────
        add(topZone, BorderLayout.NORTH);
        add(new JScrollPane(parameterForm), BorderLayout.CENTER);
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
        rebuildModeOptions();
    }

    private void rebuildModeOptions() {
        Property prop = (Property) propertyCombo.getSelectedItem();
        if (prop == null) return;

        ModelSession session = context.getActiveSession();
        ModelSession.EngineConfig engine = (session != null) 
                ? session.engineConfig 
                : ModelSession.EngineConfig.cvm();

        Mode currentSelected = (Mode) modeCombo.getSelectedItem();
        modeCombo.removeAllItems();
        List<Mode> modes = CalculationRegistry.INSTANCE.getAvailableModes(prop, engine);
        for (Mode m : modes) modeCombo.addItem(m);
        
        if (currentSelected != null && modes.contains(currentSelected)) {
            modeCombo.setSelectedItem(currentSelected);
        }
        rebuildParameterForm();
    }

    private void rebuildParameterForm() {
        parameterForm.removeAll();
        parameterFields.clear();

        Property prop = (Property) propertyCombo.getSelectedItem();
        Mode mode = (Mode) modeCombo.getSelectedItem();
        if (prop == null || mode == null) return;

        ModelSession session = context.getActiveSession();
        ModelSession.EngineConfig engine = (session != null) 
                ? session.engineConfig 
                : ModelSession.EngineConfig.cvm();

        List<Parameter> requirements = CalculationRegistry.INSTANCE.getRequirements(prop, mode, engine);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 1.0;
        gbc.insets = new Insets(4, 0, 0, 0);
        int row = 0;

        for (Parameter p : requirements) {
            gbc.gridy = row++;
            parameterForm.add(createLabel(p.name), gbc);
            
            JComponent editor = ParameterFieldFactory.createEditor(p, null);
            parameterFields.put(p, editor);
            
            gbc.gridy = row++;
            parameterForm.add(editor, gbc);
        }

        parameterForm.revalidate();
        parameterForm.repaint();
    }

    private JLabel createLabel(String text) {
        JLabel label = new JLabel(text);
        label.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        label.setForeground(LABEL_FG);
        return label;
    }

    private void updateButtonState() {
        runButton.setEnabled(context.hasActiveSession());
    }

    private void startExecution() {
        ModelSession session = context.getActiveSession();
        if (session == null) return;

        Property property = (Property) propertyCombo.getSelectedItem();
        Mode mode = (Mode) modeCombo.getSelectedItem();
        
        CalculationSpecifications specs = new CalculationSpecifications(property, mode);
        for (Map.Entry<Parameter, JComponent> entry : parameterFields.entrySet()) {
            specs.set(entry.getKey(), ParameterFieldFactory.getValue(entry.getKey(), entry.getValue()));
        }

        ModelSpecifications modelSpecs = new ModelSpecifications(
                session.systemId.elements, session.systemId.structure, session.systemId.model, session.engineConfig);

        activeWorker = new SwingWorker<Object, String>() {
            @Override
            protected Object doInBackground() throws Exception {
                publish("Executing " + property.displayName + " in " + mode.displayName + " mode...");
                if (mode == Mode.GRID_SCAN) {
                    return service.executeGridScan(modelSpecs, specs, this::publish);
                } else if (mode == Mode.LINE_SCAN) {
                    return service.executeScan(modelSpecs, specs, this::publish);
                } else {
                    return service.execute(modelSpecs, specs, this::publish);
                }
            }

            @Override
            protected void process(List<String> chunks) {
                chunks.forEach(appCtx::log);
            }

            @Override
            protected void done() {
                try {
                    Object result = get();
                    if (result instanceof ThermodynamicResult) {
                        resultSink.accept((ThermodynamicResult) result);
                    } else if (result instanceof List) {
                        // Scan result handled by logSink or specialized resultSink eventually
                    }
                    statusSink.accept("Calculation completed.");
                } catch (Exception e) {
                    if (!isCancelled()) {
                        String msg = "Calculation failed: " + e.getMessage();
                        appCtx.log(msg);
                        statusSink.accept(msg);
                    } else {
                        statusSink.accept("Calculation aborted.");
                    }
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
}
