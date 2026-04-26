package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.ModelSpecifications;
import org.ce.model.ProgressEvent;
import org.ce.model.ModelSession;
import org.ce.model.mcs.McsSuggester;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.workflow.CalculationService;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Metadata-driven calculation panel. Hosts system identity inputs (elements, structure,
 * model, engine) and dynamically builds parameter forms based on selected Property/Mode.
 * Automatically builds the session on Run — no separate Rebuild step required.
 */
public class DynamicCalculationPanel extends JPanel {

    private static final Color BG        = new Color(0x1E1E1E);
    private static final Color LABEL_FG  = new Color(0xBBBBBB);

    // ── Dot colours ───────────────────────────────────────────────────────────
    private static final Color DOT_READY = new Color(0x4EC9B0);
    private static final Color DOT_BUSY  = new Color(0xDCDCAA);
    private static final Color DOT_NONE  = new Color(0xF44747);

    private final WorkbenchContext          context;
    private final CalculationService        service;
    private final Consumer<String>          statusSink;
    private final Consumer<String>          logSink;
    private final Consumer<Object>          resultSink;
    private final Consumer<ProgressEvent>   chartSink;

    // ── System identity (migrated from SessionBar) ────────────────────────────
    private final JComboBox<String> elementsCombo  = makeEditable("Nb-Ti", "Cu-Au", "Al-Ti", "Nb-Mo");
    private final JComboBox<String> structureCombo = makeEditable("BCC_A2", "FCC_A1", "HCP_A3");
    private final JComboBox<String> modelCombo     = makeEditable("T", "T2");
    private final JComboBox<String> engineCombo    = new JComboBox<>(new String[]{"CVM", "MCS"});

    // ── Status indicator ──────────────────────────────────────────────────────
    private final DotIcon sessionDot   = new DotIcon(DOT_NONE);
    private final JLabel  sessionLabel = new JLabel("No session");
    private boolean       buildingSession = false;

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
        this.context    = context;
        this.service    = appCtx.getCalculationService();
        this.statusSink = statusSink;
        this.logSink    = logSink;
        this.resultSink = resultSink;
        this.chartSink  = chartSink;

        setBackground(BG);
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        // ── System header ─────────────────────────────────────────────────────
        add(buildSystemPanel(), BorderLayout.NORTH);

        // ── Wire combo listeners → context.setSystem() ────────────────────────
        DocumentListener pushSystem = new DocumentListener() {
            public void insertUpdate(DocumentEvent e)  { pushSystemToContext(); }
            public void removeUpdate(DocumentEvent e)  { pushSystemToContext(); }
            public void changedUpdate(DocumentEvent e) { pushSystemToContext(); }
        };
        editorDoc(elementsCombo).addDocumentListener(pushSystem);
        editorDoc(structureCombo).addDocumentListener(pushSystem);
        editorDoc(modelCombo).addDocumentListener(pushSystem);
        elementsCombo.addActionListener(e -> pushSystemToContext());
        structureCombo.addActionListener(e -> pushSystemToContext());
        modelCombo.addActionListener(e -> pushSystemToContext());

        // Engine change invalidates session status dot and rebuilds property options
        engineCombo.addActionListener(e -> onSessionChanged(null));

        // Pull system changes from other panels (e.g. CECManagementPanel)
        context.addChangeListener(this::syncCombosFromContext);

        // ── Session observer ──────────────────────────────────────────────────
        context.addSessionListener(this::onSessionChanged);

        propertyCombo.addActionListener(e -> {
            rebuildParameterForm();
            updateSelectionModel(quantityModel);
        });

        // ── Parameter Form (Scrollable) ──────────────────────────────────────
        parameterForm.setBackground(BG);
        JPanel scrollContent = new JPanel(new BorderLayout());
        scrollContent.setOpaque(false);
        scrollContent.setBackground(BG);
        scrollContent.add(parameterForm, BorderLayout.NORTH);
        JScrollPane scroll = new JScrollPane(scrollContent,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scroll.setBorder(null);
        scroll.setOpaque(false);
        scroll.getViewport().setOpaque(false);
        scroll.getViewport().setBackground(BG);

        // ── Buttons ──────────────────────────────────────────────────────────
        runButton.addActionListener(e -> startExecution());
        abortButton.addActionListener(e -> { if (activeWorker != null) activeWorker.cancel(true); });
        abortButton.setEnabled(false);
        abortButton.setForeground(new Color(0xF44747));
        runButton.setPreferredSize(new Dimension(0, 32));
        abortButton.setPreferredSize(new Dimension(0, 32));

        JPanel buttonPanel = new JPanel(new GridLayout(1, 2, 6, 0));
        buttonPanel.setBackground(BG);
        buttonPanel.add(runButton);
        buttonPanel.add(abortButton);
        buttonPanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(1, 0, 0, 0, new Color(0x3C3C3C)),
            BorderFactory.createEmptyBorder(8, 0, 0, 0)
        ));

        // ── Assemble ─────────────────────────────────────────────────────────
        add(scroll, BorderLayout.CENTER);
        add(buttonPanel, BorderLayout.SOUTH);

        // Initial state
        pushSystemToContext();
        syncCombosFromContext();
        updateButtonState();
    }

    // =========================================================================
    // System panel
    // =========================================================================

    private JPanel buildSystemPanel() {
        elementsCombo.setPreferredSize(null);
        structureCombo.setPreferredSize(null);
        modelCombo.setPreferredSize(null);
        engineCombo.setPreferredSize(null);

        JPanel p = new JPanel(new GridBagLayout());
        p.setBackground(new Color(0x252526));
        p.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(0x3C3C3C)),
            BorderFactory.createEmptyBorder(6, 8, 8, 8)
        ));

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill      = GridBagConstraints.HORIZONTAL;
        gbc.weightx   = 1.0;
        gbc.gridx     = 0;
        gbc.gridwidth = GridBagConstraints.REMAINDER;
        int row = 0;

        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 2, 0);
        p.add(sectionLabel("ELEMENTS"), gbc);
        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 6, 0);
        p.add(elementsCombo, gbc);

        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 2, 0);
        p.add(sectionLabel("STRUCTURE"), gbc);
        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 6, 0);
        p.add(structureCombo, gbc);

        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 2, 0);
        p.add(sectionLabel("MODEL"), gbc);
        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 6, 0);
        p.add(modelCombo, gbc);

        gbc.gridy = row++; gbc.insets = new Insets(0, 0, 2, 0);
        p.add(sectionLabel("ENGINE"), gbc);

        gbc.gridy = row; gbc.insets = new Insets(0, 0, 0, 0);
        JPanel engineRow = new JPanel(new BorderLayout(6, 0));
        engineRow.setOpaque(false);
        engineRow.add(engineCombo, BorderLayout.CENTER);
        sessionLabel.setForeground(new Color(0x9CDCFE));
        sessionLabel.setFont(sessionLabel.getFont().deriveFont(11f));
        JPanel statusChip = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        statusChip.setOpaque(false);
        statusChip.add(new JLabel(sessionDot));
        statusChip.add(sessionLabel);
        engineRow.add(statusChip, BorderLayout.EAST);
        p.add(engineRow, gbc);

        return p;
    }

    // =========================================================================
    // Context push/pull
    // =========================================================================

    private void pushSystemToContext() {
        String el  = editorText(elementsCombo);
        String str = editorText(structureCombo);
        String mod = editorText(modelCombo);
        if (!el.isBlank() && !str.isBlank() && !mod.isBlank()) {
            context.setSystem(el, str, mod);
        }
        updateButtonState();
    }

    private void syncCombosFromContext() {
        SystemId sys = context.getSystem();
        if (sys == null) return;
        setEditorText(elementsCombo,  sys.elements());
        setEditorText(structureCombo, sys.structure());
        setEditorText(modelCombo,     sys.model());
    }

    // =========================================================================
    // Session observer
    // =========================================================================

    private void onSessionChanged(ModelSession session) {
        ModelSession.EngineConfig engine = (session != null)
                ? session.engineConfig
                : ModelSession.EngineConfig.valueOf((String) engineCombo.getSelectedItem());

        if (session != null) {
            sessionDot.setColor(DOT_READY);
            sessionLabel.setText("Ready: " + session.label());
            engineCombo.setSelectedItem(session.engineConfig.name());
        } else if (buildingSession) {
            sessionDot.setColor(DOT_BUSY);
            sessionLabel.setText("Building…");
        } else {
            sessionDot.setColor(DOT_NONE);
            sessionLabel.setText("No session");
        }
        repaint();
        rebuildPropertyOptions(engine);
        updateButtonState();
    }

    // =========================================================================
    // Property / parameter form
    // =========================================================================

    private void rebuildPropertyOptions(ModelSession.EngineConfig engine) {
        propertyCombo.removeAllItems();
        List<Property> props = Registry.getAvailableProperties(engine);
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

        gbc.insets = new Insets(12, 0, 0, 0);

        Property prop = (Property) propertyCombo.getSelectedItem();
        if (prop == null) {
            updateUIState();
            return;
        }

        Mode mode = Mode.ANALYSIS;
        ModelSession session = context.getActiveSession();
        ModelSession.EngineConfig engine = (session != null)
                ? session.engineConfig
                : ModelSession.EngineConfig.valueOf((String) engineCombo.getSelectedItem());
        List<Parameter> requirements = Registry.getRequirements(prop, mode, engine);

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

                String elemStr = (session != null)
                        ? session.systemId.elements()
                        : editorText(elementsCombo);
                String[] elements = elemStr.split("-");
                for (int n = 1; n < elements.length; n++) {
                    parameterForm.add(createCompRangeRow(elements[n]), gbc);
                    if (n < elements.length - 1) gbc.gridy = row++;
                }
                i += 3;
            } else if (p == Parameter.MCS_L && i + 2 < requirements.size() &&
                       requirements.get(i+1) == Parameter.MCS_NEQUIL && requirements.get(i+2) == Parameter.MCS_NAVG) {

                parameterForm.add(createMcsGroup(), gbc);
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
            gbc.insets = new Insets(10, 0, 0, 0);
        }

        // Vertical spacer at bottom
        gbc.gridy = row++;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
        parameterForm.add(new JPanel() {{ setOpaque(false); }}, gbc);

        updateUIState();
    }

    // =========================================================================
    // Execution
    // =========================================================================

    private void startExecution() {
        String el  = editorText(elementsCombo);
        String str = editorText(structureCombo);
        String mod = editorText(modelCombo);
        String eng = (String) engineCombo.getSelectedItem();

        if (el.isBlank() || str.isBlank() || mod.isBlank()) {
            logSink.accept("Error: Elements, Structure, and Model must be specified.");
            return;
        }

        Property prop = (Property) propertyCombo.getSelectedItem();
        Mode mode = Mode.ANALYSIS;

        ModelSpecifications modelSpecs = new ModelSpecifications(
                el, str, mod, ModelSession.EngineConfig.valueOf(eng));

        JobSpecifications specs = new JobSpecifications(prop, mode);

        parameterFields.forEach((p, editor) -> {
            if (editor instanceof JSpinner) {
                specs.set(p, ((JSpinner) editor).getValue());
            }
        });

        if (!compSpinners.isEmpty()) {
            double[] starts = new double[compSpinners.size()];
            double[] ends   = new double[compSpinners.size()];
            double[] steps  = new double[compSpinners.size()];

            String[] elements = el.split("-");
            for (int n = 1; n < elements.length; n++) {
                JSpinner[] spinners = compSpinners.get(elements[n]);
                if (spinners != null) {
                    starts[n-1] = (Double) spinners[0].getValue();
                    ends[n-1]   = (Double) spinners[1].getValue();
                    steps[n-1]  = (Double) spinners[2].getValue();
                }
            }
            specs.set(Parameter.X_STARTS, starts);
            specs.set(Parameter.X_ENDS,   ends);
            specs.set(Parameter.X_STEPS,  steps);
        }

        buildingSession = true;
        sessionDot.setColor(DOT_BUSY);
        sessionLabel.setText("Building…");
        runButton.setEnabled(false);
        abortButton.setEnabled(true);
        repaint();

        activeWorker = new SwingWorker<Object, String>() {
            @Override
            protected Object doInBackground() throws Exception {
                // service.execute() calls getOrBuildSession() internally —
                // cluster identification runs once per unique (el, str, mod, eng)
                return service.execute(modelSpecs, specs, this::publish, chartSink);
            }

            @Override
            protected void process(List<String> chunks) {
                chunks.forEach(logSink::accept);
            }

            @Override
            protected void done() {
                buildingSession = false;
                activeWorker = null;
                try {
                    Object result = get();
                    ModelSession built = service.getLastCachedSession();
                    if (built != null) context.setActiveSession(built);
                    if (result != null) resultSink.accept(result);
                } catch (Exception e) {
                    sessionDot.setColor(DOT_NONE);
                    sessionLabel.setText("Build/run failed");
                    repaint();
                    statusSink.accept("Error: " + e.getLocalizedMessage());
                    logSink.accept("Error: " + e.getLocalizedMessage());
                } finally {
                    updateButtonState();
                    abortButton.setEnabled(false);
                }
            }
        };
        activeWorker.execute();
    }

    // =========================================================================
    // UI helpers
    // =========================================================================

    private void updateUIState() {
        updateButtonState();
        parameterForm.revalidate();
        parameterForm.repaint();
    }

    private void updateButtonState() {
        String el  = editorText(elementsCombo);
        String str = editorText(structureCombo);
        String mod = editorText(modelCombo);
        boolean systemValid = !el.isBlank() && !str.isBlank() && !mod.isBlank();
        runButton.setEnabled(systemValid && !buildingSession && activeWorker == null);
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
        grid.add(createCustomMiniField("End",   sEnd));
        grid.add(createCustomMiniField("Step",  sStep));

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
        spinner.setMaximumSize(new Dimension(Integer.MAX_VALUE, 22));
        spinner.setPreferredSize(new Dimension(1, 22));
        pnl.add(spinner, BorderLayout.CENTER);
        return pnl;
    }

    private JPanel createMiniField(Parameter p) {
        JSpinner spinner = (JSpinner) ParameterFieldFactory.createEditor(p, null);
        parameterFields.put(p, spinner);
        return createCustomMiniField(p.name.split(" ")[1], spinner);
    }

    private JPanel createMcsGroup() {
        McsSuggester.Suggestion def = McsSuggester.defaultSuggestion();

        JSpinner lSpinner     = (JSpinner) ParameterFieldFactory.createEditor(Parameter.MCS_L,      null);
        JSpinner equilSpinner = (JSpinner) ParameterFieldFactory.createEditor(Parameter.MCS_NEQUIL, null);
        JSpinner avgSpinner   = (JSpinner) ParameterFieldFactory.createEditor(Parameter.MCS_NAVG,   null);

        lSpinner.setValue(def.L());
        equilSpinner.setValue(def.nEquil());
        avgSpinner.setValue(def.nAvg());

        parameterFields.put(Parameter.MCS_L,      lSpinner);
        parameterFields.put(Parameter.MCS_NEQUIL, equilSpinner);
        parameterFields.put(Parameter.MCS_NAVG,   avgSpinner);

        JLabel hintLabel = new JLabel(McsSuggester.hint(def.L()));
        hintLabel.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
        hintLabel.setForeground(new Color(0x888888));

        // When L changes, update equil/avg suggestions
        lSpinner.addChangeListener(e -> {
            int L = (Integer) lSpinner.getValue();
            McsSuggester.Suggestion s = McsSuggester.suggest(L);
            equilSpinner.setValue(s.nEquil());
            avgSpinner.setValue(s.nAvg());
            hintLabel.setText(McsSuggester.hint(L));
        });

        JPanel outer = new JPanel(new GridBagLayout());
        outer.setOpaque(false);
        GridBagConstraints g = new GridBagConstraints();
        g.fill    = GridBagConstraints.HORIZONTAL;
        g.weightx = 1.0;
        g.gridx   = 0;

        g.gridy  = 0; g.insets = new Insets(0, 0, 2, 0);
        outer.add(createLabel("MCS Parameters"), g);

        // L / Equil / Avg in one row
        JPanel row = new JPanel(new GridLayout(1, 3, 4, 0));
        row.setOpaque(false);
        row.add(createCustomMiniField("L",     lSpinner));
        row.add(createCustomMiniField("Equil", equilSpinner));
        row.add(createCustomMiniField("Avg",   avgSpinner));

        g.gridy  = 1; g.insets = new Insets(0, 0, 3, 0);
        outer.add(row, g);

        g.gridy  = 2; g.insets = new Insets(0, 0, 0, 0);
        outer.add(hintLabel, g);

        return outer;
    }

    private JLabel createLabel(String text) {
        JLabel label = new JLabel(text);
        label.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        label.setForeground(LABEL_FG);
        return label;
    }

    // =========================================================================
    // Static helpers (mirrored from SessionBar)
    // =========================================================================

    private static JComboBox<String> makeEditable(String... items) {
        JComboBox<String> cb = new JComboBox<>(items);
        cb.setEditable(true);
        cb.setPreferredSize(new Dimension(90, 24));
        return cb;
    }

    private static JLabel sectionLabel(String text) {
        JLabel l = new JLabel(text);
        l.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        l.setForeground(new Color(0x858585));
        return l;
    }

    private static javax.swing.text.Document editorDoc(JComboBox<String> cb) {
        return ((JTextField) cb.getEditor().getEditorComponent()).getDocument();
    }

    private static String editorText(JComboBox<String> cb) {
        return ((JTextField) cb.getEditor().getEditorComponent()).getText().trim();
    }

    private static void setEditorText(JComboBox<String> cb, String text) {
        JTextField tf = (JTextField) cb.getEditor().getEditorComponent();
        if (!tf.getText().equals(text)) tf.setText(text);
    }

    // =========================================================================
    // Status dot
    // =========================================================================

    private static final class DotIcon implements javax.swing.Icon {
        private Color color;
        DotIcon(Color c) { this.color = c; }
        void setColor(Color c) { this.color = c; }
        @Override public int getIconWidth()  { return 12; }
        @Override public int getIconHeight() { return 12; }
        @Override public void paintIcon(Component c, Graphics g, int x, int y) {
            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setColor(color);
            g2.fillOval(x + 1, y + 1, 10, 10);
            g2.dispose();
        }
    }
}
