package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.CalculationService;
import org.ce.calculation.workflow.TernaryGridScan;
import org.ce.calculation.workflow.TernaryPlotRenderer;
import org.ce.model.ModelSession;
import org.ce.model.ProgressEvent;
import org.ce.model.storage.Workspace.SystemId;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;
import java.util.function.Consumer;

/**
 * Ternary isothermal section — control panel (explorer column).
 *
 * <p>Inputs only; the rendered plot itself is shown in the shared
 * {@link OutputPanel} (right-hand results area), matching the app's
 * convention that the explorer column holds parameters and the wide output
 * panel holds results/charts.</p>
 *
 * <p>Architecture: {@link TernaryGridScan} (Java, in-process, session-cached)
 * computes the composition grid; the result is written as JSON and handed to
 * {@code scripts/isothermal_section.py} (mpltern) for rendering, since no
 * maintained Java library offers ternary contour plotting — see CLAUDE.md.
 * Java remains the single source of truth for the physics; Python only turns
 * numbers into pixels.</p>
 */
public class TernaryPlotPanel extends JPanel {

    private static final Color BG       = new Color(0x1E1E1E);
    private static final Color LABEL_FG = new Color(0xBBBBBB);

    private final WorkbenchContext          context;
    private final CalculationService        service;
    private final Consumer<String>          statusSink;
    private final Consumer<String>          logSink;
    private final Consumer<BufferedImage>   plotSink;
    private final Consumer<String>          plotErrorSink;

    private final JComboBox<String> elementsCombo  = makeEditable("Nb-Ti-V", "Nb-Ti-Zr", "Nb-V-Zr", "Ti-V-Zr");
    private final JComboBox<String> structureCombo = makeEditable("BCC_A2", "FCC_A1", "HCP_A3");
    private final JComboBox<String> modelCombo     = makeEditable("T", "T2");
    private final JLabel systemStatusLabel = new JLabel(" ");
    private final JSpinner temperatureSpinner =
            new JSpinner(new SpinnerNumberModel(1273.0, 1.0, 10000.0, 10.0));
    // "GIBBS_ENERGY", "ENTHALPY", "ENTROPY", or "SRO (1NN pair)" -- the last
    // routes through pairCombo instead of naming a Property directly.
    private static final String SRO_OPTION = "SRO (1NN pair)";
    private final JComboBox<String> quantityCombo =
            new JComboBox<>(new String[] { "GIBBS_ENERGY", "ENTHALPY", "ENTROPY", SRO_OPTION });
    private final JComboBox<String> pairCombo = new JComboBox<>();
    private final JSpinner resolutionSpinner = new JSpinner(new SpinnerNumberModel(30, 5, 100, 5));

    private final JButton runButton = new JButton("Compute & Plot");
    private final JProgressBar progressBar = new JProgressBar(0, 100);

    private SwingWorker<File, Integer> activeWorker = null;

    public TernaryPlotPanel(org.ce.CEWorkbenchContext appCtx, WorkbenchContext context,
                             Consumer<String> statusSink, Consumer<String> logSink,
                             Consumer<BufferedImage> plotSink, Consumer<String> plotErrorSink) {
        this.context       = context;
        this.service       = appCtx.getCalculationService();
        this.statusSink    = statusSink;
        this.logSink       = logSink;
        this.plotSink      = plotSink;
        this.plotErrorSink = plotErrorSink;

        setBackground(BG);
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        add(buildForm(), BorderLayout.NORTH);

        // Push identity edits into the shared context (same pattern as
        // DynamicCalculationPanel), and pull context changes made by other
        // panels back into these combos.
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

        context.addChangeListener(this::syncCombosFromContext);
        syncCombosFromContext();

        runButton.addActionListener(e -> runScan());
    }

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setOpaque(false);
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(4, 4, 4, 4);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.fill = GridBagConstraints.HORIZONTAL;

        gbc.gridy = 0; gbc.gridx = 0;
        form.add(styledLabel("Elements:"), gbc);
        gbc.gridx = 1;
        form.add(elementsCombo, gbc);

        gbc.gridy = 1; gbc.gridx = 0;
        form.add(styledLabel("Structure:"), gbc);
        gbc.gridx = 1;
        form.add(structureCombo, gbc);

        gbc.gridy = 2; gbc.gridx = 0;
        form.add(styledLabel("Model:"), gbc);
        gbc.gridx = 1;
        form.add(modelCombo, gbc);

        gbc.gridy = 3; gbc.gridx = 0; gbc.gridwidth = 2;
        styleLabel(systemStatusLabel);
        form.add(systemStatusLabel, gbc);
        gbc.gridwidth = 1;

        gbc.gridy = 4; gbc.gridx = 0;
        form.add(styledLabel("Temperature (K):"), gbc);
        gbc.gridx = 1;
        form.add(temperatureSpinner, gbc);

        gbc.gridy = 5; gbc.gridx = 0;
        form.add(styledLabel("Quantity:"), gbc);
        gbc.gridx = 1;
        form.add(quantityCombo, gbc);

        gbc.gridy = 6; gbc.gridx = 0;
        JLabel pairLabel = styledLabel("Pair:");
        form.add(pairLabel, gbc);
        gbc.gridx = 1;
        form.add(pairCombo, gbc);
        pairLabel.setVisible(false);
        pairCombo.setVisible(false);
        quantityCombo.addActionListener(e -> {
            boolean sro = SRO_OPTION.equals(quantityCombo.getSelectedItem());
            pairLabel.setVisible(sro);
            pairCombo.setVisible(sro);
        });

        gbc.gridy = 7; gbc.gridx = 0;
        form.add(styledLabel("Grid resolution (n):"), gbc);
        gbc.gridx = 1;
        form.add(resolutionSpinner, gbc);

        gbc.gridy = 8; gbc.gridx = 0; gbc.gridwidth = 2;
        form.add(runButton, gbc);

        gbc.gridy = 9;
        progressBar.setStringPainted(true);
        progressBar.setVisible(false);
        form.add(progressBar, gbc);

        JPanel wrapper = new JPanel(new BorderLayout());
        wrapper.setOpaque(false);
        wrapper.add(form, BorderLayout.NORTH);
        return wrapper;
    }

    private JLabel styledLabel(String text) {
        JLabel l = new JLabel(text);
        styleLabel(l);
        return l;
    }

    private void styleLabel(JLabel l) {
        l.setForeground(LABEL_FG);
    }

    private void pushSystemToContext() {
        String el  = editorText(elementsCombo);
        String str = editorText(structureCombo);
        String mod = editorText(modelCombo);
        if (!el.isBlank() && !str.isBlank() && !mod.isBlank()) {
            context.setSystem(el, str, mod);
        }
        updateStatus();
    }

    private void syncCombosFromContext() {
        SystemId sys = context.getSystem();
        if (sys != null) {
            setEditorText(elementsCombo,  sys.elements());
            setEditorText(structureCombo, sys.structure());
            setEditorText(modelCombo,     sys.model());
        }
        updateStatus();
    }

    private void updateStatus() {
        String el = editorText(elementsCombo);
        List<String> elements = el.isBlank() ? List.of() : List.of(el.split("-"));
        if (elements.size() != 3) {
            systemStatusLabel.setText("<html><span style='color:#F44747'>Requires exactly 3 elements"
                    + (el.isBlank() ? "" : " (got " + elements.size() + ")") + "</span></html>");
            runButton.setEnabled(false);
        } else {
            systemStatusLabel.setText(" ");
            runButton.setEnabled(true);
        }
        refreshPairCombo(elements);
    }

    private void refreshPairCombo(List<String> elements) {
        List<String> pairs = new java.util.ArrayList<>();
        for (int i = 0; i < elements.size(); i++) {
            for (int j = i + 1; j < elements.size(); j++) {
                pairs.add(elements.get(i) + "-" + elements.get(j));
            }
        }
        String current = (String) pairCombo.getSelectedItem();
        pairCombo.removeAllItems();
        for (String p : pairs) pairCombo.addItem(p);
        if (current != null && pairs.contains(current)) pairCombo.setSelectedItem(current);
    }

    private void runScan() {
        if (activeWorker != null && !activeWorker.isDone()) return;

        String elStr  = editorText(elementsCombo);
        String struct = editorText(structureCombo);
        String model  = editorText(modelCombo);
        List<String> elements = List.of(elStr.split("-"));
        if (elements.size() != 3) return;
        SystemId sys = new SystemId(elStr, struct, model);

        double temperature = (Double) temperatureSpinner.getValue();
        int n = (Integer) resolutionSpinner.getValue();

        String quantitySelection = (String) quantityCombo.getSelectedItem();
        TernaryGridScan.Quantity quantity;
        if (SRO_OPTION.equals(quantitySelection)) {
            String pairSelection = (String) pairCombo.getSelectedItem();
            if (pairSelection == null) return;
            String[] parts = pairSelection.split("-");
            quantity = new TernaryGridScan.PairSroQuantity(parts[0], parts[1]);
        } else {
            quantity = new TernaryGridScan.PropertyQuantity(Property.valueOf(quantitySelection));
        }

        runButton.setEnabled(false);
        progressBar.setVisible(true);
        progressBar.setValue(0);
        progressBar.setString("Computing grid...");

        activeWorker = new SwingWorker<>() {
            @Override
            protected File doInBackground() throws Exception {
                var modelSpecs = new org.ce.calculation.CalculationDescriptor.ModelSpecifications(
                        sys.elements(), sys.structure(), sys.model(), org.ce.model.ModelSession.EngineConfig.CVM);
                ModelSession session = service.getOrBuildSession(modelSpecs, logSink);

                Consumer<ProgressEvent> progressSink = evt -> {
                    if (evt instanceof ProgressEvent.ScanPoint sp) {
                        publish((int) (100.0 * sp.index / sp.total));
                    }
                };
                TernaryGridScan.Result result = TernaryGridScan.run(
                        service, session, elements, temperature, quantity, n, progressSink);

                return TernaryPlotRenderer.render(
                        elements, sys.structure(), sys.model(), temperature, result);
            }

            @Override
            protected void process(List<Integer> chunks) {
                if (!chunks.isEmpty()) {
                    int pct = chunks.get(chunks.size() - 1);
                    progressBar.setValue(pct);
                    progressBar.setString("Computing grid... " + pct + "%");
                }
            }

            @Override
            protected void done() {
                progressBar.setVisible(false);
                runButton.setEnabled(true);
                try {
                    File pngFile = get();
                    BufferedImage img = ImageIO.read(pngFile);
                    plotSink.accept(img);
                    if (statusSink != null) statusSink.accept("Ternary plot rendered.");
                } catch (Exception e) {
                    Throwable cause = e.getCause() != null ? e.getCause() : e;
                    plotErrorSink.accept("Plot failed: " + cause.getMessage());
                    if (logSink != null) logSink.accept("[Ternary Plot] Error: " + cause.getMessage());
                }
            }
        };
        activeWorker.execute();
    }

    // =========================================================================
    // Static helpers (mirrored from DynamicCalculationPanel)
    // =========================================================================

    private static JComboBox<String> makeEditable(String... items) {
        JComboBox<String> cb = new JComboBox<>(items);
        cb.setEditable(true);
        cb.setPreferredSize(new Dimension(120, 24));
        return cb;
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
}
