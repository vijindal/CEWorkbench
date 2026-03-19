package org.ce.ui.gui;

import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.SystemId;
import org.ce.workflow.CalculationService;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.util.function.Consumer;

/**
 * Parameter panel for single-point thermodynamic calculations (shown in the Explorer column).
 *
 * <p>Elements, Structure, and Model sync with the shared {@link WorkbenchContext}.
 * Cluster ID and Hamiltonian ID are auto-derived (read-only). On completion, the
 * result is sent to {@code resultSink} which calls
 * {@link OutputPanel#showResult(ThermodynamicResult)}.</p>
 *
 * <p>Log output is routed to {@code logSink} and one-line status updates go to
 * {@code statusSink}.</p>
 */
public class CalculationPanel extends JPanel {

    // ── VS Code dark colours ──────────────────────────────────────────────────
    private static final Color BG       = new Color(0x252526);
    private static final Color LABEL_FG = new Color(0xCCCCCC);
    private static final Color ID_FG    = new Color(0x4EC9B0);   // teal

    private final CalculationService             service;
    private final WorkbenchContext               context;
    private final Consumer<String>               statusSink;
    private final Consumer<String>               logSink;
    private final Consumer<ThermodynamicResult>  resultSink;

    // User-editable inputs
    private final JTextField elementsField  = new JTextField("A-B", 12);
    private final JTextField structureField = new JTextField("BCC_B2", 10);
    private final JTextField modelField     = new JTextField("T", 8);

    // Auto-derived IDs (read-only)
    private final JTextField clusterIdField    = new JTextField(20);
    private final JTextField hamiltonianIdField = new JTextField(20);

    private final JSpinner          temperatureField = new JSpinner(new SpinnerNumberModel(1000.0, 1.0, 10000.0, 100.0));
    private final JSpinner          xBField          = new JSpinner(new SpinnerNumberModel(0.5, 0.0, 1.0, 0.05));
    private final JComboBox<String> engineBox        = new JComboBox<>(new String[]{"CVM", "MCS"});
    private final JButton           calcButton       = new JButton("Calculate");

    private boolean updatingFromContext = false;

    public CalculationPanel(CalculationService service,
                            WorkbenchContext context,
                            Consumer<String> statusSink,
                            Consumer<String> logSink,
                            Consumer<ThermodynamicResult> resultSink) {
        this.service     = service;
        this.context     = context;
        this.statusSink  = statusSink;
        this.logSink     = logSink;
        this.resultSink  = resultSink;

        setBackground(BG);

        for (JTextField f : new JTextField[]{ clusterIdField, hamiltonianIdField }) {
            f.setEditable(false);
            f.setFont(new Font(Font.MONOSPACED, Font.BOLD, 11));
            f.setForeground(ID_FG);
        }

        DocumentListener userEdit = new DocumentListener() {
            public void insertUpdate(DocumentEvent e)  { onFieldChanged(); }
            public void removeUpdate(DocumentEvent e)  { onFieldChanged(); }
            public void changedUpdate(DocumentEvent e) { onFieldChanged(); }
        };
        elementsField.getDocument().addDocumentListener(userEdit);
        structureField.getDocument().addDocumentListener(userEdit);
        modelField.getDocument().addDocumentListener(userEdit);

        context.addChangeListener(() -> {
            if (context.hasSystem()) {
                var sys = context.getSystem();
                updatingFromContext = true;
                if (!elementsField.getText().equals(sys.elements))   elementsField.setText(sys.elements);
                if (!structureField.getText().equals(sys.structure))  structureField.setText(sys.structure);
                if (!modelField.getText().equals(sys.model))          modelField.setText(sys.model);
                updatingFromContext = false;
                refreshIds();
            }
        });

        if (context.hasSystem()) {
            var sys = context.getSystem();
            elementsField.setText(sys.elements);
            structureField.setText(sys.structure);
            modelField.setText(sys.model);
        }

        refreshIds();

        setLayout(new BorderLayout(0, 0));
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));
        add(buildForm(), BorderLayout.NORTH);
    }

    // =========================================================================
    // ID derivation + context propagation
    // =========================================================================

    private void onFieldChanged() {
        refreshIds();
        if (!updatingFromContext) {
            String elements  = elementsField.getText().trim();
            String structure = structureField.getText().trim();
            String model     = modelField.getText().trim();
            if (!elements.isBlank() && !structure.isBlank() && !model.isBlank()) {
                context.setSystem(elements, structure, model);
            }
        }
    }

    private void refreshIds() {
        String elements  = elementsField.getText().trim();
        String structure = structureField.getText().trim();
        String model     = modelField.getText().trim();
        if (elements.isBlank() || structure.isBlank() || model.isBlank()) return;
        try {
            SystemId id = new SystemId(elements, structure, model);
            clusterIdField.setText(id.clusterId());
            hamiltonianIdField.setText(id.hamiltonianId());
        } catch (IllegalArgumentException ignored) {}
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
        lc.insets  = new Insets(2, 6, 1, 6);

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill    = GridBagConstraints.HORIZONTAL;
        fc.weightx = 1.0;
        fc.gridx   = 0;
        fc.insets  = new Insets(0, 6, 4, 8);

        int row = 0;

        // Elements
        lc.gridy = row++;  form.add(makeLabel("Elements:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(elementsField, fc);

        // Structure
        lc.gridy = row++;  form.add(makeLabel("Structure:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(structureField, fc);

        // Model
        lc.gridy = row++;  form.add(makeLabel("Model:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(modelField, fc);

        // Cluster ID
        lc.gridy = row++;  form.add(makeLabel("Cluster ID:", ID_FG), lc);
        fc.gridy = row++;  form.add(clusterIdField, fc);

        // Hamiltonian ID
        lc.gridy = row++;  form.add(makeLabel("Hamiltonian ID:", ID_FG), lc);
        fc.gridy = row++;  form.add(hamiltonianIdField, fc);

        // Temperature
        lc.gridy = row++;  form.add(makeLabel("Temperature (K):", LABEL_FG), lc);
        fc.gridy = row++;  form.add(temperatureField, fc);

        // Composition
        lc.gridy = row++;  form.add(makeLabel("Composition x_B:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(xBField, fc);

        // Engine
        lc.gridy = row++;  form.add(makeLabel("Engine:", LABEL_FG), lc);
        fc.gridy = row++;  form.add(engineBox, fc);

        // Button
        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 0; bc.gridy = row;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(6, 6, 4, 6);
        calcButton.addActionListener(e -> runCalculation());
        form.add(calcButton, bc);

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

    private void runCalculation() {
        String clusterId     = clusterIdField.getText().trim();
        String hamiltonianId = hamiltonianIdField.getText().trim();
        if (clusterId.isBlank() || hamiltonianId.isBlank()) {
            logSink.accept("Fill in Elements, Structure, and Model first.");
            return;
        }

        double temperature   = ((Number) temperatureField.getValue()).doubleValue();
        double xB            = ((Number) xBField.getValue()).doubleValue();
        String engineType    = (String) engineBox.getSelectedItem();
        double[] composition = {1.0 - xB, xB};

        calcButton.setEnabled(false);
        logSink.accept("Running " + engineType + " calculation...");
        logSink.accept("  Cluster     : " + clusterId);
        logSink.accept("  Hamiltonian : " + hamiltonianId);
        logSink.accept("  T = " + temperature + " K,  x_B = " + xB);
        statusSink.accept("Running " + engineType + " at T=" + temperature + " K...");

        SwingWorker<ThermodynamicResult, Void> worker = new SwingWorker<>() {
            @Override
            protected ThermodynamicResult doInBackground() throws Exception {
                return service.runSinglePoint(clusterId, hamiltonianId, temperature, composition, engineType);
            }

            @Override
            protected void done() {
                calcButton.setEnabled(true);
                try {
                    ThermodynamicResult result = get();
                    logSink.accept("Calculation complete.");
                    logSink.accept("  G = " + String.format("%.4f", result.gibbsEnergy) + " J/mol");
                    logSink.accept("  H = " + String.format("%.4f", result.enthalpy) + " J/mol");
                    logSink.accept("  S = " + String.format("%.6f", result.entropy) + " J/mol/K");
                    statusSink.accept("Done — G = " + String.format("%.4f", result.gibbsEnergy) + " J/mol");
                    resultSink.accept(result);
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };
        worker.execute();
    }
}
