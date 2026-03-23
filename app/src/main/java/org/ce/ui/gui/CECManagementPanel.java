package org.ce.ui.gui;

import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.workflow.cec.CECManagementWorkflow;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.util.function.Consumer;

/**
 * Parameter panel for scaffolding, inspecting, and editing the Hamiltonian (ECI) database.
 *
 * <p>Shown in the Explorer column. Log output is routed to {@code logSink}
 * (displayed in {@link OutputPanel}). One-line status updates go to
 * {@code statusSink} (displayed in {@link StatusBar}).</p>
 */
public class CECManagementPanel extends JPanel {

    // ── VS Code dark colours ──────────────────────────────────────────────────
    private static final Color BG        = new Color(0x252526);
    private static final Color LABEL_FG  = new Color(0xCCCCCC);
    private static final Color ID_FG     = new Color(0x4EC9B0);   // teal

    private final CECManagementWorkflow cecWorkflow;
    private final WorkbenchContext      context;
    private final Consumer<String>      statusSink;
    private final Consumer<String>      logSink;

    // User-editable inputs
    private final JTextField elementsField  = new JTextField("Nb-Ti", 12);
    private final JTextField structureField = new JTextField("BCC_A2", 10);
    private final JTextField modelField     = new JTextField("T", 8);

    // Auto-derived IDs (read-only)
    private final JTextField clusterIdField    = new JTextField(20);
    private final JTextField hamiltonianIdField = new JTextField(20);

    private final DefaultTableModel tableModel = new DefaultTableModel(
            new Object[]{"Name", "a (J/mol)", "b (J/mol/K)"}, 0) {
        @Override
        public boolean isCellEditable(int row, int column) {
            return column > 0;
        }
    };
    private final JTable cecTable = new JTable(tableModel);

    private CECEntry currentEntry         = null;
    private String   currentHamiltonianId = null;

    public CECManagementPanel(CECManagementWorkflow cecWorkflow,
                              WorkbenchContext context,
                              Consumer<String> statusSink,
                              Consumer<String> logSink) {
        this.cecWorkflow = cecWorkflow;
        this.context     = context;
        this.statusSink  = statusSink;
        this.logSink     = logSink;

        setBackground(BG);

        // Style read-only ID fields
        for (JTextField f : new JTextField[]{ clusterIdField, hamiltonianIdField }) {
            f.setEditable(false);
            f.setFont(new Font(Font.MONOSPACED, Font.BOLD, 11));
            f.setForeground(ID_FG);
        }

        // Style table for dark theme
        cecTable.setBackground(new Color(0x252526));
        cecTable.setForeground(new Color(0xD4D4D4));
        cecTable.setGridColor(new Color(0x3C3C3C));
        cecTable.setSelectionBackground(new Color(0x094771));
        cecTable.setSelectionForeground(Color.WHITE);
        cecTable.getTableHeader().setBackground(new Color(0x2D2D2D));
        cecTable.getTableHeader().setForeground(new Color(0xBBBBBB));

        // Wire document listeners
        DocumentListener refresh = new DocumentListener() {
            public void insertUpdate(DocumentEvent e)  { onFieldChanged(); }
            public void removeUpdate(DocumentEvent e)  { onFieldChanged(); }
            public void changedUpdate(DocumentEvent e) { onFieldChanged(); }
        };
        elementsField.getDocument().addDocumentListener(refresh);
        structureField.getDocument().addDocumentListener(refresh);
        modelField.getDocument().addDocumentListener(refresh);

        if (context.hasSystem()) {
            var sys = context.getSystem();
            elementsField.setText(sys.elements);
            structureField.setText(sys.structure);
            modelField.setText(sys.model);
        }

        context.addChangeListener(() -> {
            if (context.hasSystem()) {
                var sys = context.getSystem();
                if (!elementsField.getText().equals(sys.elements))  elementsField.setText(sys.elements);
                if (!structureField.getText().equals(sys.structure)) structureField.setText(sys.structure);
                if (!modelField.getText().equals(sys.model))        modelField.setText(sys.model);
            }
        });

        refreshIds();

        setLayout(new BorderLayout(6, 6));
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));
        add(buildForm(),         BorderLayout.NORTH);
        add(buildTableSection(), BorderLayout.CENTER);
    }

    // =========================================================================
    // ID derivation + context propagation
    // =========================================================================

    private void onFieldChanged() {
        refreshIds();
        String elements  = elementsField.getText().trim();
        String structure = structureField.getText().trim();
        String model     = modelField.getText().trim();
        if (!elements.isBlank() && !structure.isBlank() && !model.isBlank()) {
            context.setSystem(elements, structure, model);
        }
    }

    private void refreshIds() {
        String elements  = elementsField.getText().trim();
        String structure = structureField.getText().trim();
        String model     = modelField.getText().trim();
        if (elements.isBlank() || structure.isBlank() || model.isBlank()) return;
        try {
            String hamiltonianId = elements + "_" + structure + "_" + model;
            hamiltonianIdField.setText(hamiltonianId);
            // Cluster ID derived from hamiltonian ID
            int sep = hamiltonianId.indexOf('_');
            String structureModel = hamiltonianId.substring(sep + 1);
            int ncomp = elements.split("-").length;
            String suffix = switch (ncomp) {
                case 2 -> "bin";
                case 3 -> "tern";
                case 4 -> "quat";
                default -> "bin";
            };
            clusterIdField.setText(structureModel + "_" + suffix);
        } catch (Exception ignored) {}
    }

    // =========================================================================
    // Form layout
    // =========================================================================

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBackground(BG);
        form.setBorder(BorderFactory.createTitledBorder("Hamiltonian Database"));

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

        // Buttons
        JPanel btnPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        btnPanel.setBackground(BG);
        JButton scaffoldBtn = new JButton("Scaffold");
        JButton loadBtn     = new JButton("Load");
        JButton saveBtn     = new JButton("Save");
        JButton transformBtn = new JButton("Transform from Binary");
        scaffoldBtn.addActionListener(e -> scaffoldCEC());
        loadBtn.addActionListener(e -> loadCEC());
        saveBtn.addActionListener(e -> saveCEC());
        transformBtn.addActionListener(e -> transformFromBinary());
        btnPanel.add(scaffoldBtn);
        btnPanel.add(loadBtn);
        btnPanel.add(saveBtn);
        btnPanel.add(transformBtn);

        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 0; bc.gridy = row;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(6, 4, 4, 6);
        form.add(btnPanel, bc);

        return form;
    }

    private JScrollPane buildTableSection() {
        cecTable.setFillsViewportHeight(true);
        cecTable.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
        JScrollPane scroll = new JScrollPane(cecTable);
        scroll.setBorder(BorderFactory.createTitledBorder("ECI Terms  (edit a and b, then Save)"));
        scroll.setBackground(new Color(0x252526));
        scroll.getViewport().setBackground(new Color(0x252526));
        return scroll;
    }

    private static JLabel makeLabel(String text, Color fg) {
        JLabel l = new JLabel(text);
        l.setForeground(fg);
        return l;
    }

    // =========================================================================
    // Actions
    // =========================================================================

    private void scaffoldCEC() {
        String hamiltonianId = hamiltonianIdField.getText().trim();
        String elements      = elementsField.getText().trim();
        String structure     = structureField.getText().trim();
        String model         = modelField.getText().trim();
        if (hamiltonianId.isBlank()) return;

        logSink.accept("Scaffolding Hamiltonian for " + hamiltonianId + "...");
        statusSink.accept("Scaffolding " + hamiltonianId + "...");

        SwingWorker<CECEntry, Void> worker = new SwingWorker<>() {
            @Override
            protected CECEntry doInBackground() throws Exception {
                return cecWorkflow.scaffoldFromClusterData(hamiltonianId, elements, structure, model);
            }

            @Override
            protected void done() {
                try {
                    CECEntry entry = get();
                    currentEntry = entry;
                    currentHamiltonianId = hamiltonianId;
                    populateTable(entry);
                    logSink.accept("Scaffolded " + entry.ncf + " terms. Saved to hamiltonians/" + hamiltonianId + "/hamiltonian.json");
                    logSink.accept("Edit a and b values in the table, then click Save.");
                    statusSink.accept("Scaffolded " + hamiltonianId + "  (" + entry.ncf + " terms)");
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };
        worker.execute();
    }

    private void loadCEC() {
        String clusterId     = clusterIdField.getText().trim();
        String hamiltonianId = hamiltonianIdField.getText().trim();
        if (hamiltonianId.isBlank()) return;

        logSink.accept("Loading Hamiltonian " + hamiltonianId + "...");
        statusSink.accept("Loading " + hamiltonianId + "...");

        SwingWorker<CECEntry, Void> worker = new SwingWorker<>() {
            @Override
            protected CECEntry doInBackground() throws Exception {
                return cecWorkflow.loadAndValidateCEC(clusterId, hamiltonianId);
            }

            @Override
            protected void done() {
                try {
                    CECEntry entry = get();
                    currentEntry = entry;
                    currentHamiltonianId = hamiltonianId;
                    populateTable(entry);
                    logSink.accept("Loaded " + entry.ncf + " terms.");
                    logSink.accept("Elements: " + entry.elements
                            + "  Structure: " + entry.structurePhase
                            + "  Units: " + entry.cecUnits);
                    logSink.accept("Edit a and b values in the table, then click Save.");
                    statusSink.accept("Loaded " + hamiltonianId + "  (" + entry.ncf + " terms)");
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };
        worker.execute();
    }

    private void saveCEC() {
        if (currentEntry == null || currentHamiltonianId == null) {
            logSink.accept("Nothing loaded. Scaffold or Load a Hamiltonian first.");
            return;
        }

        if (cecTable.isEditing()) {
            cecTable.getCellEditor().stopCellEditing();
        }

        for (int row = 0; row < tableModel.getRowCount(); row++) {
            try {
                currentEntry.cecTerms[row].a = Double.parseDouble(tableModel.getValueAt(row, 1).toString());
                currentEntry.cecTerms[row].b = Double.parseDouble(tableModel.getValueAt(row, 2).toString());
            } catch (NumberFormatException ex) {
                logSink.accept("Invalid number at row " + (row + 1) + ": " + ex.getMessage());
                return;
            }
        }

        String id    = currentHamiltonianId;
        CECEntry ent = currentEntry;
        logSink.accept("Saving Hamiltonian " + id + "...");
        statusSink.accept("Saving " + id + "...");

        SwingWorker<Void, Void> worker = new SwingWorker<>() {
            @Override
            protected Void doInBackground() throws Exception {
                cecWorkflow.store.save(id, ent);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                    logSink.accept("Saved to hamiltonians/" + id + "/hamiltonian.json");
                    statusSink.accept("Saved " + id);
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };
        worker.execute();
    }

    private void populateTable(CECEntry entry) {
        tableModel.setRowCount(0);
        if (entry.cecTerms == null) return;
        for (CECTerm term : entry.cecTerms) {
            tableModel.addRow(new Object[]{ term.name, term.a, term.b });
        }
    }

    // =========================================================================
    // Basis Transformation: Binary → Ternary
    // =========================================================================

    private void transformFromBinary() {
        String ternaryId = hamiltonianIdField.getText().trim();
        if (ternaryId.isBlank()) {
            logSink.accept("Please set target Hamiltonian ID first");
            statusSink.accept("Error: target Hamiltonian ID required");
            return;
        }

        String ternaryElements = elementsField.getText().trim();
        if (ternaryElements.isBlank() || ternaryElements.split("-").length != 3) {
            logSink.accept("Please set Ternary elements (e.g., Nb-Ti-V)");
            statusSink.accept("Error: invalid ternary element string");
            return;
        }

        // Show input dialog for binary Hamiltonian ID
        String binaryId = JOptionPane.showInputDialog(
            this,
            "Enter source binary Hamiltonian ID\n(e.g., Nb-Ti_BCC_A2_T):",
            "Nb-Ti_BCC_A2_T"
        );

        if (binaryId == null || binaryId.isBlank()) return;

        logSink.accept("Transforming " + binaryId + " → " + ternaryId + "...");
        statusSink.accept("Transforming " + binaryId + " to ternary basis...");

        SwingWorker<CECEntry, Void> worker = new SwingWorker<>() {
            @Override
            protected CECEntry doInBackground() throws Exception {
                // Build species mapping from element strings
                String[] binaryElems = binaryId.split("_")[0].split("-");
                String[] ternaryElems = ternaryElements.split("-");

                if (binaryElems.length != 2) {
                    throw new IllegalArgumentException(
                        "Binary ID must have 2 elements, got: " + binaryId);
                }

                if (ternaryElems.length != 3) {
                    throw new IllegalArgumentException(
                        "Ternary elements must have 3, got: " + ternaryElements);
                }

                // Create species mapping
                org.ce.domain.hamiltonian.NumericalCECTransformer.SpeciesMapping mapping =
                    new org.ce.domain.hamiltonian.NumericalCECTransformer.SpeciesMapping();

                // Map binary species to ternary species by name
                for (int b = 0; b < binaryElems.length; b++) {
                    for (int t = 0; t < ternaryElems.length; t++) {
                        if (binaryElems[b].equalsIgnoreCase(ternaryElems[t])) {
                            mapping.addMapping(t, b);
                        }
                    }
                }

                // Perform transformation
                return cecWorkflow.transformBinaryToTernary(
                    binaryId,
                    ternaryId,
                    ternaryElements,
                    mapping
                );
            }

            @Override
            protected void done() {
                try {
                    CECEntry entry = get();
                    currentEntry = entry;
                    currentHamiltonianId = ternaryId;
                    populateTable(entry);
                    logSink.accept("Transformed " + entry.ncf + " terms to ternary basis");
                    logSink.accept("Saved to hamiltonians/" + ternaryId + "/hamiltonian.json");
                    logSink.accept("Review the transformed values and click Save if satisfied.");
                    statusSink.accept("Transformed " + ternaryId);
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };
        worker.execute();
    }
}
