package org.ce.ui.gui;

import org.ce.model.cluster.ClusterCFIdentificationPipeline.PipelineResult;
import org.ce.model.cluster.ClusterIdentificationRequest;
import org.ce.model.cluster.StructurePhaseRegistry;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * Parameter panel for Type-1 cluster identification (shown in the Explorer
 * column).
 *
 * <p>
 * Log output is routed via the unified {@link org.ce.CEWorkbenchContext#log(String)} API.
 * One-line status updates go to {@code statusSink} (displayed in
 * {@link StatusBar}).
 * </p>
 */
public class DataPreparationPanel extends JPanel {

    // ── VS Code dark colours used in this panel ───────────────────────────────
    private static final Color BG = new Color(0x252526);
    private static final Color LABEL_FG = new Color(0xCCCCCC);
    private static final Color ORD_HDR = new Color(0x569CD6); // VS Code blue
    private static final Color DIS_HDR = new Color(0xCE9178); // VS Code string orange
    private static final Color SYS_ID_FG = new Color(0x4EC9B0); // teal

    private final org.ce.CEWorkbenchContext appCtx;
    private final WorkbenchContext context;
    private final Path inputsDir;
    private final Consumer<String> statusSink;

    // Ordered phase (target) — the only file the user selects; everything
    // else (ordered symmetry, disordered structure+cluster, disordered
    // symmetry) is derived from this via StructurePhaseRegistry.
    private final JComboBox<String> orderedClusterCombo;

    // Derived — read-only, shown for transparency, not user-editable.
    private final JTextField orderedSymField = new JTextField();
    private final JTextField disorderedClusterField = new JTextField();
    private final JTextField disorderedSymField = new JTextField();
    private final JLabel derivationWarningLabel = new JLabel(" ");

    private final JTextField elementsField  = new JTextField("Nb-Ti", 12);
    private final JSpinner numCompSpinner = new JSpinner(new SpinnerNumberModel(2, 2, 20, 1));
    private final JTextField systemIdField = new JTextField(24);

    private final JButton runBtn;

    /** Holds the fully-resolved derivation, or null if resolution failed. */
    private ClusterIdentificationRequest derived;

    public DataPreparationPanel(org.ce.CEWorkbenchContext appCtx,
            WorkbenchContext context,
            Consumer<String> statusSink) {
        this.appCtx = appCtx;
        this.context = context;
        this.inputsDir = new org.ce.model.storage.Workspace().inputsDir();
        this.statusSink = statusSink;

        setBackground(BG);

        String[] clusFiles = scanInputsDir(inputsDir, "clus", ".txt", "", false);

        orderedClusterCombo = makeCombo(clusFiles);

        for (JTextField f : new JTextField[] { orderedSymField, disorderedClusterField, disorderedSymField }) {
            f.setEditable(false);
            f.setForeground(LABEL_FG);
        }
        derivationWarningLabel.setForeground(new Color(0xF48771)); // VS Code error orange

        systemIdField.setEditable(false);
        systemIdField.setFont(new Font(Font.MONOSPACED, Font.BOLD, 12));
        systemIdField.setForeground(SYS_ID_FG);

        orderedClusterCombo.addActionListener(e -> refreshDerivedFields());
        numCompSpinner.addChangeListener(e -> refreshSystemId());
        // Keep numComp in sync with elements field
        elementsField.getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            public void insertUpdate(javax.swing.event.DocumentEvent e)  { syncNumComp(); }
            public void removeUpdate(javax.swing.event.DocumentEvent e)  { syncNumComp(); }
            public void changedUpdate(javax.swing.event.DocumentEvent e) { syncNumComp(); }
            private void syncNumComp() {
                int n = elementsField.getText().trim().split("-").length;
                if (n >= 2) numCompSpinner.setValue(n);
            }
        });
        refreshDerivedFields();

        runBtn = new JButton("Run Identification");

        setLayout(new BorderLayout(0, 0));
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));
        add(buildForm(), BorderLayout.NORTH);
    }

    // =========================================================================
    // Derivation: ordered cluster file -> {ordered sym, disordered structure,
    // disordered cluster, disordered sym}
    //
    // Delegates entirely to ClusterIdentificationRequest.Builder — the same
    // model-layer code path the CLI uses — so the GUI's preview and the
    // actual pipeline dispatch can never diverge.
    // =========================================================================

    /**
     * Attempts to build a {@link ClusterIdentificationRequest} from a single
     * ordered cluster-file selection, using {@link ClusterIdentificationRequest.Builder}'s
     * derivation of the disordered parent (structure, model, cluster file,
     * symmetry group) via {@link StructurePhaseRegistry}.
     *
     * @param orderedClusterSelection bare filename as shown in the dropdown
     *        (e.g. {@code "BCC_B2-T.txt"}), not the {@code "clus/"}-prefixed
     *        path {@link ClusterIdentificationRequest} expects.
     * @return the built request, or null if derivation/validation failed
     *         (a message is set on {@link #derivationWarningLabel} in that case).
     */
    private ClusterIdentificationRequest deriveFromOrderedCluster(String orderedClusterSelection) {
        derivationWarningLabel.setText(" ");
        try {
            return ClusterIdentificationRequest.builder()
                    .orderedClusterFile("clus/" + orderedClusterSelection)
                    .numComponents((int) numCompSpinner.getValue())
                    .build();
        } catch (Exception e) {
            derivationWarningLabel.setText(e.getMessage());
            return null;
        }
    }

    private void refreshDerivedFields() {
        Object sel = orderedClusterCombo.getSelectedItem();
        if (sel == null || sel.toString().isBlank()) {
            derived = null;
            orderedSymField.setText("");
            disorderedClusterField.setText("");
            disorderedSymField.setText("");
            return;
        }
        derived = deriveFromOrderedCluster(sel.toString().trim());
        if (derived != null) {
            orderedSymField.setText(derived.getOrderedSymmetryGroup());
            disorderedClusterField.setText(derived.getDisorderedClusterFile());
            disorderedSymField.setText(derived.getDisorderedSymmetryGroup());
        } else {
            orderedSymField.setText("");
            disorderedClusterField.setText("");
            disorderedSymField.setText("");
        }
        refreshSystemId();
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private static JComboBox<String> makeCombo(String[] items) {
        JComboBox<String> cb = new JComboBox<>(items);
        cb.setEditable(true);
        return cb;
    }

    private static String[] scanInputsDir(Path inputsDir, String subDir, String extension,
            String prefix, boolean stripExt) {
        try {
            File dir = inputsDir.resolve(subDir).toFile();
            File[] files = dir.listFiles((d, n) -> n.endsWith(extension));
            if (files != null && files.length > 0) {
                Arrays.sort(files);
                String[] items = new String[files.length];
                for (int i = 0; i < files.length; i++) {
                    String name = files[i].getName();
                    if (stripExt)
                        name = name.substring(0, name.length() - extension.length());
                    items[i] = prefix + name;
                }
                return items;
            }
        } catch (Exception ignored) {
        }
        return new String[0];
    }

    // =========================================================================
    // System-ID generation
    // =========================================================================

    private void refreshSystemId() {
        if (derived == null) {
            systemIdField.setText("");
            return;
        }
        systemIdField.setText(derived.getStructurePhase() + "_" + derived.getModel());
    }

    // =========================================================================
    // Form layout
    // =========================================================================

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBackground(BG);
        form.setBorder(BorderFactory.createTitledBorder("Cluster Identification"));

        // Full-width header constraint
        GridBagConstraints hc = new GridBagConstraints();
        hc.gridx = 0;
        hc.gridwidth = 2;
        hc.fill = GridBagConstraints.HORIZONTAL;
        hc.anchor = GridBagConstraints.WEST;
        hc.insets = new Insets(8, 6, 2, 6);

        // Label (full width)
        GridBagConstraints lc = new GridBagConstraints();
        lc.gridx = 0;
        lc.gridwidth = 2;
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(2, 6, 1, 6);

        // Combo (expands)
        GridBagConstraints fc = new GridBagConstraints();
        fc.gridx = 0;
        fc.fill = GridBagConstraints.HORIZONTAL;
        fc.weightx = 1.0;
        fc.insets = new Insets(0, 6, 3, 2);

        // Browse button (fixed)
        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 1;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(0, 2, 3, 6);

        int row = 0;

        // ── Ordered phase (the only user selection) ──
        JLabel ordLabel = new JLabel("── Ordered phase (target) ──");
        ordLabel.setForeground(ORD_HDR);
        ordLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 11));
        hc.gridy = row++;
        form.add(ordLabel, hc);

        lc.gridy = row++;
        form.add(makeLabel("Cluster file:", LABEL_FG), lc);
        fc.gridy = row;
        bc.gridy = row++;
        form.add(orderedClusterCombo, fc);
        JButton browseOrdClus = new JButton("Browse");
        browseOrdClus.addActionListener(e -> browseFile(orderedClusterCombo, "clus", "", false));
        form.add(browseOrdClus, bc);

        // ── Derived (read-only) ──
        JLabel derivedLabel = new JLabel("── Derived automatically ──");
        derivedLabel.setForeground(DIS_HDR);
        derivedLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 11));
        hc.gridy = row++;
        form.add(derivedLabel, hc);

        lc.gridy = row++;
        form.add(makeLabel("Ordered symmetry group:", LABEL_FG), lc);
        fc.gridy = row++;
        form.add(orderedSymField, fc);

        lc.gridy = row++;
        form.add(makeLabel("Disordered parent cluster file:", LABEL_FG), lc);
        fc.gridy = row++;
        form.add(disorderedClusterField, fc);

        lc.gridy = row++;
        form.add(makeLabel("Disordered parent symmetry group:", LABEL_FG), lc);
        fc.gridy = row++;
        form.add(disorderedSymField, fc);

        hc.gridy = row++;
        form.add(derivationWarningLabel, hc);

        // ── Shared ──
        lc.gridy = row++;
        form.add(makeLabel("Elements (e.g. Nb-Ti):", LABEL_FG), lc);
        fc.gridy = row++;
        form.add(elementsField, fc);

        lc.gridy = row++;
        form.add(makeLabel("Num components:", LABEL_FG), lc);
        fc.gridy = row;
        bc.gridy = row++;
        form.add(numCompSpinner, fc);

        lc.gridy = row++;
        form.add(makeLabel("System ID (auto):", SYS_ID_FG), lc);
        fc.gridy = row++;
        form.add(systemIdField, fc);

        // Buttons row
        JPanel btnRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
        btnRow.setBackground(BG);
        runBtn.addActionListener(e -> runIdentification());
        btnRow.add(runBtn);

        GridBagConstraints rbc = new GridBagConstraints();
        rbc.gridx = 0;
        rbc.gridy = row;
        rbc.gridwidth = 2;
        rbc.anchor = GridBagConstraints.WEST;
        rbc.insets = new Insets(10, 6, 4, 6);
        form.add(btnRow, rbc);

        return form;
    }

    private static JLabel makeLabel(String text, Color fg) {
        JLabel l = new JLabel(text);
        l.setForeground(fg);
        return l;
    }

    // =========================================================================
    // Browse helpers
    // =========================================================================

    private void browseFile(JComboBox<String> combo, String subDir, String prefix, boolean stripExt) {
        JFileChooser chooser = new JFileChooser(inputsDir.resolve(subDir).toFile());
        chooser.setFileFilter(new FileNameExtensionFilter("Text files (*.txt)", "txt"));
        chooser.setDialogTitle("Select file");
        if (chooser.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
            String name = chooser.getSelectedFile().getName();
            if (stripExt && name.endsWith(".txt"))
                name = name.substring(0, name.length() - 4);
            String item = prefix + name;
            addIfAbsent(combo, item);
            combo.setSelectedItem(item);
        }
    }

    private static void addIfAbsent(JComboBox<String> combo, String item) {
        for (int i = 0; i < combo.getItemCount(); i++) {
            if (item.equals(combo.getItemAt(i)))
                return;
        }
        combo.addItem(item);
    }

    // =========================================================================
    // Identification workflow
    // =========================================================================

    private void runIdentification() {
        Object ordClusSel = orderedClusterCombo.getSelectedItem();
        if (ordClusSel == null || ordClusSel.toString().isBlank())
            return;

        // Re-derive (rather than trusting a possibly-stale `derived` field)
        // so the user always gets an up-to-date check right before running.
        ClusterIdentificationRequest config = deriveFromOrderedCluster(ordClusSel.toString().trim());
        if (config == null) {
            statusSink.accept("Cannot run: " + derivationWarningLabel.getText());
            appCtx.log("Error: " + derivationWarningLabel.getText());
            return;
        }

        String systemId = systemIdField.getText().trim();

        appCtx.log("System ID      : " + systemId);
        appCtx.log("Ordered cluster: " + config.getOrderedClusterFile());
        appCtx.log("Ordered sym    : " + config.getOrderedSymmetryGroup());
        appCtx.log("Disordered clus: " + config.getDisorderedClusterFile());
        appCtx.log("Disordered sym : " + config.getDisorderedSymmetryGroup());
        appCtx.log("Components     : " + config.getNumComponents());
        appCtx.clearLog();
        statusSink.accept("Running cluster identification for " + systemId + "...");

        final String resolvedStructure = config.getStructurePhase();
        final String resolvedModel = config.getModel();

        runBtn.setEnabled(false);

        SwingWorker<PipelineResult, String> worker = new SwingWorker<PipelineResult, String>() {
            @Override
            protected PipelineResult doInBackground() throws Exception {
                publish("Stage 1-2: Cluster + CF identification...");

                return org.ce.model.cluster.ClusterCFIdentificationPipeline.runFullWorkflow(config, this::publish);
            }

            @Override
            protected void process(List<String> chunks) {
                for (String msg : chunks)
                    appCtx.log(msg);
            }

            @Override
            protected void done() {
                runBtn.setEnabled(true);
                try {
                    get();
                    appCtx.log("\nIdentification complete.");
                    appCtx.log("Click [Rebuild] in the Session Bar to build a session for this system.");

                    // Update system identity — SessionBar picks up the change and syncs its combos.
                    String elements = elementsField.getText().trim();
                    if (!elements.isBlank()) {
                        context.setSystem(elements, resolvedStructure, resolvedModel);
                    }
                    statusSink.accept("Identification done — click Rebuild in Session Bar to continue.");
                } catch (Exception ex) {
                    appCtx.log("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };

        worker.execute();
    }
}
