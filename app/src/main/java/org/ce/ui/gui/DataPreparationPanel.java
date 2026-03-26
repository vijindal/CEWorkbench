package org.ce.ui.gui;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrixBuilder;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.Vector3D;
import org.ce.domain.cluster.cvcf.BccA2CvCfTransformations;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.InputLoader;
import org.ce.storage.Workspace;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;

/**
 * Parameter panel for Type-1 cluster identification (shown in the Explorer column).
 *
 * <p>The user selects two sets of files:</p>
 * <ol>
 *   <li>The <em>ordered</em> cluster file + symmetry group (e.g. BCC_B2-T + BCC_B2-SG).</li>
 *   <li>The <em>disordered</em> (parent/HSP) cluster file + symmetry group (e.g. BCC_A2-T).</li>
 * </ol>
 *
 * <p>Log output is routed to {@code logSink} (displayed in {@link OutputPanel}).
 * One-line status updates go to {@code statusSink} (displayed in {@link StatusBar}).</p>
 */
public class DataPreparationPanel extends JPanel {

    // ── VS Code dark colours used in this panel ───────────────────────────────
    private static final Color BG         = new Color(0x252526);
    private static final Color LABEL_FG   = new Color(0xCCCCCC);
    private static final Color ORD_HDR    = new Color(0x569CD6);   // VS Code blue
    private static final Color DIS_HDR    = new Color(0xCE9178);   // VS Code string orange
    private static final Color SYS_ID_FG  = new Color(0x4EC9B0);   // teal

    private final ClusterDataStore  clusterStore;
    private final Path              inputsDir;
    private final WorkbenchContext  context;
    private final Consumer<String>  statusSink;
    private final Consumer<String>  logSink;

    // Ordered phase (target) — determines system ID
    private final JComboBox<String> orderedClusterCombo;
    private final JComboBox<String> orderedSymCombo;

    // Disordered phase (parent/HSP) — determines ncf for Hamiltonian scaffold
    private final JComboBox<String> disorderedClusterCombo;
    private final JComboBox<String> disorderedSymCombo;

    private final JSpinner   numCompSpinner = new JSpinner(new SpinnerNumberModel(2, 2, 20, 1));
    private final JTextField systemIdField  = new JTextField(24);

    public DataPreparationPanel(ClusterDataStore clusterStore,
                                WorkbenchContext context,
                                Consumer<String> statusSink,
                                Consumer<String> logSink) {
        this.clusterStore = clusterStore;
        this.inputsDir    = new Workspace().inputsDir();
        this.context      = context;
        this.statusSink   = statusSink;
        this.logSink      = logSink;

        setBackground(BG);

        String[] clusFiles = scanInputsDir(inputsDir, "clus", ".txt", "clus/", false);
        String[] symGroups = scanInputsDir(inputsDir, "sym",  ".txt", "", true);

        orderedClusterCombo    = makeCombo(clusFiles);
        orderedSymCombo        = makeCombo(symGroups);
        disorderedClusterCombo = makeCombo(clusFiles);
        disorderedSymCombo     = makeCombo(symGroups);

        systemIdField.setEditable(false);
        systemIdField.setFont(new Font(Font.MONOSPACED, Font.BOLD, 12));
        systemIdField.setForeground(SYS_ID_FG);

        orderedClusterCombo.addActionListener(e -> refreshSystemId());
        numCompSpinner.addChangeListener(e -> refreshSystemId());
        refreshSystemId();

        setLayout(new BorderLayout(0, 0));
        setBorder(BorderFactory.createEmptyBorder(8, 10, 10, 10));
        add(buildForm(), BorderLayout.NORTH);
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
                    if (stripExt) name = name.substring(0, name.length() - extension.length());
                    items[i] = prefix + name;
                }
                return items;
            }
        } catch (Exception ignored) {}
        return new String[0];
    }

    // =========================================================================
    // System-ID generation
    // =========================================================================

    private void refreshSystemId() {
        Object sel = orderedClusterCombo.getSelectedItem();
        if (sel == null || sel.toString().isBlank()) return;
        systemIdField.setText(generateSystemId(sel.toString(), (int) numCompSpinner.getValue()));
    }

    private static String generateSystemId(String clusterFile, int numComp) {
        String base = clusterFile;
        int slash = base.lastIndexOf('/');
        if (slash >= 0) base = base.substring(slash + 1);
        if (base.toLowerCase().endsWith(".txt")) base = base.substring(0, base.length() - 4);
        base = base.replace('-', '_');

        String suffix = switch (numComp) {
            case 2  -> "bin";
            case 3  -> "tern";
            case 4  -> "quat";
            default -> numComp + "comp";
        };
        return base + "_" + suffix;
    }

    private static String[] parseStructureModel(String systemId) {
        String stripped = systemId.replaceAll("_(bin|tern|quat|\\d+comp)$", "");
        int last = stripped.lastIndexOf('_');
        if (last < 0) return null;
        return new String[]{ stripped.substring(0, last), stripped.substring(last + 1) };
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
        hc.gridx = 0; hc.gridwidth = 2; hc.fill = GridBagConstraints.HORIZONTAL;
        hc.anchor = GridBagConstraints.WEST;
        hc.insets = new Insets(8, 6, 2, 6);

        // Label (full width)
        GridBagConstraints lc = new GridBagConstraints();
        lc.gridx = 0; lc.gridwidth = 2; lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(2, 6, 1, 6);

        // Combo (expands)
        GridBagConstraints fc = new GridBagConstraints();
        fc.gridx = 0; fc.fill = GridBagConstraints.HORIZONTAL; fc.weightx = 1.0;
        fc.insets = new Insets(0, 6, 3, 2);

        // Browse button (fixed)
        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 1; bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(0, 2, 3, 6);

        int row = 0;

        // ── Ordered phase ──
        JLabel ordLabel = new JLabel("── Ordered phase (target) ──");
        ordLabel.setForeground(ORD_HDR);
        ordLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 11));
        hc.gridy = row++;
        form.add(ordLabel, hc);

        lc.gridy = row++;  form.add(makeLabel("Cluster file:", LABEL_FG), lc);
        fc.gridy = row; bc.gridy = row++;
        form.add(orderedClusterCombo, fc);
        JButton browseOrdClus = new JButton("Browse");
        browseOrdClus.addActionListener(e -> browseFile(orderedClusterCombo, "clus", "clus/", false));
        form.add(browseOrdClus, bc);

        lc.gridy = row++;  form.add(makeLabel("Symmetry group:", LABEL_FG), lc);
        fc.gridy = row; bc.gridy = row++;
        form.add(orderedSymCombo, fc);
        JButton browseOrdSym = new JButton("Browse");
        browseOrdSym.addActionListener(e -> browseFile(orderedSymCombo, "sym", "", true));
        form.add(browseOrdSym, bc);

        // ── Disordered (parent/HSP) phase ──
        JLabel disLabel = new JLabel("── Disordered parent (HSP) ──");
        disLabel.setForeground(DIS_HDR);
        disLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 11));
        hc.gridy = row++;
        form.add(disLabel, hc);

        lc.gridy = row++;  form.add(makeLabel("Cluster file:", LABEL_FG), lc);
        fc.gridy = row; bc.gridy = row++;
        form.add(disorderedClusterCombo, fc);
        JButton browseDisClus = new JButton("Browse");
        browseDisClus.addActionListener(e -> browseFile(disorderedClusterCombo, "clus", "clus/", false));
        form.add(browseDisClus, bc);

        lc.gridy = row++;  form.add(makeLabel("Symmetry group:", LABEL_FG), lc);
        fc.gridy = row; bc.gridy = row++;
        form.add(disorderedSymCombo, fc);
        JButton browseDisSym = new JButton("Browse");
        browseDisSym.addActionListener(e -> browseFile(disorderedSymCombo, "sym", "", true));
        form.add(browseDisSym, bc);

        // ── Shared ──
        lc.gridy = row++;  form.add(makeLabel("Num components:", LABEL_FG), lc);
        fc.gridy = row; bc.gridy = row++;
        form.add(numCompSpinner, fc);

        lc.gridy = row++;  form.add(makeLabel("System ID (auto):", SYS_ID_FG), lc);
        fc.gridy = row++;
        form.add(systemIdField, fc);

        // Run button
        GridBagConstraints rbc = new GridBagConstraints();
        rbc.gridx = 0; rbc.gridy = row; rbc.gridwidth = 2;
        rbc.anchor = GridBagConstraints.WEST;
        rbc.insets = new Insets(10, 6, 4, 6);
        JButton runBtn = new JButton("Run Identification & Save");
        runBtn.addActionListener(e -> runIdentification());
        form.add(runBtn, rbc);

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
            if (stripExt && name.endsWith(".txt")) name = name.substring(0, name.length() - 4);
            String item = prefix + name;
            addIfAbsent(combo, item);
            combo.setSelectedItem(item);
        }
    }

    private static void addIfAbsent(JComboBox<String> combo, String item) {
        for (int i = 0; i < combo.getItemCount(); i++) {
            if (item.equals(combo.getItemAt(i))) return;
        }
        combo.addItem(item);
    }

    // =========================================================================
    // Identification workflow
    // =========================================================================

    private void runIdentification() {
        Object ordClusSel = orderedClusterCombo.getSelectedItem();
        Object ordSymSel  = orderedSymCombo.getSelectedItem();
        Object disClusSel = disorderedClusterCombo.getSelectedItem();
        Object disSymSel  = disorderedSymCombo.getSelectedItem();

        if (ordClusSel == null || ordSymSel == null
                || disClusSel == null || disSymSel == null) return;

        String ordClus  = ordClusSel.toString().trim();
        String ordSym   = ordSymSel.toString().trim();
        String disClus  = disClusSel.toString().trim();
        String disSym   = disSymSel.toString().trim();
        int    numComp  = (int) numCompSpinner.getValue();
        String systemId = systemIdField.getText().trim();

        logSink.accept("System ID      : " + systemId);
        logSink.accept("Ordered cluster: " + ordClus);
        logSink.accept("Ordered sym    : " + ordSym);
        logSink.accept("Disordered clus: " + disClus + "  (parent — determines ncf)");
        logSink.accept("Disordered sym : " + disSym);
        logSink.accept("Components     : " + numComp);
        statusSink.accept("Running cluster identification for " + systemId + "...");

        SwingWorker<AllClusterData, String> worker = new SwingWorker<>() {
            @Override
            protected AllClusterData doInBackground() throws Exception {
                publish("Stage 1-2: Cluster + CF identification...");

                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .orderedClusterFile(ordClus)
                        .orderedSymmetryGroup(ordSym)
                        .disorderedClusterFile(disClus)
                        .disorderedSymmetryGroup(disSym)
                        .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})
                        .translationVector(new Vector3D(0,0,0))
                        .numComponents(numComp)
                        .build();

                AllClusterData partial = ClusterIdentificationWorkflow.identify(config);

                publish("Stage 3: Building C-matrix...");
                List<Cluster> maxClusters = InputLoader.parseClusterFile(ordClus);
                CMatrixResult cmatrix = CMatrixBuilder.build(
                        partial.getDisorderedClusterResult(),
                        partial.getDisorderedCFResult(),
                        maxClusters,
                        numComp,
                        BccA2CvCfTransformations.binaryBasis()
                );

                publish("Saving AllClusterData...");
                AllClusterData full = new AllClusterData(
                        partial.getDisorderedClusterResult(),
                        partial.getOrderedClusterResult(),
                        partial.getDisorderedCFResult(),
                        partial.getOrderedCFResult(),
                        cmatrix
                );
                clusterStore.save(systemId, full);
                return full;
            }

            @Override
            protected void process(List<String> chunks) {
                for (String msg : chunks) logSink.accept(msg);
            }

            @Override
            protected void done() {
                try {
                    AllClusterData result = get();
                    logSink.accept("Done!  " + result.getSummary());
                    logSink.accept("Saved to ~/CEWorkbench/cluster-data/" + systemId + "/cluster_data.json");
                    logSink.accept("  disorderedCFResult ncf = "
                            + result.getDisorderedCFResult().getNcf() + "  (from " + disClus + ")");
                    logSink.accept("  orderedCFResult    ncf = "
                            + result.getOrderedCFResult().getNcf()    + "  (from " + ordClus + ")");
                    statusSink.accept("Identification complete — saved " + systemId);

                    String[] sm = parseStructureModel(systemId);
                    if (sm != null) {
                        String elements = context.hasSystem() ? context.getSystem().elements : "";
                        if (!elements.isBlank()) {
                            context.setSystem(elements, sm[0], sm[1]);
                        }
                    }
                } catch (Exception ex) {
                    logSink.accept("Error: " + ex.getMessage());
                    statusSink.accept("Error: " + ex.getMessage());
                }
            }
        };

        worker.execute();
    }
}
