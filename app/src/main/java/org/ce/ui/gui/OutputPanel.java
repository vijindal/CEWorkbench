package org.ce.ui.gui;

import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.SystemId;

import javax.swing.*;
import java.awt.*;

/**
 * VS Code-style main output panel (third column, fills remaining width).
 *
 * <p>Divided vertically by a draggable splitter:</p>
 * <ul>
 *   <li><b>Top — RESULTS</b>: displays the latest {@link ThermodynamicResult} with
 *       a dark card layout. Subscribes to {@link WorkbenchContext} to show the
 *       active system identity.</li>
 *   <li><b>Bottom — OUTPUT</b>: dark terminal-style {@link JTextArea} where all
 *       panels append log lines via {@link #appendLog(String)}.</li>
 * </ul>
 *
 * <p>Visual spec (VS Code Dark+):</p>
 * <ul>
 *   <li>Background:       {@code #1E1E1E}</li>
 *   <li>Section headers:  {@code #2D2D2D} bar, {@code #BBBBBB} text</li>
 *   <li>Label dim:        {@code #858585}</li>
 *   <li>Value text:       {@code #D4D4D4} monospaced</li>
 *   <li>Teal accent:      {@code #4EC9B0} (system identity)</li>
 *   <li>Log bg/fg:        {@code #1E1E1E} / {@code #CCCCCC}</li>
 * </ul>
 */
public class OutputPanel extends JPanel {

    // ── colours ───────────────────────────────────────────────────────────────
    private static final Color BG           = new Color(0x1E1E1E);
    private static final Color HDR_BG       = new Color(0x2D2D2D);
    private static final Color HDR_FG       = new Color(0xBBBBBB);
    private static final Color BORDER_CLR   = new Color(0x3C3C3C);
    private static final Color LABEL_DIM    = new Color(0x858585);
    private static final Color VALUE_FG     = new Color(0xD4D4D4);
    private static final Color TEAL         = new Color(0x4EC9B0);
    private static final Color LOG_FG       = new Color(0xCCCCCC);

    // ── result labels ─────────────────────────────────────────────────────────
    private final JLabel sysLabel   = makeValueLabel("— no system selected —");
    private final JLabel tempLabel  = makeValueLabel("—");
    private final JLabel compLabel  = makeValueLabel("—");
    private final JLabel gibbsLabel = makeValueLabel("—");
    private final JLabel enthlLabel = makeValueLabel("—");
    private final JLabel entrLabel  = makeValueLabel("—");

    // ── log ───────────────────────────────────────────────────────────────────
    private final JTextArea logArea = new JTextArea();

    public OutputPanel(WorkbenchContext context) {
        setLayout(new BorderLayout());
        setBackground(BG);

        JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                buildResultsSection(context),
                buildLogSection());
        split.setDividerLocation(300);
        split.setDividerSize(4);
        split.setContinuousLayout(true);
        split.setBorder(null);
        split.setBackground(BG);

        add(split, BorderLayout.CENTER);

        // Update system label when context changes
        context.addChangeListener(() -> {
            if (context.hasSystem()) {
                SystemId s = context.getSystem();
                sysLabel.setText(s.elements + "  /  " + s.structure + "  /  " + s.model);
                sysLabel.setForeground(TEAL);
            } else {
                sysLabel.setText("— no system selected —");
                sysLabel.setForeground(LABEL_DIM);
            }
        });
    }

    // =========================================================================
    // Results section (top)
    // =========================================================================

    private JPanel buildResultsSection(WorkbenchContext context) {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(BG);
        panel.add(makeSectionHeader("RESULTS"), BorderLayout.NORTH);

        JPanel grid = new JPanel(new GridBagLayout());
        grid.setBackground(BG);
        grid.setBorder(BorderFactory.createEmptyBorder(16, 20, 16, 20));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(6, 0, 6, 20);

        GridBagConstraints vc = new GridBagConstraints();
        vc.anchor  = GridBagConstraints.WEST;
        vc.fill    = GridBagConstraints.HORIZONTAL;
        vc.weightx = 1.0;
        vc.insets  = new Insets(6, 0, 6, 0);

        sysLabel.setForeground(LABEL_DIM);

        addRow(grid, lc, vc, 0, "System",            sysLabel);
        addRow(grid, lc, vc, 1, "Temperature  (K)",  tempLabel);
        addRow(grid, lc, vc, 2, "Composition  x_B",  compLabel);
        addRow(grid, lc, vc, 3, "Gibbs  (J/mol)",    gibbsLabel);
        addRow(grid, lc, vc, 4, "Enthalpy  (J/mol)", enthlLabel);
        addRow(grid, lc, vc, 5, "Entropy  (J/mol/K)", entrLabel);

        // Push rows to top
        GridBagConstraints gc = new GridBagConstraints();
        gc.gridx = 0; gc.gridy = 6; gc.weighty = 1.0;
        gc.gridwidth = 2;
        grid.add(Box.createGlue(), gc);

        JScrollPane scroll = new JScrollPane(grid);
        scroll.setBorder(null);
        scroll.setBackground(BG);
        scroll.getViewport().setBackground(BG);
        panel.add(scroll, BorderLayout.CENTER);

        return panel;
    }

    private void addRow(JPanel p, GridBagConstraints lc, GridBagConstraints vc,
                        int row, String name, JLabel value) {
        lc.gridx = 0; lc.gridy = row;
        vc.gridx = 1; vc.gridy = row;
        JLabel lbl = new JLabel(name);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        lbl.setForeground(LABEL_DIM);
        p.add(lbl, lc);
        p.add(value, vc);
    }

    private static JLabel makeValueLabel(String text) {
        JLabel l = new JLabel(text);
        l.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 13));
        l.setForeground(VALUE_FG);
        return l;
    }

    // =========================================================================
    // Log/Output section (bottom)
    // =========================================================================

    private JPanel buildLogSection() {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(BG);
        panel.add(makeSectionHeader("OUTPUT"), BorderLayout.NORTH);

        logArea.setEditable(false);
        logArea.setBackground(BG);
        logArea.setForeground(LOG_FG);
        logArea.setCaretColor(LOG_FG);
        logArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        logArea.setLineWrap(true);
        logArea.setWrapStyleWord(true);
        logArea.setBorder(BorderFactory.createEmptyBorder(6, 12, 6, 12));
        logArea.setText("Ready — select a section from the activity bar to begin.\n");

        JScrollPane scroll = new JScrollPane(logArea);
        scroll.setBorder(null);
        scroll.setBackground(BG);
        scroll.getViewport().setBackground(BG);
        panel.add(scroll, BorderLayout.CENTER);

        return panel;
    }

    // =========================================================================
    // Section header helper
    // =========================================================================

    private JPanel makeSectionHeader(String title) {
        JPanel hdr = new JPanel(new BorderLayout());
        hdr.setBackground(HDR_BG);
        hdr.setPreferredSize(new Dimension(0, 28));
        hdr.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 0, 1, 0, BORDER_CLR),
                BorderFactory.createEmptyBorder(0, 12, 0, 0)));
        JLabel lbl = new JLabel(title);
        lbl.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        lbl.setForeground(HDR_FG);
        hdr.add(lbl, BorderLayout.WEST);
        return hdr;
    }

    // =========================================================================
    // Public API
    // =========================================================================

    /**
     * Appends a line to the Output log. Thread-safe — can be called from any thread.
     *
     * @param text text to append (newline appended automatically if not present)
     */
    public void appendLog(String text) {
        String line = text.endsWith("\n") ? text : text + "\n";
        SwingUtilities.invokeLater(() -> {
            logArea.append(line);
            logArea.setCaretPosition(logArea.getDocument().getLength());
        });
    }

    /** Clears the Output log. Thread-safe. */
    public void clearLog() {
        SwingUtilities.invokeLater(() -> logArea.setText(""));
    }

    /**
     * Populates the Results section with a completed calculation result.
     * Thread-safe — can be called from any thread.
     */
    public void showResult(ThermodynamicResult result) {
        SwingUtilities.invokeLater(() -> {
            tempLabel.setText(String.format("%.1f", result.temperature));
            compLabel.setText(result.composition != null && result.composition.length > 1
                    ? String.format("%.4f", result.composition[1]) : "—");
            gibbsLabel.setText(String.format("%.6f", result.gibbsEnergy));
            enthlLabel.setText(String.format("%.6f", result.enthalpy));
            entrLabel.setText(String.format("%.6f", result.entropy));
        });
    }

    /** Resets all result fields to placeholder dashes. Thread-safe. */
    public void clearResult() {
        SwingUtilities.invokeLater(() -> {
            tempLabel.setText("—");
            compLabel.setText("—");
            gibbsLabel.setText("—");
            enthlLabel.setText("—");
            entrLabel.setText("—");
        });
    }
}
