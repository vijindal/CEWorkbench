package org.ce.ui.gui;

import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace.SystemId;
import org.ce.calculation.CalculationDescriptor.ModelSpecifications;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.util.List;
import java.util.function.Consumer;

/**
 * Persistent 40 px session toolbar shown between the HeaderBar and the main content area.
 *
 * <p>Single owner of {@link ModelSession.Builder#build}. All calculation-mode panels
 * (Single Point, Line Scan, Map) rely on the session built here — they only read
 * {@link WorkbenchContext#getActiveSession()} and listen to
 * {@link WorkbenchContext#addSessionListener}.</p>
 *
 * <pre>
 *  [Elements▼] [Structure▼] [Model▼]  Engine:[CVM▼]  [Rebuild]  ● Ready: Nb-Ti/BCC_A2/T
 * </pre>
 *
 * <p>The status dot colour:
 * <ul>
 *   <li>Red   — no session / session invalidated</li>
 *   <li>Yellow — build in progress</li>
 *   <li>Teal  — session ready</li>
 * </ul>
 */
public class SessionBar extends JPanel {

    // ── colours ──────────────────────────────────────────────────────────────
    private static final Color BG          = new Color(0x2D2D2D);
    private static final Color ACCENT_BLUE = new Color(0x0E639C);
    private static final Color DOT_READY   = new Color(0x4EC9B0);
    private static final Color DOT_BUSY    = new Color(0xDCDCAA);
    private static final Color DOT_NONE    = new Color(0xF44747);

    // ── widgets ───────────────────────────────────────────────────────────────
    private final JComboBox<String> elementsCombo  = makeEditable("Nb-Ti", "Cu-Au", "Al-Ti", "Nb-Mo");
    private final JComboBox<String> structureCombo = makeEditable("BCC_A2", "FCC_A1", "HCP_A3");
    private final JComboBox<String> modelCombo     = makeEditable("T", "T2");
    private final JComboBox<String> engineCombo    = new JComboBox<>(new String[]{"CVM", "MCS"});

    private final JButton rebuildBtn = new JButton("Rebuild");
    private final JLabel  statusLabel = new JLabel("No session");
    private final DotIcon dotIcon    = new DotIcon(DOT_NONE);

    // ── state ─────────────────────────────────────────────────────────────────
    private final org.ce.CEWorkbenchContext appCtx;
    private final WorkbenchContext          context;
    private final Consumer<String>          statusSink;
    private final Consumer<String>          logSink;
    private boolean                         building = false;

    public SessionBar(org.ce.CEWorkbenchContext appCtx,
                      WorkbenchContext context,
                      Consumer<String> statusSink,
                      Consumer<String> logSink) {
        this.appCtx     = appCtx;
        this.context    = context;
        this.statusSink = statusSink;
        this.logSink    = logSink;

        setBackground(BG);
        setPreferredSize(new Dimension(0, 40));
        setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(0x444444)));
        setLayout(new FlowLayout(FlowLayout.LEFT, 6, 6));

        // ── label helpers ─────────────────────────────────────────────────────
        add(dimLabel("Elements:"));
        add(elementsCombo);
        add(dimLabel("Structure:"));
        add(structureCombo);
        add(dimLabel("Model:"));
        add(modelCombo);
        add(sep());
        add(dimLabel("Engine:"));
        add(engineCombo);
        add(sep());

        styleButton(rebuildBtn);
        rebuildBtn.addActionListener(e -> rebuildSession());
        add(rebuildBtn);
        add(sep());

        JLabel dotLabel = new JLabel(dotIcon);
        add(dotLabel);
        statusLabel.setForeground(new Color(0x9CDCFE));
        statusLabel.setFont(statusLabel.getFont().deriveFont(11f));
        add(statusLabel);

        // ── push combo edits → WorkbenchContext ───────────────────────────────
        DocumentListener pushSystem = new DocumentListener() {
            public void insertUpdate(DocumentEvent e) { pushToContext(); }
            public void removeUpdate(DocumentEvent e) { pushToContext(); }
            public void changedUpdate(DocumentEvent e) { pushToContext(); }
        };
        editorDoc(elementsCombo).addDocumentListener(pushSystem);
        editorDoc(structureCombo).addDocumentListener(pushSystem);
        editorDoc(modelCombo).addDocumentListener(pushSystem);

        elementsCombo.addActionListener(e -> pushToContext());
        structureCombo.addActionListener(e -> pushToContext());
        modelCombo.addActionListener(e -> pushToContext());

        // ── pull system changes from other panels ─────────────────────────────
        context.addChangeListener(this::syncCombosFromContext);

        // ── react to session lifecycle ─────────────────────────────────────────
        context.addSessionListener(this::onSessionChanged);

        // Push initial defaults
        pushToContext();
    }

    // =========================================================================
    // Session rebuild
    // =========================================================================

    private void rebuildSession() {
        if (building) return;
        String el   = editorText(elementsCombo);
        String str  = editorText(structureCombo);
        String mod  = editorText(modelCombo);
        String eng  = (String) engineCombo.getSelectedItem();
        
        if (el.isBlank() || str.isBlank() || mod.isBlank()) {
            logSink.accept("SessionBar: Elements, Structure, and Model must be specified.");
            return;
        }

        ModelSpecifications specs = new ModelSpecifications(el, str, mod, EngineConfig.valueOf(eng));

        building = true;
        rebuildBtn.setEnabled(false);
        dotIcon.setColor(DOT_BUSY);
        statusLabel.setText("Building…");
        repaint();
        statusSink.accept("Building session [" + eng + "]…");

        SwingWorker<ModelSession, String> worker = new SwingWorker<>() {
            @Override
            protected ModelSession doInBackground() throws Exception {
                // Calculation Layer Role: Model Construction
                return appCtx.getCalculationService().getOrBuildSession(specs, this::publish);
            }
            @Override
            protected void process(List<String> chunks) { chunks.forEach(logSink); }
            @Override
            protected void done() {
                building = false;
                rebuildBtn.setEnabled(true);
                try {
                    ModelSession session = get();
                    context.setActiveSession(session);
                    statusSink.accept("Session ready — " + session.label());
                } catch (Exception ex) {
                    String msg = ex.getCause() != null ? ex.getCause().getMessage() : ex.getMessage();
                    logSink.accept("Session build failed: " + msg);
                    statusSink.accept("Error: " + msg);
                    dotIcon.setColor(DOT_NONE);
                    statusLabel.setText("Build failed");
                    repaint();
                }
            }
        };
        worker.execute();
    }

    // =========================================================================
    // Listeners
    // =========================================================================

    private void pushToContext() {
        String el  = editorText(elementsCombo);
        String str = editorText(structureCombo);
        String mod = editorText(modelCombo);
        if (!el.isBlank() && !str.isBlank() && !mod.isBlank()) {
            context.setSystem(el, str, mod);
        }
    }

    private void syncCombosFromContext() {
        SystemId sys = context.getSystem();
        if (sys == null) return;
        setEditorText(elementsCombo,  sys.elements);
        setEditorText(structureCombo, sys.structure);
        setEditorText(modelCombo,     sys.model);
    }

    private void onSessionChanged(ModelSession session) {
        if (session != null) {
            dotIcon.setColor(DOT_READY);
            statusLabel.setText("Ready: " + session.label());
            // Sync engine combo to actual session engine
            engineCombo.setSelectedItem(session.engineConfig.name());
        } else {
            dotIcon.setColor(DOT_NONE);
            statusLabel.setText("No session");
        }
        repaint();
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private static JComboBox<String> makeEditable(String... items) {
        JComboBox<String> cb = new JComboBox<>(items);
        cb.setEditable(true);
        cb.setPreferredSize(new Dimension(90, 24));
        return cb;
    }

    private static JLabel dimLabel(String text) {
        JLabel l = new JLabel(text);
        l.setForeground(new Color(0x858585));
        l.setFont(l.getFont().deriveFont(11f));
        return l;
    }

    private static JSeparator sep() {
        JSeparator s = new JSeparator(SwingConstants.VERTICAL);
        s.setPreferredSize(new Dimension(1, 24));
        s.setForeground(new Color(0x444444));
        return s;
    }

    private static void styleButton(JButton btn) {
        btn.setBackground(ACCENT_BLUE);
        btn.setForeground(Color.WHITE);
        btn.setFocusPainted(false);
        btn.setFont(btn.getFont().deriveFont(11f));
        btn.setPreferredSize(new Dimension(80, 24));
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
