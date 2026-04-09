package org.ce.ui.gui;

import static org.ce.model.cluster.ClusterPrimitives.*;

import org.ce.calculation.QuantityDescriptor;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECTerm;
import org.ce.calculation.engine.ProgressEvent;
import org.ce.model.result.ThermodynamicResult;
import org.ce.model.storage.Workspace.SystemId;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * VS Code-style main output panel (third column, fills remaining width).
 *
 * <p>Divided vertically by a draggable splitter:</p>
 * <ul>
 *   <li><b>Top — RESULTS</b>: live chart of engine progress plus a scalar result strip.
 *       MCS shows energy-per-site vs sweep; CVM shows log₁₀|∇G| vs N-R iteration.
 *       Correlation functions (CFs) at each step are carried in the event and also
 *       reported in the log.</li>
 *   <li><b>Bottom — OUTPUT</b>: dark terminal-style log for text progress messages.</li>
 * </ul>
 */
public class OutputPanel extends JPanel {

    // ── colours ───────────────────────────────────────────────────────────────
    private static final Color BG         = new Color(0x1E1E1E);
    private static final Color HDR_BG     = new Color(0x2D2D2D);
    private static final Color HDR_FG     = new Color(0xBBBBBB);
    private static final Color BORDER_CLR = new Color(0x2D2D2D);
    private static final Color LOG_FG     = new Color(0xCCCCCC);

    private final ResultChartPanel chartPanel;
    private final CECEditorPanel cecEditorPanel;
    private final JPanel resultsBody;
    private final CardLayout resultsBodyLayout;
    private static final String CARD_CHART = "chart";
    private static final String CARD_CEC = "cec";
    private final JTextArea        logArea = new JTextArea();

    public OutputPanel(WorkbenchContext context) {
        setLayout(new BorderLayout());
        setBackground(BG);

        chartPanel = new ResultChartPanel(context);
        cecEditorPanel = new CECEditorPanel();
        resultsBodyLayout = new CardLayout();
        resultsBody = new JPanel(resultsBodyLayout);
        resultsBody.setBackground(BG);
        resultsBody.add(chartPanel, CARD_CHART);
        resultsBody.add(cecEditorPanel, CARD_CEC);

        JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                buildResultsSection(),
                buildLogSection());
        split.setDividerLocation(320);
        split.setDividerSize(1);
        split.setContinuousLayout(true);
        split.setBorder(null);
        split.setBackground(BG);

        add(split, BorderLayout.CENTER);
    }

    // =========================================================================
    // Results section (top) — chart + scalar strip
    // =========================================================================

    private JPanel buildResultsSection() {
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBackground(BG);
        panel.add(makeSectionHeader("RESULTS"), BorderLayout.NORTH);
        panel.add(resultsBody, BorderLayout.CENTER);
        return panel;
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

    /** Appends a line to the Output log. Thread-safe. */
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
     * Routes a structured engine event to the chart panel and logs iteration details.
     * Called on the EDT via SwingWorker.process() — no extra invokeLater needed.
     */
    public void onChartEvent(ProgressEvent event) {
        resultsBodyLayout.show(resultsBody, CARD_CHART);
        // Log per-iteration/sweep details with H, S, and CFs
        if (event instanceof ProgressEvent.CvmIteration) {
            ProgressEvent.CvmIteration ci = (ProgressEvent.CvmIteration) event;
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("  iter %3d  |∇G| = %.3e  G = %.4f  H = %.4f  S = %.6f",
                                   ci.iteration, ci.gradientNorm, ci.gibbsEnergy,
                                   ci.enthalpy, ci.entropy));
            if (ci.cfs != null && ci.cfs.length > 0) {
                sb.append("  CFs: [");
                for (int i = 0; i < ci.cfs.length; i++) {
                    if (i > 0) sb.append(", ");
                    sb.append(String.format("%.4f", ci.cfs[i]));
                }
                sb.append("]");
            }
            appendLog(sb.toString());
        } else if (event instanceof ProgressEvent.McSweep) {
            ProgressEvent.McSweep ms = (ProgressEvent.McSweep) event;
            // Only log averaging sweeps with CFs; equilibration is too noisy
            if (!ms.equilibration && ms.cfs != null && ms.cfs.length > 0) {
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("  sweep %3d  E/site = %.4f  acc = %.3f",
                                       ms.step, ms.energyPerSite, ms.acceptanceRate));
                sb.append("  CFs: [");
                for (int i = 0; i < ms.cfs.length; i++) {
                    if (i > 0) sb.append(", ");
                    sb.append(String.format("%.4f", ms.cfs[i]));
                }
                sb.append("]");
                appendLog(sb.toString());
            }
        }
        chartPanel.onEvent(event);
    }

    /**
     * Populates the Results strip with a completed calculation result.
     * Thread-safe — can be called from any thread.
     */
    public void showResult(ThermodynamicResult result) {
        SwingUtilities.invokeLater(() -> {
            resultsBodyLayout.show(resultsBody, CARD_CHART);
            chartPanel.showResult(result);
        });
    }

    /**
     * Shows the results of a line scan (T-sweep or x-sweep) in the chart panel.
     * Thread-safe — can be called from any thread.
     *
     * @param results  list of thermodynamic results (one per scan point)
     * @param scanType "T" for temperature scan, "X" for composition scan
     */
    public void showLineScanResult(List<ThermodynamicResult> results, String scanType) {
        SwingUtilities.invokeLater(() -> {
            resultsBodyLayout.show(resultsBody, CARD_CHART);
            chartPanel.showLineScanResult(results, scanType);
        });
    }

    /**
     * Shows the results of a 2D grid scan (T × x) as a heatmap.
     * Thread-safe — can be called from any thread.
     *
     * @param grid rows indexed by temperature (ascending), each row by composition (ascending)
     */
    public void showMapResult(List<List<ThermodynamicResult>> grid) {
        SwingUtilities.invokeLater(() -> {
            resultsBodyLayout.show(resultsBody, CARD_CHART);
            chartPanel.showMapResult(grid);
        });
    }

    /** Resets the chart. Thread-safe. */
    public void clearResult() {
        SwingUtilities.invokeLater(() -> {
            chartPanel.clear();
            resultsBodyLayout.show(resultsBody, CARD_CHART);
        });
    }

    /** Shows orthogonal and CVCF CECs loaded directly from JSON files. */
    public void showCECResult(CECEntry orthEntry, CECEntry cvcfEntry) {
        SwingUtilities.invokeLater(() -> {
            cecEditorPanel.showEntries(orthEntry, cvcfEntry);
            resultsBodyLayout.show(resultsBody, CARD_CEC);
        });
    }

    /** Applies current edits from the CEC results editor back into the entry. */
    public boolean applyCECEdits(CECEntry entry) {
        if (entry == null) return false;
        try {
            cecEditorPanel.applyEdits(entry);
            return true;
        } catch (Exception ex) {
            appendLog("Error: " + ex.getMessage());
            return false;
        }
    }

    private static final class CECEditorPanel extends JPanel {
        private final DefaultTableModel cvcfModel = new DefaultTableModel(
                new Object[]{"Name", "a (J/mol)", "b (J/mol/K)"}, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return column > 0;
            }
        };
        private final DefaultTableModel orthoModel = new DefaultTableModel(
                new Object[]{"Name", "a (J/mol)", "b (J/mol/K)"}, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return column > 0;
            }
        };
        private final JTable cvcfTable = new JTable(cvcfModel);
        private final JTable orthoTable = new JTable(orthoModel);
        private final JLabel title = new JLabel("CEC Editor");
        private CECEntry currentOrthEntry = null;
        private CECEntry currentCvcfEntry = null;

        CECEditorPanel() {
            setLayout(new BorderLayout(6, 6));
            setBackground(new Color(0x1E1E1E));

            title.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 12));
            title.setForeground(new Color(0xBBBBBB));
            title.setBorder(BorderFactory.createEmptyBorder(6, 8, 0, 0));
            add(title, BorderLayout.NORTH);

            configureTable(cvcfTable);
            configureTable(orthoTable);

            JTabbedPane tabs = new JTabbedPane();
            tabs.addTab("CVCF Basis (Editable)", new JScrollPane(cvcfTable));
            tabs.addTab("Orthogonal Basis (Derived)", new JScrollPane(orthoTable));
            add(tabs, BorderLayout.CENTER);
        }

        private void configureTable(JTable table) {
            table.setFillsViewportHeight(true);
            table.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
            table.setBackground(new Color(0x252526));
            table.setForeground(new Color(0xD4D4D4));
            table.setGridColor(new Color(0x3C3C3C));
            table.setSelectionBackground(new Color(0x094771));
            table.setSelectionForeground(Color.WHITE);
            table.getTableHeader().setBackground(new Color(0x2D2D2D));
            table.getTableHeader().setForeground(new Color(0xBBBBBB));
        }

        void showEntries(CECEntry orthEntry, CECEntry cvcfEntry) {
            this.currentOrthEntry = orthEntry;
            this.currentCvcfEntry = cvcfEntry;
            cvcfModel.setRowCount(0);
            orthoModel.setRowCount(0);
            if (orthEntry == null && cvcfEntry == null) {
                title.setText("CEC Editor");
                return;
            }
            CECEntry ref = cvcfEntry != null ? cvcfEntry : orthEntry;
            title.setText("CEC Editor  |  " + ref.elements + " / " + ref.structurePhase + " / " + ref.model);

            if (orthEntry != null && orthEntry.cecTerms != null) {
                for (CECTerm term : orthEntry.cecTerms) {
                    orthoModel.addRow(new Object[]{term.name, term.a, term.b});
                }
            }
            if (cvcfEntry != null && cvcfEntry.cecTerms != null) {
                for (CECTerm term : cvcfEntry.cecTerms) {
                    cvcfModel.addRow(new Object[]{term.name, term.a, term.b});
                }
            }
        }

        void applyEdits(CECEntry entry) {
            if (entry == null || entry.cecTerms == null) {
                return;
            }
            boolean useCvcf = (entry == currentCvcfEntry) || (currentCvcfEntry != null && currentOrthEntry == null);
            JTable table = useCvcf ? cvcfTable : orthoTable;
            DefaultTableModel model = useCvcf ? cvcfModel : orthoModel;
            if (table.isEditing()) {
                table.getCellEditor().stopCellEditing();
            }
            int n = Math.min(entry.cecTerms.length, model.getRowCount());
            for (int i = 0; i < n; i++) {
                Object aVal = model.getValueAt(i, 1);
                Object bVal = model.getValueAt(i, 2);
                try {
                    entry.cecTerms[i].a = Double.parseDouble(aVal.toString());
                    entry.cecTerms[i].b = Double.parseDouble(bVal.toString());
                } catch (Exception ex) {
                    throw new IllegalArgumentException("Invalid number at row " + (i + 1));
                }
            }
        }
    }

    // =========================================================================
    // Inner class: ResultChartPanel
    // =========================================================================

    private static class ResultChartPanel extends JPanel {

        // ── chart colours ─────────────────────────────────────────────────────
        private static final Color BG        = new Color(0x1E1E1E);
        private static final Color PLOT_BG   = new Color(0x252526);
        private static final Color GRID_CLR  = new Color(0x2D2D2D);
        private static final Color AXIS_FG   = new Color(0x858585);
        private static final Color EQUIL_CLR = new Color(0x569CD6);  // VS Code blue
        private static final Color AVG_CLR   = new Color(0x4EC9B0);  // teal
        private static final Color CVM_CLR   = new Color(0xDCDCAA);  // yellow
        private static final Color BORDER    = new Color(0x3C3C3C);
        private static final Color STRIP_BG  = new Color(0x252526);
        private static final Color LABEL_DIM = new Color(0x858585);
        private static final Color VALUE_FG  = new Color(0xD4D4D4);

        // ── chart margins ─────────────────────────────────────────────────────
        private static final int ML = 64, MB = 34, MT = 22, MR = 16;
        private static final int STRIP_H = 54;

        private enum Mode { IDLE, MCS, CVM, LINESCAN, MAP }

        private Mode   mode = Mode.IDLE;

        // MCS data: [step, energyPerSite]
        private final List<double[]> equilData = new ArrayList<>();
        private final List<double[]> avgData   = new ArrayList<>();
        // CVM data: [iteration, log10(gradNorm)]
        private final List<double[]> cvmData   = new ArrayList<>();
        // CVM thermodynamic properties on secondary axis: [iteration, value]
        private final List<double[]> cvmGData    = new ArrayList<>();
        private final List<double[]> cvmHData    = new ArrayList<>();
        private final List<double[]> cvmNegTSData = new ArrayList<>();

        private ThermodynamicResult lastResult = null;
        private String systemText = "no system selected";
        private String textResult = null;

        // Line-scan data: list of [xValue, G, H, S] per scan point
        private final List<double[]> lsData = new ArrayList<>();
        private String lsScanType = "T";  // "T" = temperature scan, "X" = composition scan

        // Map data: rows indexed by temperature (ascending), columns by composition (ascending)
        private List<List<ThermodynamicResult>> mapGrid = new ArrayList<>();

        ResultChartPanel(WorkbenchContext context) {
            setBackground(BG);
            context.addChangeListener(() -> {
                if (context.hasSystem()) {
                    SystemId s = context.getSystem();
                    systemText = s.elements + "  /  " + s.structure + "  /  " + s.model;
                } else {
                    systemText = "no system selected";
                }
                repaint();
            });
        }

        // ── public API ────────────────────────────────────────────────────────

        void clear() {
            mode = Mode.IDLE;
            equilData.clear();
            avgData.clear();
            cvmData.clear();
            cvmGData.clear();
            cvmHData.clear();
            cvmNegTSData.clear();
            lsData.clear();
            mapGrid = new ArrayList<>();
            lastResult = null;
            textResult = null;
            repaint();
        }

        void onEvent(ProgressEvent evt) {
            if (evt instanceof ProgressEvent.EngineStart) {
                ProgressEvent.EngineStart es = (ProgressEvent.EngineStart) evt;
                clear();
                mode = "MCS".equals(es.engineType) ? Mode.MCS : Mode.CVM;
            } else if (evt instanceof ProgressEvent.McSweep) {
                ProgressEvent.McSweep ms = (ProgressEvent.McSweep) evt;
                if (ms.equilibration) equilData.add(new double[]{ms.step, ms.energyPerSite});
                else                  avgData.add(new double[]{ms.step, ms.energyPerSite});
            } else if (evt instanceof ProgressEvent.CvmIteration) {
                ProgressEvent.CvmIteration ci = (ProgressEvent.CvmIteration) evt;
                double logNorm = ci.gradientNorm > 0 ? Math.log10(ci.gradientNorm) : -15.0;
                cvmData.add(new double[]{ci.iteration, logNorm});
                cvmGData.add(new double[]{ci.iteration, ci.gibbsEnergy});
                cvmHData.add(new double[]{ci.iteration, ci.enthalpy});
                cvmNegTSData.add(new double[]{ci.iteration, ci.entropy});
            }
            repaint();
        }

        void showResult(ThermodynamicResult result) {
            textResult = null;
            lastResult = result;
            repaint();
        }

        void showTextResult(String text) {
            mode = Mode.IDLE;
            equilData.clear();
            avgData.clear();
            cvmData.clear();
            cvmGData.clear();
            cvmHData.clear();
            cvmNegTSData.clear();
            lsData.clear();
            lastResult = null;
            textResult = text;
            repaint();
        }

        void showLineScanResult(List<ThermodynamicResult> results, String scanType) {
            lsData.clear();
            lsScanType = scanType;
            for (ThermodynamicResult r : results) {
                double xVal = "T".equals(scanType)
                        ? r.temperature
                        : (r.composition.length > 1 ? r.composition[1] : r.composition[0]);
                lsData.add(new double[]{xVal, r.gibbsEnergy, r.enthalpy, r.entropy});
            }
            mode = Mode.LINESCAN;
            lastResult = results.isEmpty() ? null : results.get(results.size() - 1);
            repaint();
        }

        void showMapResult(List<List<ThermodynamicResult>> grid) {
            mapGrid = grid;
            mode = Mode.MAP;
            // Set lastResult to the centre point for the strip display
            if (!grid.isEmpty()) {
                List<ThermodynamicResult> midRow = grid.get(grid.size() / 2);
                lastResult = midRow.isEmpty() ? null : midRow.get(midRow.size() / 2);
            }
            repaint();
        }

        // ── painting ──────────────────────────────────────────────────────────

        @Override
        protected void paintComponent(Graphics g0) {
            Graphics2D g = (Graphics2D) g0;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,      RenderingHints.VALUE_ANTIALIAS_ON);
            g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int w = getWidth(), h = getHeight();
            g.setColor(BG);
            g.fillRect(0, 0, w, h);

            int chartH = h - STRIP_H;
            drawChart(g, w, chartH);
            drawStrip(g, 0, chartH, w, STRIP_H);
        }

        // ── chart area ────────────────────────────────────────────────────────

        private void drawChart(Graphics2D g, int w, int h) {
            // For CVM, show secondary axis with G, H, -TS; for MCS, no secondary axis
            boolean hasThermo = mode == Mode.CVM && !cvmGData.isEmpty();
            int mr = hasThermo ? 60 : MR;

            int px = ML, py = MT, pw = w - ML - mr, ph = h - MT - MB;
            if (pw <= 0 || ph <= 0) return;

            // Plot background
            g.setColor(PLOT_BG);
            g.fillRect(px, py, pw, ph);

            // ── Line scan mode ────────────────────────────────────────────────
            if (mode == Mode.LINESCAN) {
                drawLineScanChart(g, px, py, pw, ph, w);
                g.setColor(BORDER);
                g.drawRect(px, py, pw, ph);
                return;
            }

            // ── Map mode ──────────────────────────────────────────────────────
            if (mode == Mode.MAP) {
                drawMapChart(g, px, py, pw, ph);
                g.setColor(BORDER);
                g.drawRect(px, py, pw, ph);
                return;
            }

            if (mode == Mode.IDLE) {
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
                g.setColor(LABEL_DIM);
                if (textResult == null || textResult.isBlank()) {
                    FontMetrics fm = g.getFontMetrics();
                    String msg = "Run a calculation to see the convergence chart";
                    g.drawString(msg, px + (pw - fm.stringWidth(msg)) / 2, py + ph / 2 + 4);
                } else {
                    g.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
                    FontMetrics fm = g.getFontMetrics();
                    int y = py + 16;
                    for (String line : textResult.split("\\R")) {
                        if (y > py + ph - 6) break;
                        g.drawString(line, px + 10, y);
                        y += fm.getHeight();
                    }
                }
                g.setColor(BORDER);
                g.drawRect(px, py, pw, ph);
                return;
            }

            // Gather data
            List<double[]> primary   = mode == Mode.MCS ? equilData : new ArrayList<>();
            List<double[]> secondary = mode == Mode.MCS ? avgData   : cvmData;
            List<double[]> all       = new ArrayList<>(primary);
            all.addAll(secondary);

            if (all.isEmpty()) {
                g.setColor(BORDER);
                g.drawRect(px, py, pw, ph);
                return;
            }

            // Compute axis ranges
            double xMin = all.get(0)[0], xMax = all.get(all.size() - 1)[0];
            double yMin = Double.MAX_VALUE, yMax = -Double.MAX_VALUE;
            for (double[] pt : all) {
                if (pt[1] < yMin) yMin = pt[1];
                if (pt[1] > yMax) yMax = pt[1];
            }
            if (xMax == xMin) xMax = xMin + 1;
            double yPad = (yMax - yMin) * 0.12;
            if (yPad == 0) yPad = 0.001;
            yMin -= yPad; yMax += yPad;

            // Compute secondary Y axis range (CVM Gibbs energy only) if available
            double thermoMin = Double.MAX_VALUE, thermoMax = -Double.MAX_VALUE;
            if (hasThermo) {
                for (double[] pt : cvmGData) {
                    if (pt[1] < thermoMin) thermoMin = pt[1];
                    if (pt[1] > thermoMax) thermoMax = pt[1];
                }
                double thermoPad = (thermoMax - thermoMin) * 0.12;
                if (thermoPad == 0) thermoPad = 1.0;
                thermoMin -= thermoPad; thermoMax += thermoPad;
            }

            // Grid
            g.setStroke(new BasicStroke(1f));
            g.setColor(GRID_CLR);
            int nGrid = 5;
            for (int i = 0; i <= nGrid; i++) {
                double frac = (double) i / nGrid;
                g.drawLine(px, py + (int)(frac * ph), px + pw, py + (int)(frac * ph));
                g.drawLine(px + (int)(frac * pw), py, px + (int)(frac * pw), py + ph);
            }

            // Y-axis labels (left)
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            g.setColor(AXIS_FG);
            FontMetrics fm = g.getFontMetrics();
            for (int i = 0; i <= nGrid; i++) {
                double val = yMax - (double) i / nGrid * (yMax - yMin);
                int    gy  = py + (int)((double) i / nGrid * ph);
                String lbl = formatAxisVal(val);
                g.drawString(lbl, px - fm.stringWidth(lbl) - 4, gy + fm.getAscent() / 2);
            }

            // Y-axis labels (right, for CVM thermodynamic secondary axis)
            if (hasThermo) {
                for (int i = 0; i <= nGrid; i++) {
                    double val = thermoMax - (double) i / nGrid * (thermoMax - thermoMin);
                    int    gy  = py + (int)((double) i / nGrid * ph);
                    String lbl = formatAxisVal(val);
                    g.drawString(lbl, px + pw + 4, gy + fm.getAscent() / 2);
                }
            }

            // X-axis labels
            for (int i = 0; i <= nGrid; i++) {
                double val = xMin + (double) i / nGrid * (xMax - xMin);
                int    gx  = px + (int)((double) i / nGrid * pw);
                String lbl = String.format("%.0f", val);
                g.drawString(lbl, gx - fm.stringWidth(lbl) / 2, py + ph + 14);
            }

            // Axis titles
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = g.getFontMetrics();
            String xTitle = mode == Mode.MCS ? "Sweep" : "Iteration";
            String yTitle = mode == Mode.MCS ? "E / site" : "log\u2081\u2080 |\u2207G|";
            g.drawString(xTitle, px + (pw - fm.stringWidth(xTitle)) / 2, py + ph + 28);

            // Rotated Y title (left)
            Graphics2D g2 = (Graphics2D) g.create();
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g2.translate(10, py + ph / 2);
            g2.rotate(-Math.PI / 2);
            g2.setColor(AXIS_FG);
            g2.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = g2.getFontMetrics();
            g2.drawString(yTitle, -fm.stringWidth(yTitle) / 2, fm.getAscent() / 2);
            g2.dispose();

            // Rotated Y title (right, for CVM thermodynamic secondary axis)
            if (hasThermo) {
                Graphics2D g3 = (Graphics2D) g.create();
                g3.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                g3.translate(w - 10, py + ph / 2);
                g3.rotate(Math.PI / 2);
                g3.setColor(AXIS_FG);
                g3.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
                fm = g3.getFontMetrics();
                g3.drawString("J/mol", -fm.stringWidth("J/mol") / 2, fm.getAscent() / 2);
                g3.dispose();
            }

            // Chart title (top-centre)
            g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
            g.setColor(AXIS_FG);
            fm = g.getFontMetrics();
            String title = mode == Mode.MCS
                    ? "Monte Carlo — Energy per Site  [" + systemText + "]"
                    : "CVM Newton-Raphson Convergence  [" + systemText + "]";
            String shortened = title.length() > 72 ? title.substring(0, 69) + "..." : title;
            g.drawString(shortened, px + (pw - fm.stringWidth(shortened)) / 2, py - 6);

            // MCS phase boundary (dashed vertical line)
            if (mode == Mode.MCS && !avgData.isEmpty()) {
                int bx = toX(avgData.get(0)[0], xMin, xMax, px, pw);
                g.setColor(new Color(0x555555));
                float[] dash = {4f, 4f};
                g.setStroke(new BasicStroke(1f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10f, dash, 0f));
                g.drawLine(bx, py, bx, py + ph);
                g.setStroke(new BasicStroke(1f));
            }

            // Clip to plot area for polylines
            Shape oldClip = g.getClip();
            g.setClip(px, py, pw, ph);

            drawSeries(g, primary,   EQUIL_CLR, px, py, pw, ph, xMin, xMax, yMin, yMax);
            drawSeries(g, secondary, mode == Mode.MCS ? AVG_CLR : CVM_CLR, px, py, pw, ph, xMin, xMax, yMin, yMax);

            // Draw CVM Gibbs energy (G) on secondary Y axis
            if (hasThermo) {
                drawSeriesOnSecondaryAxis(g, cvmGData, CVM_CLR, px, py, pw, ph, xMin, xMax, thermoMin, thermoMax);
            }

            g.setClip(oldClip);

            // Plot border
            g.setStroke(new BasicStroke(1f));
            g.setColor(BORDER);
            g.drawRect(px, py, pw, ph);

            // Legend (top-right inside plot)
            drawLegend(g, px + pw - 8, py + 8);
        }

        private void drawMapChart(Graphics2D g, int px, int py, int pw, int ph) {
            if (mapGrid == null || mapGrid.isEmpty()) {
                g.setColor(LABEL_DIM);
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
                FontMetrics fm = g.getFontMetrics();
                String msg = "No map data";
                g.drawString(msg, px + (pw - fm.stringWidth(msg)) / 2, py + ph / 2);
                return;
            }

            int nT = mapGrid.size();
            int nX = mapGrid.get(0).size();
            if (nX == 0) return;

            // Collect all Gibbs energies to find global min/max
            double gMin = Double.MAX_VALUE, gMax = -Double.MAX_VALUE;
            for (List<ThermodynamicResult> row : mapGrid) {
                for (ThermodynamicResult r : row) {
                    if (r.isFreeEnergyValid()) {
                        if (r.gibbsEnergy < gMin) gMin = r.gibbsEnergy;
                        if (r.gibbsEnergy > gMax) gMax = r.gibbsEnergy;
                    }
                }
            }
            if (gMax == gMin) gMax = gMin + 1;

            // T and x axis ranges
            ThermodynamicResult first = mapGrid.get(0).get(0);
            ThermodynamicResult last  = mapGrid.get(nT - 1).get(nX - 1);
            double tMin = first.temperature, tMax = last.temperature;
            double xMin = first.composition.length > 1 ? first.composition[1] : 0;
            double xMax = last.composition.length  > 1 ? last.composition[1]  : 1;

            // Reserve right margin for colour bar
            int cbW  = 18;  // colour bar width
            int cbGap = 8;
            int plotW = pw - cbW - cbGap - 36; // 36 for cb labels
            if (plotW < 10) plotW = pw / 2;

            int cellW = Math.max(1, plotW / nX);
            int cellH = Math.max(1, ph / nT);

            // Draw heatmap cells (row 0 = lowest T at bottom)
            Shape oldClip = g.getClip();
            g.setClip(px, py, plotW, ph);
            for (int ti = 0; ti < nT; ti++) {
                List<ThermodynamicResult> row = mapGrid.get(ti);
                int cy = py + ph - (ti + 1) * cellH;  // low T at bottom
                for (int xi = 0; xi < row.size(); xi++) {
                    ThermodynamicResult r = row.get(xi);
                    double g_val = r.isFreeEnergyValid() ? r.gibbsEnergy : gMin;
                    double t = (g_val - gMin) / (gMax - gMin);  // 0=min(blue), 1=max(red)
                    g.setColor(heatColor(t));
                    g.fillRect(px + xi * cellW, cy, cellW, cellH);
                }
            }
            g.setClip(oldClip);

            // Axis labels
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            g.setColor(AXIS_FG);
            FontMetrics fm = g.getFontMetrics();
            int nGrid = 5;
            for (int i = 0; i <= nGrid; i++) {
                // X axis (composition)
                double xv = xMin + (double) i / nGrid * (xMax - xMin);
                int gx = px + (int)((double) i / nGrid * plotW);
                String xl = String.format("%.2f", xv);
                g.drawString(xl, gx - fm.stringWidth(xl) / 2, py + ph + 14);
                // Y axis (temperature)
                double tv = tMin + (double) i / nGrid * (tMax - tMin);
                int gy = py + ph - (int)((double) i / nGrid * ph);
                String tl = String.format("%.0f", tv);
                g.drawString(tl, px - fm.stringWidth(tl) - 4, gy + fm.getAscent() / 2);
            }

            // Axis titles
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = g.getFontMetrics();
            String xTitle = "x\u2082  (composition)";
            g.drawString(xTitle, px + (plotW - fm.stringWidth(xTitle)) / 2, py + ph + 28);
            Graphics2D gT = (Graphics2D) g.create();
            gT.translate(10, py + ph / 2);
            gT.rotate(-Math.PI / 2);
            gT.setColor(AXIS_FG);
            gT.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = gT.getFontMetrics();
            gT.drawString("T (K)", -fm.stringWidth("T (K)") / 2, fm.getAscent() / 2);
            gT.dispose();

            // Chart title
            g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
            g.setColor(AXIS_FG);
            fm = g.getFontMetrics();
            String title = "G (J/mol)  vs  T \u00d7 x  [" + systemText + "]";
            String shortened = title.length() > 72 ? title.substring(0, 69) + "..." : title;
            g.drawString(shortened, px + (plotW - fm.stringWidth(shortened)) / 2, py - 6);

            // Colour bar
            int cbX = px + plotW + cbGap;
            for (int i = 0; i < ph; i++) {
                double t = 1.0 - (double) i / ph;
                g.setColor(heatColor(t));
                g.fillRect(cbX, py + i, cbW, 1);
            }
            g.setColor(BORDER);
            g.drawRect(cbX, py, cbW, ph);

            // Colour bar labels
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            g.setColor(AXIS_FG);
            fm = g.getFontMetrics();
            for (int i = 0; i <= 4; i++) {
                double val = gMax - (double) i / 4 * (gMax - gMin);
                int ly = py + (int)((double) i / 4 * ph) + fm.getAscent() / 2;
                g.drawString(formatAxisVal(val), cbX + cbW + 3, ly);
            }
        }

        /** Maps t∈[0,1] → blue→white→red colour gradient. */
        private static Color heatColor(double t) {
            t = Math.max(0, Math.min(1, t));
            if (t < 0.5) {
                // blue → white
                float f = (float)(t * 2);
                return new Color(f, f, 1f);
            } else {
                // white → red
                float f = (float)((1 - t) * 2);
                return new Color(1f, f, f);
            }
        }

        private static final Color LS_G_CLR = new Color(0x4EC9B0);  // teal — G
        private static final Color LS_H_CLR = new Color(0xCE9178);  // orange — H
        private static final Color LS_S_CLR = new Color(0x9CDCFE);  // sky — S

        private void drawLineScanChart(Graphics2D g, int px, int py, int pw, int ph, int totalW) {
            if (lsData.isEmpty()) {
                g.setColor(LABEL_DIM);
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
                FontMetrics fm = g.getFontMetrics();
                String msg = "No line-scan data";
                g.drawString(msg, px + (pw - fm.stringWidth(msg)) / 2, py + ph / 2);
                return;
            }

            // Axis ranges for primary Y (G, H in J/mol) — left axis
            double xMin = lsData.get(0)[0], xMax = lsData.get(lsData.size()-1)[0];
            if (xMax == xMin) xMax = xMin + 1;

            double yMin = Double.MAX_VALUE, yMax = -Double.MAX_VALUE;
            for (double[] pt : lsData) {
                if (pt[1] < yMin) yMin = pt[1];  // G
                if (pt[1] > yMax) yMax = pt[1];
                if (pt[2] < yMin) yMin = pt[2];  // H
                if (pt[2] > yMax) yMax = pt[2];
            }
            double pad = (yMax - yMin) * 0.12;
            if (pad == 0) pad = 1.0;
            yMin -= pad; yMax += pad;

            // Secondary Y: entropy (S) — right axis
            double sMin = Double.MAX_VALUE, sMax = -Double.MAX_VALUE;
            for (double[] pt : lsData) {
                if (pt[3] < sMin) sMin = pt[3];
                if (pt[3] > sMax) sMax = pt[3];
            }
            double sPad = (sMax - sMin) * 0.12;
            if (sPad == 0) sPad = 0.001;
            sMin -= sPad; sMax += sPad;

            int nGrid = 5;

            // Grid
            g.setStroke(new BasicStroke(1f));
            g.setColor(GRID_CLR);
            for (int i = 0; i <= nGrid; i++) {
                double frac = (double) i / nGrid;
                g.drawLine(px, py + (int)(frac * ph), px + pw, py + (int)(frac * ph));
                g.drawLine(px + (int)(frac * pw), py, px + (int)(frac * pw), py + ph);
            }

            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            g.setColor(AXIS_FG);
            FontMetrics fm = g.getFontMetrics();

            // Left Y-axis labels (G, H)
            for (int i = 0; i <= nGrid; i++) {
                double val = yMax - (double) i / nGrid * (yMax - yMin);
                int gy = py + (int)((double) i / nGrid * ph);
                String lbl = formatAxisVal(val);
                g.drawString(lbl, px - fm.stringWidth(lbl) - 4, gy + fm.getAscent() / 2);
            }

            // Right Y-axis labels (S)
            int rightLabelX = px + pw + 4;
            for (int i = 0; i <= nGrid; i++) {
                double val = sMax - (double) i / nGrid * (sMax - sMin);
                int gy = py + (int)((double) i / nGrid * ph);
                String lbl = formatAxisVal(val);
                g.drawString(lbl, rightLabelX, gy + fm.getAscent() / 2);
            }

            // X-axis labels
            String xLabel = "T".equals(lsScanType) ? "T (K)" : "x₂";
            for (int i = 0; i <= nGrid; i++) {
                double val = xMin + (double) i / nGrid * (xMax - xMin);
                int gx = px + (int)((double) i / nGrid * pw);
                String lbl = "T".equals(lsScanType) ? String.format("%.0f", val) : String.format("%.3f", val);
                g.drawString(lbl, gx - fm.stringWidth(lbl) / 2, py + ph + 14);
            }

            // Axis titles
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = g.getFontMetrics();
            g.drawString(xLabel, px + (pw - fm.stringWidth(xLabel)) / 2, py + ph + 28);

            Graphics2D gL = (Graphics2D) g.create();
            gL.translate(10, py + ph / 2);
            gL.rotate(-Math.PI / 2);
            gL.setColor(AXIS_FG);
            gL.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = gL.getFontMetrics();
            String leftTitle = "G, H  (J/mol)";
            gL.drawString(leftTitle, -fm.stringWidth(leftTitle) / 2, fm.getAscent() / 2);
            gL.dispose();

            Graphics2D gR = (Graphics2D) g.create();
            gR.translate(totalW - 10, py + ph / 2);
            gR.rotate(Math.PI / 2);
            gR.setColor(AXIS_FG);
            gR.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            fm = gR.getFontMetrics();
            String rightTitle = "S  (J/mol/K)";
            gR.drawString(rightTitle, -fm.stringWidth(rightTitle) / 2, fm.getAscent() / 2);
            gR.dispose();

            // Chart title
            g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
            g.setColor(AXIS_FG);
            fm = g.getFontMetrics();
            String title = "T".equals(lsScanType)
                    ? "G, H, S  vs  Temperature  [" + systemText + "]"
                    : "G, H, S  vs  Composition  [" + systemText + "]";
            String shortened = title.length() > 72 ? title.substring(0, 69) + "..." : title;
            g.drawString(shortened, px + (pw - fm.stringWidth(shortened)) / 2, py - 6);

            // Draw series clipped to plot area
            Shape oldClip = g.getClip();
            g.setClip(px, py, pw, ph);

            // Build point arrays for G, H on left axis; S on right axis
            int n = lsData.size();
            int[] xs = new int[n], yG = new int[n], yH = new int[n], yS = new int[n];
            for (int i = 0; i < n; i++) {
                xs[i] = toX(lsData.get(i)[0], xMin, xMax, px, pw);
                yG[i] = toY(lsData.get(i)[1], yMin, yMax, py, ph);
                yH[i] = toY(lsData.get(i)[2], yMin, yMax, py, ph);
                yS[i] = toYR(lsData.get(i)[3], sMin, sMax, py, ph);
            }

            g.setStroke(new BasicStroke(1.8f));
            g.setColor(LS_H_CLR);
            if (n == 1) g.fillOval(xs[0]-2, yH[0]-2, 4, 4);
            else        g.drawPolyline(xs, yH, n);

            g.setColor(LS_G_CLR);
            if (n == 1) g.fillOval(xs[0]-2, yG[0]-2, 4, 4);
            else        g.drawPolyline(xs, yG, n);

            g.setColor(LS_S_CLR);
            g.setStroke(new BasicStroke(1.4f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
                    10f, new float[]{5f, 3f}, 0f));
            if (n == 1) g.fillOval(xs[0]-2, yS[0]-2, 4, 4);
            else        g.drawPolyline(xs, yS, n);
            g.setStroke(new BasicStroke(1f));

            g.setClip(oldClip);

            // Legend (top-right inside plot)
            drawLineScanLegend(g, px + pw - 8, py + 8);
        }

        private void drawLineScanLegend(Graphics2D g, int rightX, int topY) {
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            int lineW = 14, lineH = 12, y = topY;

            g.setStroke(new BasicStroke(1.8f));
            g.setColor(LS_G_CLR);
            g.drawLine(rightX - 60, y + 6, rightX - 60 + lineW, y + 6);
            g.setColor(AXIS_FG);
            g.drawString("G", rightX - 43, y + 9);
            y += lineH;

            g.setStroke(new BasicStroke(1.8f));
            g.setColor(LS_H_CLR);
            g.drawLine(rightX - 60, y + 6, rightX - 60 + lineW, y + 6);
            g.setColor(AXIS_FG);
            g.drawString("H", rightX - 43, y + 9);
            y += lineH;

            g.setStroke(new BasicStroke(1.4f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
                    10f, new float[]{5f, 3f}, 0f));
            g.setColor(LS_S_CLR);
            g.drawLine(rightX - 60, y + 6, rightX - 60 + lineW, y + 6);
            g.setColor(AXIS_FG);
            g.drawString("S (→)", rightX - 43, y + 9);
            g.setStroke(new BasicStroke(1f));
        }

        private void drawSeries(Graphics2D g, List<double[]> data, Color color,
                                int px, int py, int pw, int ph,
                                double xMin, double xMax, double yMin, double yMax) {
            if (data.isEmpty()) return;
            g.setColor(color);
            g.setStroke(new BasicStroke(1.5f));
            if (data.size() == 1) {
                int sx = toX(data.get(0)[0], xMin, xMax, px, pw);
                int sy = toY(data.get(0)[1], yMin, yMax, py, ph);
                g.fillOval(sx - 2, sy - 2, 4, 4);
                return;
            }
            int[] xs = new int[data.size()], ys = new int[data.size()];
            for (int i = 0; i < data.size(); i++) {
                xs[i] = toX(data.get(i)[0], xMin, xMax, px, pw);
                ys[i] = toY(data.get(i)[1], yMin, yMax, py, ph);
            }
            g.drawPolyline(xs, ys, data.size());
        }

        private void drawSeriesOnSecondaryAxis(Graphics2D g, List<double[]> data, Color color,
                                              int px, int py, int pw, int ph,
                                              double xMin, double xMax, double yMin, double yMax) {
            if (data.isEmpty()) return;
            g.setColor(color);
            g.setStroke(new BasicStroke(1.5f));
            if (data.size() == 1) {
                int sx = toX(data.get(0)[0], xMin, xMax, px, pw);
                int sy = toYR(data.get(0)[1], yMin, yMax, py, ph);
                g.fillOval(sx - 2, sy - 2, 4, 4);
                return;
            }
            int[] xs = new int[data.size()], ys = new int[data.size()];
            for (int i = 0; i < data.size(); i++) {
                xs[i] = toX(data.get(i)[0], xMin, xMax, px, pw);
                ys[i] = toYR(data.get(i)[1], yMin, yMax, py, ph);
            }
            g.drawPolyline(xs, ys, data.size());
        }

        private void drawLegend(Graphics2D g, int rightX, int topY) {
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 9));
            g.setStroke(new BasicStroke(1.5f));
            int lineW = 14;
            int y = topY;
            int lineH = 12;

            if (mode == Mode.MCS) {
                g.setColor(EQUIL_CLR);
                g.drawLine(rightX - 100, y + 6, rightX - 100 + lineW, y + 6);
                g.setColor(AXIS_FG);
                g.drawString("Equilibration", rightX - 83, y + 9);
                y += lineH;
                g.setColor(AVG_CLR);
                g.drawLine(rightX - 100, y + 6, rightX - 100 + lineW, y + 6);
                g.setColor(AXIS_FG);
                g.drawString("Averaging", rightX - 83, y + 9);
                y += lineH;
            } else if (mode == Mode.CVM) {
                g.setColor(CVM_CLR);
                g.drawLine(rightX - 80, y + 6, rightX - 80 + lineW, y + 6);
                g.setColor(AXIS_FG);
                g.drawString("|\u2207G| (log)", rightX - 63, y + 9);
                y += lineH;
                // CVM secondary axis: Gibbs energy only
                g.setColor(CVM_CLR);
                g.drawLine(rightX - 80, y + 6, rightX - 80 + lineW, y + 6);
                g.setColor(AXIS_FG);
                g.drawString("G", rightX - 63, y + 9);
            }
        }

        // ── scalar result strip ───────────────────────────────────────────────

        private void drawStrip(Graphics2D g, int x, int y, int w, int h) {
            g.setColor(STRIP_BG);
            g.fillRect(x, y, w, h);
            g.setColor(BORDER);
            g.drawLine(x, y, x + w, y);

            if (lastResult == null) {
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
                g.setColor(LABEL_DIM);
                String msg = (textResult == null || textResult.isBlank())
                        ? "— awaiting result —"
                        : "— showing loaded Type-1b result —";
                g.drawString(msg, x + 16, y + h / 2 + 4);
                return;
            }

            int colW = (w - 24) / 6;
            int x0   = x + 12;
            int row1 = y + 16, row2 = y + 36;

            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            g.setColor(LABEL_DIM);
            g.drawString("T (K)",       x0,               row1);
            g.drawString("x\u2082",     x0 + colW,        row1);
            g.drawString("G  (J/mol)",  x0 + 2 * colW,    row1);
            g.drawString("H  (J/mol)",  x0 + 3 * colW,    row1);
            g.drawString("Cv",          x0 + 4 * colW,    row1);
            g.drawString("σH",          x0 + 5 * colW,    row1);

            g.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
            g.setColor(VALUE_FG);
            g.drawString(String.format("%.1f",  lastResult.temperature),  x0,             row2);
            g.drawString(lastResult.composition.length > 1
                    ? String.format("%.4f", lastResult.composition[1]) : "\u2014",
                    x0 + colW,  row2);
            g.drawString(lastResult.isFreeEnergyValid()
                    ? String.format("%.4f", lastResult.gibbsEnergy) : "\u2014",
                    x0 + 2 * colW, row2);
            g.drawString(String.format("%.4f", lastResult.enthalpy),    x0 + 3 * colW, row2);
            g.drawString(!Double.isNaN(lastResult.heatCapacity)
                    ? String.format("%.5f", lastResult.heatCapacity) : "\u2014",
                    x0 + 4 * colW, row2);
            g.drawString(!Double.isNaN(lastResult.stdEnthalpy)
                    ? "±" + String.format("%.5f", lastResult.stdEnthalpy) : "\u2014",
                    x0 + 5 * colW, row2);
        }

        // ── helpers ───────────────────────────────────────────────────────────

        private static int toX(double v, double vMin, double vMax, int px, int pw) {
            return px + (int)((v - vMin) / (vMax - vMin) * pw);
        }

        private static int toY(double v, double vMin, double vMax, int py, int ph) {
            return py + ph - (int)((v - vMin) / (vMax - vMin) * ph);
        }

        private static int toYR(double v, double vMin, double vMax, int py, int ph) {
            return py + ph - (int)((v - vMin) / (vMax - vMin) * ph);
        }

        private static String formatAxisVal(double v) {
            double abs = Math.abs(v);
            if (abs == 0) return "0";
            if (abs >= 100 || (abs > 0 && abs < 0.01)) return String.format("%.2e", v);
            return String.format("%.3f", v);
        }
    }
}
