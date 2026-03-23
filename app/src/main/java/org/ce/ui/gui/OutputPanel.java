package org.ce.ui.gui;

import org.ce.domain.engine.ProgressEvent;
import org.ce.domain.result.ThermodynamicResult;
import org.ce.storage.SystemId;

import javax.swing.*;
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
    private final JTextArea        logArea = new JTextArea();

    public OutputPanel(WorkbenchContext context) {
        setLayout(new BorderLayout());
        setBackground(BG);

        chartPanel = new ResultChartPanel(context);

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
        panel.add(chartPanel, BorderLayout.CENTER);
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
        SwingUtilities.invokeLater(() -> chartPanel.showResult(result));
    }

    /** Resets the chart. Thread-safe. */
    public void clearResult() {
        SwingUtilities.invokeLater(() -> chartPanel.clear());
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

        private enum Mode { IDLE, MCS, CVM }

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
            lastResult = null;
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
            lastResult = result;
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

            if (mode == Mode.IDLE) {
                g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
                g.setColor(LABEL_DIM);
                FontMetrics fm = g.getFontMetrics();
                String msg = "Run a calculation to see the convergence chart";
                g.drawString(msg, px + (pw - fm.stringWidth(msg)) / 2, py + ph / 2 + 4);
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
                g.drawString("— awaiting result —", x + 16, y + h / 2 + 4);
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
            g.drawString(String.format("%.4f", lastResult.gibbsEnergy), x0 + 2 * colW, row2);
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
