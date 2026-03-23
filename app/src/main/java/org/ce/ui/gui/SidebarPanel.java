package org.ce.ui.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.List;

/**
 * Sidebar navigation panel for the CE Workbench.
 *
 * <p>Displays four workflow steps as a vertical list. Each step shows a numbered
 * badge, a bold label, and a small subtitle. Clicking a step fires the navigation
 * callback supplied by {@link MainWindow}.</p>
 *
 * <p>Visual spec:</p>
 * <ul>
 *   <li>Background: {@code #2B3A55} (navy)</li>
 *   <li>Active left border: {@code #4A90D9} (blue, 3 px)</li>
 *   <li>Active text: white, bold</li>
 *   <li>Inactive text: {@code #A0B4CC}</li>
 *   <li>Hover background: {@code #374A6A}</li>
 * </ul>
 */
public class SidebarPanel extends JPanel {

    // ── colours ───────────────────────────────────────────────────────────────
    static final Color BG           = new Color(0x2B3A55);
    static final Color BG_HOVER     = new Color(0x374A6A);
    static final Color ACCENT       = new Color(0x4A90D9);
    static final Color TEXT_ACTIVE  = Color.WHITE;
    static final Color TEXT_INACTIVE = new Color(0xA0B4CC);
    static final Color BADGE_BG     = new Color(0x1A2840);

    // ── step definitions ──────────────────────────────────────────────────────
    private static final String[][] STEPS = {
        { "1a", "Data Prep",    "Cluster ID"       },
        { "1b", "Hamiltonian",  "ECI / CEC"        },
        { "2",  "Calculate",    "CVM / MCS"        },
        { "",   "Results",      "Equilibrium state" },
    };

    // ── state ─────────────────────────────────────────────────────────────────
    private final List<StepButton> buttons = new ArrayList<>();
    private int activeIndex = 0;

    public SidebarPanel(Runnable[] navigateCallbacks) {
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setBackground(BG);
        setPreferredSize(new Dimension(178, 0));
        setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1, new Color(0x1A2840)));

        for (int i = 0; i < STEPS.length; i++) {
            StepButton btn = new StepButton(i, STEPS[i][0], STEPS[i][1], STEPS[i][2],
                    navigateCallbacks[i]);
            buttons.add(btn);
            add(btn);
        }
        add(Box.createVerticalGlue());

        setActive(0);
    }

    /** Highlights the step at {@code index} as the current active step. */
    public void setActive(int index) {
        activeIndex = index;
        buttons.forEach(b -> b.setActive(b.index == activeIndex));
        revalidate();
        repaint();
    }

    /** Returns the index of the currently active step. */
    public int getActiveIndex() {
        return activeIndex;
    }

    // =========================================================================
    // Inner: StepButton
    // =========================================================================

    private static class StepButton extends JPanel {

        final int index;
        private final String badge;
        private final String label;
        private final String subtitle;
        private final Runnable onClick;

        private boolean active  = false;
        private boolean hovered = false;

        StepButton(int index, String badge, String label, String subtitle, Runnable onClick) {
            this.index    = index;
            this.badge    = badge;
            this.label    = label;
            this.subtitle = subtitle;
            this.onClick  = onClick;

            setOpaque(false);
            setPreferredSize(new Dimension(178, 62));
            setMaximumSize(new Dimension(Integer.MAX_VALUE, 62));
            setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

            addMouseListener(new MouseAdapter() {
                @Override public void mouseEntered(MouseEvent e) { hovered = true;  repaint(); }
                @Override public void mouseExited(MouseEvent e)  { hovered = false; repaint(); }
                @Override public void mouseClicked(MouseEvent e) { if (onClick != null) onClick.run(); }
            });
        }

        void setActive(boolean active) {
            this.active = active;
            repaint();
        }

        @Override
        protected void paintComponent(Graphics g0) {
            Graphics2D g = (Graphics2D) g0.create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int w = getWidth(), h = getHeight();

            // Background
            g.setColor(hovered && !active ? BG_HOVER : BG);
            g.fillRect(0, 0, w, h);

            // Active left accent bar
            if (active) {
                g.setColor(ACCENT);
                g.fillRect(0, 0, 3, h);
            }

            // Badge circle (left margin)
            int cx = 26, cy = h / 2;
            int r  = 14;
            g.setColor(active ? ACCENT : BADGE_BG);
            g.fillOval(cx - r, cy - r, r * 2, r * 2);

            if (!badge.isEmpty()) {
                g.setColor(active ? Color.WHITE : TEXT_INACTIVE);
                g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, badge.length() > 1 ? 9 : 11));
                FontMetrics fm = g.getFontMetrics();
                g.drawString(badge, cx - fm.stringWidth(badge) / 2,
                        cy + fm.getAscent() / 2 - 1);
            }

            // Step label
            int tx = cx + r + 10;
            g.setColor(active ? TEXT_ACTIVE : TEXT_INACTIVE);
            g.setFont(new Font(Font.SANS_SERIF, active ? Font.BOLD : Font.PLAIN, 13));
            g.drawString(label, tx, cy - 3);

            // Subtitle
            g.setColor(active ? new Color(0xCCDDEE) : new Color(0x6B829E));
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            g.drawString(subtitle, tx, cy + 12);

            g.dispose();
        }
    }
}
