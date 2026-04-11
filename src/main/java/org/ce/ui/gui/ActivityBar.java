package org.ce.ui.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.List;

/**
 * VS Code-style thin activity bar (left strip).
 *
 * <p>Shows three step icons — Data Prep, Hamiltonian, Calculate. Clicking an icon
 * fires the corresponding navigation callback and highlights the active item with
 * a blue left accent bar.</p>
 *
 * <p>Visual spec:</p>
 * <ul>
 *   <li>Background:       {@code #333333}</li>
 *   <li>Active border:    {@code #007ACC} (2 px, left)</li>
 *   <li>Active icon:      white</li>
 *   <li>Inactive icon:    {@code #858585}</li>
 *   <li>Hover background: {@code #2A2D2E}</li>
 * </ul>
 */
public class ActivityBar extends JPanel {

    // ── colours ───────────────────────────────────────────────────────────────
    static final Color BG           = new Color(0x333333);
    static final Color ACTIVE_FG    = Color.WHITE;
    static final Color INACTIVE_FG  = new Color(0x858585);
    static final Color ACTIVE_BDR   = new Color(0x007ACC);
    static final Color HOVER_BG     = new Color(0x2A2D2E);

    // ── step definitions: { badge, short label } ──────────────────────────────
    private static final String[][] ITEMS = {
        { "1a", "Data Prep"  },
        { "1b", "Hamilt."    },
        { "TH", "Thermo"     },
    };

    private final List<ActivityItem> items = new ArrayList<>();
    private int activeIndex = 0;

    public ActivityBar(Runnable[] callbacks) {
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setBackground(BG);
        setPreferredSize(new Dimension(52, 0));
        setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1, new Color(0x000000)));

        for (int i = 0; i < ITEMS.length; i++) {
            ActivityItem item = new ActivityItem(i, ITEMS[i][0], ITEMS[i][1], callbacks[i]);
            items.add(item);
            add(item);
        }
        add(Box.createVerticalGlue());
        setActive(0);
    }

    /** Highlights item at {@code index}. */
    public void setActive(int index) {
        activeIndex = index;
        items.forEach(it -> it.setActive(it.index == activeIndex));
    }

    // =========================================================================
    // Inner: ActivityItem
    // =========================================================================

    private static class ActivityItem extends JPanel {

        final int    index;
        final String badge;
        final String label;
        final Runnable onClick;

        boolean active;
        boolean hovered;

        ActivityItem(int index, String badge, String label, Runnable onClick) {
            this.index   = index;
            this.badge   = badge;
            this.label   = label;
            this.onClick = onClick;

            setOpaque(false);
            setPreferredSize(new Dimension(52, 56));
            setMaximumSize(new Dimension(52, 56));
            setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

            addMouseListener(new MouseAdapter() {
                @Override public void mouseEntered(MouseEvent e) { hovered = true;  repaint(); }
                @Override public void mouseExited(MouseEvent e)  { hovered = false; repaint(); }
                @Override public void mouseClicked(MouseEvent e) { if (onClick != null) onClick.run(); }
            });
        }

        void setActive(boolean a) { this.active = a; repaint(); }

        @Override
        protected void paintComponent(Graphics g0) {
            Graphics2D g = (Graphics2D) g0.create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,      RenderingHints.VALUE_ANTIALIAS_ON);
            g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            int w = getWidth(), h = getHeight();

            // Background
            g.setColor(hovered && !active ? HOVER_BG : BG);
            g.fillRect(0, 0, w, h);

            // Active left accent (2 px)
            if (active) {
                g.setColor(ACTIVE_BDR);
                g.fillRect(0, h / 4, 2, h / 2);
            }

            Color fg = active ? ACTIVE_FG : INACTIVE_FG;

            // Badge — larger centred text
            g.setColor(fg);
            g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, badge.length() > 1 ? 11 : 14));
            FontMetrics fm = g.getFontMetrics();
            g.drawString(badge, (w - fm.stringWidth(badge)) / 2, h / 2 - 1);

            // Short label below badge
            g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 8));
            fm = g.getFontMetrics();
            String shortened = label.length() > 7 ? label.substring(0, 7) : label;
            g.drawString(shortened, (w - fm.stringWidth(shortened)) / 2, h / 2 + 11);

            g.dispose();
        }
    }
}
