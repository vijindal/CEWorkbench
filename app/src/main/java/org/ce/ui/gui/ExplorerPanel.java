package org.ce.ui.gui;

import javax.swing.*;
import java.awt.*;

/**
 * VS Code-style explorer/sidebar panel (second column).
 *
 * <p>Fixed width (~290 px), dark background. A title bar at the top shows the
 * name of the currently active section in small-caps style. The body uses
 * {@link CardLayout} to show exactly one parameter panel at a time.</p>
 *
 * <p>Visual spec:</p>
 * <ul>
 *   <li>Background:      {@code #252526}</li>
 *   <li>Title bar:       {@code #252526}, text {@code #BBBBBB}</li>
 *   <li>Right border:    {@code #3C3C3C}</li>
 * </ul>
 */
public class ExplorerPanel extends JPanel {

    private static final Color BG         = new Color(0x252526);
    private static final Color TITLE_FG   = new Color(0xBBBBBB);
    private static final Color BORDER_CLR = new Color(0x2D2D2D);

    private static final String[] TITLES = {
        "DATA PREPARATION",
        "HAMILTONIAN  /  ECI",
        "CALCULATION",
    };

    private static final String[] CARDS = { "DataPrep", "Hamiltonian", "Calculate" };

    private final CardLayout cardLayout = new CardLayout();
    private final JPanel     cardPanel  = new JPanel(cardLayout);
    private final JLabel     titleLabel;

    public ExplorerPanel() {
        setLayout(new BorderLayout());
        setBackground(BG);
        setPreferredSize(new Dimension(290, 0));
        setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1, BORDER_CLR));

        // ── title bar ─────────────────────────────────────────────────────────
        JPanel titleBar = new JPanel(new BorderLayout());
        titleBar.setBackground(BG);
        titleBar.setPreferredSize(new Dimension(0, 28));
        titleBar.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createMatteBorder(0, 0, 1, 0, BORDER_CLR),
                BorderFactory.createEmptyBorder(0, 12, 0, 0)));

        titleLabel = new JLabel(TITLES[0]);
        titleLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 10));
        titleLabel.setForeground(TITLE_FG);
        titleBar.add(titleLabel, BorderLayout.WEST);

        // ── card body ─────────────────────────────────────────────────────────
        cardPanel.setBackground(BG);

        add(titleBar,  BorderLayout.NORTH);
        add(cardPanel, BorderLayout.CENTER);
    }

    /** Registers a parameter panel at the given slot index (0=DataPrep, 1=Hamiltonian, 2=Calc). */
    public void addCard(JPanel panel, int index) {
        cardPanel.add(panel, CARDS[index]);
    }

    /** Shows the card at {@code index} and updates the title bar. */
    public void showCard(int index) {
        if (index < 0 || index >= CARDS.length) return;
        titleLabel.setText(TITLES[index]);
        cardLayout.show(cardPanel, CARDS[index]);
    }
}
