package org.ce.ui.gui;

import javax.swing.*;
import java.awt.*;

/**
 * Thin footer strip that shows one-line operation status messages.
 *
 * <p>Panels post updates via the {@code postStatus} method exposed by
 * {@link MainWindow} rather than writing directly here, keeping them
 * decoupled from this component.</p>
 */
public class StatusBar extends JPanel {

    private static final Color BG  = new Color(0x1A1A1A);
    private static final Color FG  = new Color(0x858585);
    private static final Color TOP = new Color(0x000000);

    private final JLabel label;

    public StatusBar() {
        setBackground(BG);
        setBorder(BorderFactory.createMatteBorder(1, 0, 0, 0, TOP));
        setPreferredSize(new Dimension(0, 22));
        setLayout(new BorderLayout());

        label = new JLabel("  Ready");
        label.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        label.setForeground(FG);
        add(label, BorderLayout.WEST);
    }

    /** Updates the status bar text (must be called on the Event Dispatch Thread). */
    public void post(String message) {
        label.setText("  " + message);
    }
}
