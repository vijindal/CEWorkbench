package org.ce.ui.gui;

import org.ce.storage.SystemId;

import javax.swing.*;
import java.awt.*;

/**
 * Slim title bar at the top of the window (VS Code dark style).
 *
 * <p>Shows the application name on the left and the active system identity
 * (elements · structure · model) on the right. Subscribes to
 * {@link WorkbenchContext} and updates automatically when the system changes.</p>
 *
 * <p>Visual spec:</p>
 * <ul>
 *   <li>Background:       {@code #323233}</li>
 *   <li>App name:         white, bold</li>
 *   <li>System identity:  {@code #9CDCFE} (VS Code light-blue) / {@code #555555} when absent</li>
 *   <li>Height:           28 px</li>
 * </ul>
 */
public class HeaderBar extends JPanel {

    private static final Color BG          = new Color(0x1A1A1A);
    private static final Color TEXT_TITLE  = new Color(0xCCCCCC);
    private static final Color TEXT_SYSTEM = new Color(0x9CDCFE);
    private static final Color TEXT_NONE   = new Color(0x555555);
    private static final Color SEPARATOR   = new Color(0x333333);

    private final WorkbenchContext context;
    private final JLabel systemLabel;

    public HeaderBar(WorkbenchContext context) {
        this.context = context;

        setBackground(BG);
        setPreferredSize(new Dimension(0, 28));
        setLayout(new BorderLayout(0, 0));
        setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(0x000000)));

        JLabel titleLabel = new JLabel("  CE WORKBENCH");
        titleLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 12));
        titleLabel.setForeground(TEXT_TITLE);
        titleLabel.setBorder(BorderFactory.createEmptyBorder(0, 8, 0, 0));
        add(titleLabel, BorderLayout.WEST);

        systemLabel = new JLabel();
        systemLabel.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
        systemLabel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 14));
        systemLabel.setHorizontalAlignment(SwingConstants.RIGHT);
        add(systemLabel, BorderLayout.EAST);

        context.addChangeListener(this::refreshSystem);
        refreshSystem();
    }

    private void refreshSystem() {
        if (context.hasSystem()) {
            SystemId s = context.getSystem();
            systemLabel.setText(s.elements + "  ·  " + s.structure + "  ·  " + s.model);
            systemLabel.setForeground(TEXT_SYSTEM);
        } else {
            systemLabel.setText("no system selected");
            systemLabel.setForeground(TEXT_NONE);
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        int sep = getWidth() - systemLabel.getPreferredSize().width - 28;
        g2.setColor(SEPARATOR);
        g2.drawLine(sep, 6, sep, getHeight() - 6);
    }
}
