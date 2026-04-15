package org.ce.ui.gui;

import org.ce.model.ThermodynamicResult;

import javax.swing.*;
import java.awt.*;

/**
 * Panel that displays the latest thermodynamic calculation result.
 *
 * <p>Shows the system identity (from {@link WorkbenchContext}) at the top so
 * results are always associated with the system they belong to.</p>
 */
public class ResultsPanel extends JPanel {

    private final WorkbenchContext context;

    private final JLabel systemLabel      = new JLabel("— no system —");
    private final JLabel temperatureLabel = makeValueLabel();
    private final JLabel compositionLabel = makeValueLabel();
    private final JLabel gibbsLabel       = makeValueLabel();
    private final JLabel enthalpyLabel    = makeValueLabel();
    private final JLabel entropyLabel     = makeValueLabel();

    public ResultsPanel(WorkbenchContext context) {
        this.context = context;

        setLayout(new BorderLayout(8, 8));
        setBorder(BorderFactory.createEmptyBorder(20, 30, 20, 30));

        add(buildSystemBanner(), BorderLayout.NORTH);
        add(buildResultGrid(),   BorderLayout.CENTER);

        context.addChangeListener(this::refreshSystemBanner);
        refreshSystemBanner();

        showEmpty();
    }

    // =========================================================================
    // System banner
    // =========================================================================

    private JPanel buildSystemBanner() {
        JPanel banner = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 4));
        banner.setOpaque(false);

        JLabel titleLabel = new JLabel("System:  ");
        titleLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 12));
        titleLabel.setForeground(new Color(0x3A5070));

        systemLabel.setFont(new Font(Font.MONOSPACED, Font.BOLD, 13));
        systemLabel.setForeground(new Color(0x1A6020));

        banner.add(titleLabel);
        banner.add(systemLabel);
        return banner;
    }

    private void refreshSystemBanner() {
        if (context.hasSystem()) {
            var sys = context.getSystem();
            systemLabel.setText(sys.elements + "  /  " + sys.structure + "  /  " + sys.model);
            systemLabel.setForeground(new Color(0x1A6020));
        } else {
            systemLabel.setText("— no system selected —");
            systemLabel.setForeground(new Color(0x888888));
        }
    }

    // =========================================================================
    // Result grid
    // =========================================================================

    private JPanel buildResultGrid() {
        JPanel grid = new JPanel(new GridBagLayout());
        grid.setBorder(BorderFactory.createTitledBorder("Equilibrium State"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(8, 12, 8, 12);

        GridBagConstraints vc = new GridBagConstraints();
        vc.anchor  = GridBagConstraints.WEST;
        vc.fill    = GridBagConstraints.HORIZONTAL;
        vc.weightx = 1.0;
        vc.insets  = new Insets(8, 0, 8, 12);

        addRow(grid, lc, vc, 0, "Temperature (K):",    temperatureLabel);
        addRow(grid, lc, vc, 1, "Composition x_B:",    compositionLabel);
        addRow(grid, lc, vc, 2, "Gibbs energy (J/mol):", gibbsLabel);
        addRow(grid, lc, vc, 3, "Enthalpy (J/mol):",   enthalpyLabel);
        addRow(grid, lc, vc, 4, "Entropy (J/mol/K):",  entropyLabel);

        return grid;
    }

    private void addRow(JPanel panel, GridBagConstraints lc, GridBagConstraints vc,
                        int row, String labelText, JLabel valueLabel) {
        lc.gridx = 0; lc.gridy = row;
        vc.gridx = 1; vc.gridy = row;
        JLabel label = new JLabel(labelText);
        label.setFont(label.getFont().deriveFont(Font.BOLD));
        panel.add(label, lc);
        panel.add(valueLabel, vc);
    }

    private static JLabel makeValueLabel() {
        JLabel label = new JLabel("—");
        label.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 14));
        return label;
    }

    // =========================================================================
    // Public API
    // =========================================================================

    /** Populate all fields from a ThermodynamicResult. */
    public void showResult(ThermodynamicResult result) {
        temperatureLabel.setText(String.format("%.1f", result.temperature));
        compositionLabel.setText(result.composition.length > 1
                ? String.format("%.4f", result.composition[1]) : "—");
        gibbsLabel.setText(String.format("%.6f", result.gibbsEnergy));
        enthalpyLabel.setText(String.format("%.6f", result.enthalpy));
        entropyLabel.setText(Double.isNaN(result.entropy) ? "—" : String.format("%.6f", result.entropy));
        refreshSystemBanner();
    }

    /** Reset all fields to placeholder dashes. */
    public void showEmpty() {
        temperatureLabel.setText("—");
        compositionLabel.setText("—");
        gibbsLabel.setText("—");
        enthalpyLabel.setText("—");
        entropyLabel.setText("—");
    }
}
