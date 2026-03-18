package org.ce.ui.gui;

import org.ce.domain.result.ThermodynamicResult;

import javax.swing.*;
import java.awt.*;

/**
 * Panel that displays the latest thermodynamic calculation result.
 */
public class ResultsPanel extends JPanel {

    private final JLabel temperatureLabel = makeValueLabel();
    private final JLabel compositionLabel = makeValueLabel();
    private final JLabel gibbsLabel       = makeValueLabel();
    private final JLabel enthalpyLabel    = makeValueLabel();
    private final JLabel entropyLabel     = makeValueLabel();

    public ResultsPanel() {
        setLayout(new BorderLayout(8, 8));
        setBorder(BorderFactory.createEmptyBorder(20, 30, 20, 30));

        add(buildResultGrid(), BorderLayout.CENTER);
        showEmpty();
    }

    private JPanel buildResultGrid() {
        JPanel grid = new JPanel(new GridBagLayout());
        grid.setBorder(BorderFactory.createTitledBorder("Equilibrium State"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(8, 12, 8, 12);

        GridBagConstraints vc = new GridBagConstraints();
        vc.anchor = GridBagConstraints.WEST;
        vc.fill   = GridBagConstraints.HORIZONTAL;
        vc.weightx = 1.0;
        vc.insets = new Insets(8, 0, 8, 12);

        addRow(grid, lc, vc, 0, "Temperature (K):",   temperatureLabel);
        addRow(grid, lc, vc, 1, "Composition x_B:",   compositionLabel);
        addRow(grid, lc, vc, 2, "Gibbs energy (J/mol):", gibbsLabel);
        addRow(grid, lc, vc, 3, "Enthalpy (J/mol):",  enthalpyLabel);
        addRow(grid, lc, vc, 4, "Entropy (J/mol/K):", entropyLabel);

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

    /** Populate all fields from a ThermodynamicResult. */
    public void showResult(ThermodynamicResult result) {
        temperatureLabel.setText(String.format("%.1f", result.temperature));
        compositionLabel.setText(result.composition.length > 1
                ? String.format("%.4f", result.composition[1]) : "—");
        gibbsLabel.setText(String.format("%.6f", result.gibbsEnergy));
        enthalpyLabel.setText(String.format("%.6f", result.enthalpy));

        double G = result.gibbsEnergy;
        double H = result.enthalpy;
        double T = result.temperature;
        double S = (T > 0) ? (H - G) / T : 0.0;
        entropyLabel.setText(String.format("%.6f", S));
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
