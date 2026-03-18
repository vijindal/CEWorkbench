package org.ce.ui.gui;

import org.ce.domain.result.ThermodynamicResult;
import org.ce.workflow.CalculationService;

import javax.swing.*;
import java.awt.*;

/**
 * Panel for running single-point thermodynamic calculations.
 */
public class CalculationPanel extends JPanel {

    private final CalculationService service;
    private final ResultsPanel       resultsPanel;
    private final JTabbedPane        tabs;

    private final JTextField clusterIdField     = new JTextField("BCC_A2_T_bin", 18);
    private final JTextField hamiltonianIdField  = new JTextField("A-B_BCC_A2_T", 18);
    private final JSpinner   temperatureField    = new JSpinner(new SpinnerNumberModel(1000.0, 1.0, 10000.0, 100.0));
    private final JSpinner   xBField             = new JSpinner(new SpinnerNumberModel(0.5, 0.0, 1.0, 0.05));
    private final JComboBox<String> engineBox    = new JComboBox<>(new String[]{"CVM", "MCS"});
    private final JTextArea  logArea             = new JTextArea(8, 60);
    private final JButton    calcButton          = new JButton("Calculate");

    public CalculationPanel(CalculationService service, ResultsPanel resultsPanel, JTabbedPane tabs) {
        this.service      = service;
        this.resultsPanel = resultsPanel;
        this.tabs         = tabs;

        setLayout(new BorderLayout(8, 8));
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        add(buildForm(), BorderLayout.NORTH);
        add(buildLog(), BorderLayout.CENTER);
    }

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBorder(BorderFactory.createTitledBorder("Calculation Parameters"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(4, 6, 4, 4);

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill = GridBagConstraints.HORIZONTAL;
        fc.insets = new Insets(4, 0, 4, 14);

        int row = 0;

        lc.gridx = 0; lc.gridy = row; fc.gridx = 1; fc.gridy = row++;
        form.add(new JLabel("Cluster ID:"), lc);
        form.add(clusterIdField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Hamiltonian ID:"), lc);
        form.add(hamiltonianIdField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Temperature (K):"), lc);
        form.add(temperatureField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Composition x_B:"), lc);
        form.add(xBField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Engine:"), lc);
        form.add(engineBox, fc);

        calcButton.addActionListener(e -> runCalculation());
        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 1; bc.gridy = row;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(8, 0, 4, 6);
        form.add(calcButton, bc);

        return form;
    }

    private JScrollPane buildLog() {
        logArea.setEditable(false);
        logArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        logArea.setText("Press Calculate to run a thermodynamic calculation.\n");
        JScrollPane scroll = new JScrollPane(logArea);
        scroll.setBorder(BorderFactory.createTitledBorder("Log"));
        return scroll;
    }

    private void runCalculation() {
        String clusterId     = clusterIdField.getText().trim();
        String hamiltonianId = hamiltonianIdField.getText().trim();
        double temperature   = ((Number) temperatureField.getValue()).doubleValue();
        double xB            = ((Number) xBField.getValue()).doubleValue();
        String engineType    = (String) engineBox.getSelectedItem();
        double[] composition = {1.0 - xB, xB};

        calcButton.setEnabled(false);
        logArea.setText("Running " + engineType + " calculation...\n");
        logArea.append("  Cluster     : " + clusterId + "\n");
        logArea.append("  Hamiltonian : " + hamiltonianId + "\n");
        logArea.append("  T = " + temperature + " K,  x_B = " + xB + "\n");

        SwingWorker<ThermodynamicResult, Void> worker = new SwingWorker<>() {
            @Override
            protected ThermodynamicResult doInBackground() throws Exception {
                return service.runSinglePoint(clusterId, hamiltonianId, temperature, composition, engineType);
            }

            @Override
            protected void done() {
                calcButton.setEnabled(true);
                try {
                    ThermodynamicResult result = get();
                    logArea.append("\nCalculation complete.\n");
                    logArea.append("  G = " + String.format("%.4f", result.gibbsEnergy) + " J/mol\n");
                    logArea.append("  H = " + String.format("%.4f", result.enthalpy) + " J/mol\n");
                    resultsPanel.showResult(result);
                    tabs.setSelectedIndex(3);
                } catch (Exception ex) {
                    logArea.append("\nError: " + ex.getMessage() + "\n");
                }
            }
        };
        worker.execute();
    }
}
