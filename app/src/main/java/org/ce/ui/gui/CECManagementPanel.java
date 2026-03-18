package org.ce.ui.gui;

import org.ce.domain.hamiltonian.CECEntry;
import org.ce.domain.hamiltonian.CECTerm;
import org.ce.workflow.cec.CECManagementWorkflow;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;

/**
 * Panel for scaffolding and inspecting the Hamiltonian (ECI) database.
 */
public class CECManagementPanel extends JPanel {

    private final CECManagementWorkflow cecWorkflow;

    private final JTextField clusterIdField     = new JTextField("BCC_A2_T_bin", 16);
    private final JTextField hamiltonianIdField  = new JTextField("A-B_BCC_A2_T", 16);
    private final JTextField elementsField       = new JTextField("A-B", 12);
    private final JTextField structPhaseField    = new JTextField("A2", 10);
    private final JTextField modelField          = new JTextField("T", 8);
    private final JTextArea  statusArea          = new JTextArea(4, 60);
    private final DefaultTableModel tableModel   = new DefaultTableModel(
            new Object[]{"Name", "a (J/mol)", "b (J/mol/K)"}, 0);
    private final JTable cecTable = new JTable(tableModel);

    public CECManagementPanel(CECManagementWorkflow cecWorkflow) {
        this.cecWorkflow = cecWorkflow;

        setLayout(new BorderLayout(8, 8));
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        add(buildForm(), BorderLayout.NORTH);
        add(buildTableSection(), BorderLayout.CENTER);
        add(buildStatus(), BorderLayout.SOUTH);
    }

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBorder(BorderFactory.createTitledBorder("Hamiltonian Database"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(4, 6, 4, 4);

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill = GridBagConstraints.HORIZONTAL;
        fc.insets = new Insets(4, 0, 4, 10);

        // Row 0: Cluster ID + Hamiltonian ID
        int col = 0;
        lc.gridx = col; lc.gridy = 0; fc.gridx = col + 1; fc.gridy = 0; col += 2;
        form.add(new JLabel("Cluster ID:"), lc);
        form.add(clusterIdField, fc);

        lc.gridx = col; fc.gridx = col + 1; col += 2;
        form.add(new JLabel("Hamiltonian ID:"), lc);
        form.add(hamiltonianIdField, fc);

        // Row 1: Elements + Structure-phase + Model
        col = 0;
        lc.gridx = col; lc.gridy = 1; fc.gridx = col + 1; fc.gridy = 1; col += 2;
        form.add(new JLabel("Elements:"), lc);
        form.add(elementsField, fc);

        lc.gridx = col; fc.gridx = col + 1; col += 2;
        form.add(new JLabel("Structure-phase:"), lc);
        form.add(structPhaseField, fc);

        lc.gridx = col; fc.gridx = col + 1;
        form.add(new JLabel("Model:"), lc);
        form.add(modelField, fc);

        // Row 2: buttons
        JPanel btnPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        JButton scaffoldBtn = new JButton("Scaffold Empty Hamiltonian");
        JButton loadBtn     = new JButton("Load Hamiltonian");
        scaffoldBtn.addActionListener(e -> scaffoldCEC());
        loadBtn.addActionListener(e -> loadCEC());
        btnPanel.add(scaffoldBtn);
        btnPanel.add(loadBtn);

        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 0; bc.gridy = 2;
        bc.gridwidth = 6;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(8, 6, 4, 6);
        form.add(btnPanel, bc);

        return form;
    }

    private JScrollPane buildTableSection() {
        cecTable.setFillsViewportHeight(true);
        cecTable.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        JScrollPane scroll = new JScrollPane(cecTable);
        scroll.setBorder(BorderFactory.createTitledBorder("ECI Terms"));
        return scroll;
    }

    private JPanel buildStatus() {
        statusArea.setEditable(false);
        statusArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 11));
        JScrollPane scroll = new JScrollPane(statusArea);
        scroll.setPreferredSize(new Dimension(0, 90));
        JPanel panel = new JPanel(new BorderLayout());
        panel.setBorder(BorderFactory.createTitledBorder("Status"));
        panel.add(scroll);
        return panel;
    }

    private void scaffoldCEC() {
        String clusterId     = clusterIdField.getText().trim();
        String hamiltonianId = hamiltonianIdField.getText().trim();
        String elements      = elementsField.getText().trim();
        String structPhase   = structPhaseField.getText().trim();
        String model         = modelField.getText().trim();

        statusArea.setText("Scaffolding Hamiltonian for " + hamiltonianId + "...\n");

        SwingWorker<CECEntry, Void> worker = new SwingWorker<>() {
            @Override
            protected CECEntry doInBackground() throws Exception {
                return cecWorkflow.scaffoldFromClusterData(hamiltonianId, elements, structPhase, model);
            }

            @Override
            protected void done() {
                try {
                    CECEntry entry = get();
                    populateTable(entry);
                    statusArea.append("Scaffolded " + entry.ncf + " terms. Saved to ~/CEWorkbench/hamiltonians/"
                            + hamiltonianId + "/hamiltonian.json\n");
                    statusArea.append("Edit hamiltonian.json to enter real ECI values.\n");
                } catch (Exception ex) {
                    statusArea.append("Error: " + ex.getMessage() + "\n");
                }
            }
        };
        worker.execute();
    }

    private void loadCEC() {
        String clusterId     = clusterIdField.getText().trim();
        String hamiltonianId = hamiltonianIdField.getText().trim();
        statusArea.setText("Loading Hamiltonian " + hamiltonianId + "...\n");

        SwingWorker<CECEntry, Void> worker = new SwingWorker<>() {
            @Override
            protected CECEntry doInBackground() throws Exception {
                return cecWorkflow.loadAndValidateCEC(clusterId, hamiltonianId);
            }

            @Override
            protected void done() {
                try {
                    CECEntry entry = get();
                    populateTable(entry);
                    statusArea.append("Loaded " + entry.ncf + " terms.\n");
                    statusArea.append("Elements: " + entry.elements
                            + "  Phase: " + entry.structurePhase
                            + "  Units: " + entry.cecUnits + "\n");
                } catch (Exception ex) {
                    statusArea.append("Error: " + ex.getMessage() + "\n");
                }
            }
        };
        worker.execute();
    }

    private void populateTable(CECEntry entry) {
        tableModel.setRowCount(0);
        if (entry.cecTerms == null) return;
        for (CECTerm term : entry.cecTerms) {
            tableModel.addRow(new Object[]{term.name, term.a, term.b});
        }
    }
}
