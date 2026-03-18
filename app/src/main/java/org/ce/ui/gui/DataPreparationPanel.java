package org.ce.ui.gui;

import org.ce.domain.cluster.AllClusterData;
import org.ce.domain.cluster.CMatrixBuilder;
import org.ce.domain.cluster.CMatrixResult;
import org.ce.domain.cluster.Cluster;
import org.ce.domain.cluster.Vector3D;
import org.ce.storage.ClusterDataStore;
import org.ce.storage.InputLoader;
import org.ce.workflow.ClusterIdentificationRequest;
import org.ce.workflow.ClusterIdentificationWorkflow;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * Panel for running Type-1 cluster identification and saving AllClusterData.
 */
public class DataPreparationPanel extends JPanel {

    private final ClusterDataStore clusterStore;

    private final JTextField systemIdField   = new JTextField("A2-demo", 20);
    private final JTextField clusterFileField = new JTextField("clus/A2-T.txt", 28);
    private final JTextField symGroupField    = new JTextField("A2-SG", 20);
    private final JSpinner   numCompSpinner   = new JSpinner(new SpinnerNumberModel(5, 2, 20, 1));
    private final JTextArea  logArea          = new JTextArea(14, 60);

    public DataPreparationPanel(ClusterDataStore clusterStore) {
        this.clusterStore = clusterStore;

        setLayout(new BorderLayout(8, 8));
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        add(buildForm(), BorderLayout.NORTH);
        add(buildLog(), BorderLayout.CENTER);
    }

    private JPanel buildForm() {
        JPanel form = new JPanel(new GridBagLayout());
        form.setBorder(BorderFactory.createTitledBorder("Cluster Identification"));

        GridBagConstraints lc = new GridBagConstraints();
        lc.anchor = GridBagConstraints.WEST;
        lc.insets = new Insets(4, 6, 4, 4);

        GridBagConstraints fc = new GridBagConstraints();
        fc.fill = GridBagConstraints.HORIZONTAL;
        fc.weightx = 1.0;
        fc.insets = new Insets(4, 0, 4, 6);

        int row = 0;

        lc.gridx = 0; lc.gridy = row; fc.gridx = 1; fc.gridy = row++;
        form.add(new JLabel("System ID:"), lc);
        form.add(systemIdField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Cluster file (classpath):"), lc);
        form.add(clusterFileField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Symmetry group:"), lc);
        form.add(symGroupField, fc);

        lc.gridy = row; fc.gridy = row++;
        form.add(new JLabel("Num components:"), lc);
        form.add(numCompSpinner, fc);

        JButton runBtn = new JButton("Run Identification & Save");
        runBtn.addActionListener(e -> runIdentification());

        GridBagConstraints bc = new GridBagConstraints();
        bc.gridx = 1; bc.gridy = row;
        bc.anchor = GridBagConstraints.WEST;
        bc.insets = new Insets(8, 0, 4, 6);
        form.add(runBtn, bc);

        return form;
    }

    private JScrollPane buildLog() {
        logArea.setEditable(false);
        logArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        logArea.setText("Output will appear here after running identification.\n");
        return new JScrollPane(logArea);
    }

    private void runIdentification() {
        String systemId   = systemIdField.getText().trim();
        String clusterFile = clusterFileField.getText().trim();
        String symGroup   = symGroupField.getText().trim();
        int    numComp    = (int) numCompSpinner.getValue();

        logArea.setText("Running cluster identification for system: " + systemId + "\n");

        SwingWorker<AllClusterData, String> worker = new SwingWorker<>() {
            @Override
            protected AllClusterData doInBackground() throws Exception {
                publish("Building configuration...");

                ClusterIdentificationRequest config = ClusterIdentificationRequest.builder()
                        .disorderedClusterFile(clusterFile)
                        .orderedClusterFile(clusterFile)
                        .disorderedSymmetryGroup(symGroup)
                        .orderedSymmetryGroup(symGroup)
                        .transformationMatrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}})
                        .translationVector(new Vector3D(0,0,0))
                        .numComponents(numComp)
                        .build();

                publish("Stage 1-2: Cluster + CF identification...");
                AllClusterData partial = ClusterIdentificationWorkflow.identify(config);

                publish("Stage 3: Building C-matrix...");
                List<Cluster> maxClusters = InputLoader.parseClusterFile(clusterFile);
                CMatrixResult cmatrix = CMatrixBuilder.build(
                        partial.getDisorderedClusterResult(),
                        partial.getDisorderedCFResult(),
                        maxClusters,
                        numComp
                );

                publish("Saving AllClusterData...");
                AllClusterData full = new AllClusterData(
                        partial.getDisorderedClusterResult(),
                        partial.getOrderedClusterResult(),
                        partial.getDisorderedCFResult(),
                        partial.getOrderedCFResult(),
                        cmatrix
                );
                clusterStore.save(systemId, full);
                return full;
            }

            @Override
            protected void process(java.util.List<String> chunks) {
                for (String msg : chunks) {
                    logArea.append(msg + "\n");
                }
            }

            @Override
            protected void done() {
                try {
                    AllClusterData result = get();
                    logArea.append("\nDone! " + result.getSummary());
                    logArea.append("\nSaved to ~/CEWorkbench/cluster-data/" + systemId + "/cluster_data.json\n");
                } catch (Exception ex) {
                    logArea.append("\nError: " + ex.getMessage() + "\n");
                }
            }
        };

        worker.execute();
    }
}
