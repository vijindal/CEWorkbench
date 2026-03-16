package org.ce.domain.engine;

import org.ce.domain.cluster.ClusterData;

/**
 * Input parameters for CVM thermodynamic calculations.
 */
public class CVMInput {

    private final ClusterData clusterData;
    private final double[] eci;
    private final double temperature;
    private final double[] composition;

    private final double tolerance;

    public CVMInput(
            ClusterData clusterData,
            double[] eci,
            double temperature,
            double[] composition,
            double tolerance) {

        this.clusterData = clusterData;
        this.eci = eci.clone();
        this.temperature = temperature;
        this.composition = composition.clone();
        this.tolerance = tolerance;
    }

    public ClusterData getClusterData() {
        return clusterData;
    }

    public double[] getEci() {
        return eci.clone();
    }

    public double getTemperature() {
        return temperature;
    }

    public double[] getComposition() {
        return composition.clone();
    }

    public double getTolerance() {
        return tolerance;
    }
}
