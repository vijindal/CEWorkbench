package org.ce.domain.engine;

import org.ce.domain.cluster.ClusterData;

/**
 * Input parameters for Monte Carlo thermodynamic calculations.
 */
public class MCSInput {

    private final ClusterData clusterData;
    private final double[] eci;
    private final double temperature;
    private final double[] composition;

    private final int supercellSize;
    private final int equilibrationSteps;
    private final int averagingSteps;

    public MCSInput(
            ClusterData clusterData,
            double[] eci,
            double temperature,
            double[] composition,
            int supercellSize,
            int equilibrationSteps,
            int averagingSteps) {

        this.clusterData = clusterData;
        this.eci = eci.clone();
        this.temperature = temperature;
        this.composition = composition.clone();
        this.supercellSize = supercellSize;
        this.equilibrationSteps = equilibrationSteps;
        this.averagingSteps = averagingSteps;
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

    public int getSupercellSize() {
        return supercellSize;
    }

    public int getEquilibrationSteps() {
        return equilibrationSteps;
    }

    public int getAveragingSteps() {
        return averagingSteps;
    }
}
