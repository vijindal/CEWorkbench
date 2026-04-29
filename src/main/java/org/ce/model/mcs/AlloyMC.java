package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.Random;
import java.util.function.Consumer;
import java.util.logging.Logger;

/**
 * A simplified Monte Carlo engine for alloy equilibrium simulations.
 * Implements the Metropolis algorithm for canonical ensemble (atom swaps).
 *
 * This version manages its own Random number generator internally.
 */
public class AlloyMC {
    private static final Logger LOG = Logger.getLogger(AlloyMC.class.getName());

    private final MCSGeometry geo;
    private final LatticeConfig config;
    private final CECEntry cecEntry;
    private final Random rng;
    private final double R; // Gas constant
    private final int L;    // Supercell dimension
    private final int numComp;
    private final int ncf;

    private final Embeddings.DeltaScratch scratch;
    private double[] correlationFunctions; // v_cvcf (Correlation Functions in CVCF basis)

    private double temperature = Double.NaN;
    private double[] composition;
    private double[] eciCvcf;
    
    private long attempts = 0;
    private long accepted = 0;

    /**
     * Initializes the MC engine and builds the expensive geometry.
     *
     * @param session      The model session containing system parameters.
     * @param L            The supercell dimension.
     * @param progressSink Sink for geometry construction logs.
     */
    public AlloyMC(ModelSession session, int L, Consumer<String> progressSink) {
        this.L = L;
        this.cecEntry = session.cecEntry;
        this.rng = new Random(); // Internal initialization
        this.R = PhysicsConstants.R_GAS;

        // 1. Initiate geometry (Expensive - done once)
        LOG.info("AlloyMC: building geometry for L=" + L);
        this.geo = MCSGeometry.build(session, L, progressSink);
        
        // 2. Initialize lattice configuration based on geometry
        this.numComp = session.numComponents();
        this.config = new LatticeConfig(geo.nSites(), numComp);
        this.ncf = geo.getNcf();

        // 3. Pre-allocate scratch space for ΔE (O(1) allocation)
        int totalEmb = Embeddings.totalCfEmbeddingCount(geo.cfEmbeddings());
        this.scratch = new Embeddings.DeltaScratch(ncf, totalEmb);
        
        // 4. Initialize Correlation Functions state array
        // Layout: [uOrthNonPoint (ncf) | uPoint (numComp-1) | Empty (1)]
        this.correlationFunctions = new double[ncf + numComp];

        LOG.info(String.format("AlloyMC initialized: L=%d, sites=%d, ncf=%d", 
                L, config.getN(), ncf));
    }

    /**
     * Sets the seed for the internal random number generator to ensure reproducibility.
     */
    public void setSeed(long seed) {
        this.rng.setSeed(seed);
        LOG.info("AlloyMC: random seed set to " + seed);
    }

    // --- State Accessors ---

    public double getTemperature() { return temperature; }

    /**
     * Sets the simulation temperature and re-evaluates internal ECIs.
     */
    public void setTemperature(double T) {
        setTemperature(T, null);
    }

    public void setTemperature(double T, Consumer<String> progressSink) {
        if (Double.compare(T, this.temperature) == 0) return;
        this.temperature = T;
        if (cecEntry != null && geo.basis != null) {
            LOG.fine("AlloyMC: evaluating ECIs at T=" + T);
            this.eciCvcf = CECEvaluator.evaluate(cecEntry, T, geo.basis, "MCS", progressSink);
        } else {
            LOG.warning("AlloyMC: cannot evaluate ECIs (cecEntry or basis is null)");
            this.eciCvcf = new double[0];
        }
    }

    public double[] getComposition() {
        return composition != null ? composition.clone() : null;
    }

    /**
     * Sets the simulation composition.
     */
    public void setComposition(double[] x) {
        if (x == null || x.length != config.getNumComp()) {
            throw new IllegalArgumentException("Invalid composition array length");
        }
        this.composition = x.clone();
    }

    /**
     * Runs the MC simulation for the current state (T, x) and specified parameters.
     * 
     * @param nEquil       Number of equilibration sweeps.
     * @param nAvg         Number of averaging sweeps.
     */
    public void run(int nEquil, int nAvg) {
        if (Double.isNaN(temperature)) throw new IllegalStateException("Temperature not set");
        if (composition == null)       throw new IllegalStateException("Composition not set");

        LOG.info(String.format("Starting AlloyMC run: T=%.2f K, x=%s, nEquil=%d, nAvg=%d", 
                temperature, java.util.Arrays.toString(composition), nEquil, nAvg));

        // 1. Initialise lattice occupation based on composition
        config.randomise(composition, rng);

        // 2. Initial state energy
        double currentEnergy = calculateTotalEnergy();

        // 3. Equilibration phase
        for (int s = 0; s < nEquil; s++) {
            currentEnergy += runSweep();
        }

        resetCounters();

        // 4. Averaging phase
        for (int s = 0; s < nAvg; s++) {
            currentEnergy += runSweep();
            sampleProperties(currentEnergy);
        }

        LOG.info(String.format("AlloyMC run complete. Final energy: %.6f, Acceptance rate: %.2f%%", 
                currentEnergy, getAcceptanceRate() * 100));
    }

    /**
     * Executes one sweep (N trial moves).
     * @return Total change in energy for the sweep.
     */
    private double runSweep() {
        double sweepDeltaE = 0;
        int N = config.getN();
        for (int i = 0; i < N; i++) {
            sweepDeltaE += attemptSwap();
        }
        return sweepDeltaE;
    }

    /**
     * Attempts a single canonical (swap) move between two random sites.
     * @return Change in energy if the move is accepted, 0.0 otherwise.
     */
    private double attemptSwap() {
        int N = config.getN();
        int i = rng.nextInt(N);
        int j = rng.nextInt(N);

        int occI = config.getOccupation(i);
        int occJ = config.getOccupation(j);

        if (occI == occJ) return 0.0;

        attempts++;

        double dE = calculateDeltaE(i, j, occI, occJ);

        // Metropolis acceptance criterion
        boolean accept = (dE <= 0) || (rng.nextDouble() < Math.exp(-dE / (R * temperature)));

        if (accept) {
            config.setOccupation(i, occJ);
            config.setOccupation(j, occI);
            accepted++;
            return dE;
        }

        return 0.0;
    }

    /**
     * Calculates the total energy of the current configuration.
     * E = N_sites * sum( ECI_l * CF_l )
     */
    public double calculateTotalEnergy() {
        if (eciCvcf == null) return 0.0;
        
        // 1. Update CFs to match the current configuration
        updateCorrelationFunctions();
        
        // 2. Get the CVCF-basis CFs for dot product with ECIs (Delegated to Geometry)
        double[] vFull = geo.getCvcfCorrelationFunctions(correlationFunctions, composition);
        
        // 3. Compute Dot Product: Energy per site
        double energyPerSite = 0.0;
        for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
            energyPerSite += eciCvcf[l] * vFull[l];
        }
        
        return energyPerSite * getNSites();
    }

    /**
     * Recalculates the Correlation Functions for the current lattice state.
     * Populates the internal state array with [uOrthNonPoint | uPoint | Empty].
     */
    public void updateCorrelationFunctions() {
        if (geo.cfEmbeddings() == null) return;

        // 1. Measure orthogonal non-point CFs (uOrth)
        double[] uOrthNonPoint = Embeddings.measureCVsFromConfig(
                config,
                geo.cfEmbeddings(),
                geo.getFlatBasisMatrix(),
                ncf,
                numComp);
        
        // 2. Populate the state array: [uOrthNonPoint | uPoint | Empty]
        System.arraycopy(uOrthNonPoint, 0, correlationFunctions, 0, ncf);
        
        // Point CFs: Average of basis functions over the lattice
        double[] x = config.composition();
        double[] basisSeq = org.ce.model.cluster.ClusterMath.buildBasis(numComp);
        for (int k = 1; k < numComp; k++) {
            double phiK = 0.0;
            for (int s = 0; s < numComp; s++) {
                phiK += x[s] * Math.pow(basisSeq[s], k);
            }
            correlationFunctions[ncf + k - 1] = phiK;
        }
        
        // Empty cluster is always 1.0
        correlationFunctions[ncf + numComp - 1] = 1.0;
    }
    
    public double[] getCvcfCorrelationFunctions() {
        return geo.getCvcfCorrelationFunctions(correlationFunctions, composition);
    }

    /**
     * Calculates the change in energy for swapping occupations at sites i and j.
     */
    private double calculateDeltaE(int i, int j, int oldOccI, int oldOccJ) {
        // TODO: Implement incremental update logic
        return 0.0;
    }

    /**
     * Samples thermodynamic properties during the averaging phase.
     */
    private void sampleProperties(double energy) {
        // TODO: Implement property accumulation
    }

    private void resetCounters() {
        attempts = 0;
        accepted = 0;
    }

    public double getAcceptanceRate() {
        return attempts == 0 ? 0.0 : (double) accepted / attempts;
    }

    public double[] getEciCvcf() { return eciCvcf; }
    public LatticeConfig getConfig() { return config; }
    public MCSGeometry getGeo() { return geo; }
    public Random getRng() { return rng; }

    public long getAttempts() { return attempts; }
    public long getAccepted() { return accepted; }
    public CECEntry getCecEntry() { return cecEntry; }
    public double[] getCorrelationFunctions() { return correlationFunctions; }
    
    public void setCorrelationFunctions(double[] correlationFunctions) {
        if (correlationFunctions.length != this.correlationFunctions.length) {
            throw new IllegalArgumentException("CF vector size mismatch. Expected: " + this.correlationFunctions.length);
        }
        this.correlationFunctions = correlationFunctions.clone();
    }
    public int getL() { return L; }
    public int getNSites() { return geo != null ? geo.nSites() : 0; }
    public int getNumComp() { return numComp; }
    public int getNcf() { return ncf; }
    public double getR_Gas() { return R; }
}
