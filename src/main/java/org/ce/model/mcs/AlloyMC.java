package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.PhysicsConstants;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.hamiltonian.CECEvaluator;

import java.util.List;
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
    private final int nSites;
    private final int numComp;
    private final int ncf;

    private final Embeddings.DeltaScratch scratch;
    private double[] correlationFunctions;     // u (Orthogonal Basis)
    private double[] cvcfCorrelationFunctions; // v (CVCF Basis)

    private double temperature = Double.NaN;
    private double[] composition;
    private double[] eciCvcf;
    private int maxEmbPerCol;
    
    private long attempts = 0;
    private long accepted = 0;

    // --- Averaging Accumulators ---
    private double   sumEnergy = 0;
    private double   sumEnergySq = 0;
    private double[] sumCvcf;
    private double[] sumCvcfSq;
    private int      nSamples = 0;

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
        this.nSites = geo.nSites();
        this.config = new LatticeConfig(nSites, numComp);
        this.ncf = geo.getNcf();

        // 3. Pre-allocate scratch space for ΔE (O(1) allocation)
        this.maxEmbPerCol = Embeddings.maxEmbPerCfColumn(geo.cfEmbeddings());
        int scratchSize = ncf * maxEmbPerCol;
        this.scratch = new Embeddings.DeltaScratch(ncf, scratchSize);
        
        // 4. Initialize Correlation Functions state array
        // Layout: [uOrthNonPoint (ncf) | uPoint (numComp-1) | Empty (1)]
        this.correlationFunctions = new double[ncf + numComp];
        this.cvcfCorrelationFunctions = new double[ncf + numComp];

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
        updateCorrelationFunctions();
        
        System.out.println("\n--- Stage: Initial Randomisation ---");
        System.out.println("Measured Orthogonal CFs (uOrth):");
        for (int i = 0; i < correlationFunctions.length; i++) {
             System.out.println(String.format("  Col %2d: %.4f", i, correlationFunctions[i]));
        }
        
        System.out.println("\nRecovered CVCF Basis Vector (vFull):");
        double[] cvcf = getCvcfCorrelationFunctions();
        List<String> names = geo.getBasis().cfNames;
        List<String> defs = geo.getBasis().cfDefinitions;
        for (int i = 0; i < cvcf.length; i++) {
            String name = (i < names.size()) ? names.get(i) : "?";
            String def = (i < defs.size()) ? defs.get(i) : "";
            System.out.println(String.format("  Col %2d [%-7s]: %-40s = %.4f", i, name, def, cvcf[i]));
        }

        // 2. Initial state energy
        double currentEnergy = calculateTotalEnergy();
        System.out.println("\n--- Stage: Initial Energy ---");
        System.out.println(String.format("  Total Energy: %.6f", currentEnergy));
        System.out.println(String.format("  Energy per site: %.6f", currentEnergy / getNSites()));

        // 3. Equilibration phase
        System.out.println(String.format("\n--- Stage: Equilibration (%d sweeps) ---", nEquil));
        for (int s = 0; s < nEquil; s++) {
            currentEnergy += runSweep();
            if (s % 100 == 0 || s == nEquil - 1) {
                printState(s, currentEnergy / nSites);
            }
        }
        
        // Final refresh after equilibration to reset any drift
        updateCorrelationFunctions();
        currentEnergy = calculateTotalEnergy(); 
        System.out.println(String.format("  Energy after equilibration (refreshed): %.6f", currentEnergy));

        resetCounters();

        // 4. Averaging phase
        System.out.println(String.format("\n--- Stage: Averaging (%d sweeps) ---", nAvg));
        int sampleInterval = 20;
        for (int s = 0; s < nAvg; s++) {
            currentEnergy += runSweep();
            
            // Sub-sample to avoid expensive full O(N) updates every sweep
            if (s % sampleInterval == 0 || s == nAvg - 1) {
                updateCorrelationFunctions();
                sampleProperties(currentEnergy);
            }
            
            if (s % 100 == 0 || s == nAvg - 1) {
                printState(s, currentEnergy / nSites);
            }
        }

        LOG.info(String.format("AlloyMC run complete. Final energy: %.6f, Acceptance rate: %.2f%%", 
                currentEnergy, getAcceptanceRate() * 100));
    }

    /**
     * Prints the current simulation state including energy, acceptance rate, and CFs.
     */
    private void printState(int sweep, double energyPerSite) {
        // Ensure CFs are fresh for printing
        updateCorrelationFunctions();
        
        System.out.println(String.format("  Sweep %4d | E/site: %.6f | Acc: %.2f%%", 
                sweep, energyPerSite, getAcceptanceRate() * 100));
        
        double[] cvcf = getCvcfCorrelationFunctions();
        System.out.print("    CVCFs (first 5): ");
        for (int i = 0; i < Math.min(5, ncf); i++) {
            System.out.print(String.format("%.4f ", cvcf[i]));
        }
        System.out.println();
    }

    /**
     * Executes one sweep (N trial moves).
     * @return Total change in energy for the sweep.
     */
    private double runSweep() {
        double sweepDeltaE = 0;
        for (int i = 0; i < nSites; i++) {
            sweepDeltaE += attemptSwap();
        }
        return sweepDeltaE;
    }

    /**
     * Attempts a single canonical (swap) move between two random sites.
     * @return Change in energy if the move is accepted, 0.0 otherwise.
     */
    private double attemptSwap() {
        int i = rng.nextInt(nSites);
        int j = rng.nextInt(nSites);

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
        
        System.out.println("\n--- Energy Audit: Dot Product Calculation ---");
        System.out.println(String.format("  %-6s | %-12s | %-12s | %-12s", "Col", "ECI", "CVCF", "Product"));
        
        // 2. Compute Dot Product: Energy per site using cached CVCFs
        double energyPerSite = 0.0;
        for (int l = 0; l < ncf && l < eciCvcf.length; l++) {
            double term = eciCvcf[l] * cvcfCorrelationFunctions[l];
            energyPerSite += term;
            System.out.println(String.format("  Col %2d: %12.4f * %12.6f = %12.6f", 
                    l, eciCvcf[l], cvcfCorrelationFunctions[l], term));
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
        
        // Point CFs: Average of basis functions over the lattice.
        // We use the geometry's metadata (cfBasisIndices) to ensure the powers (k)
        // are mapped to the correct columns in the state vector.
        double[] x = config.composition();
        double[] basisSeq = org.ce.model.cluster.ClusterMath.buildBasis(numComp);
        int[][] cfBasisIndices = geo.getCfBasisIndices();
        for (int i = 0; i < numComp - 1; i++) {
            int col = ncf + i;
            int k = cfBasisIndices[col][0]; // Point CFs have single-index basis [k]
            double phiK = 0.0;
            for (int s = 0; s < numComp; s++) {
                phiK += x[s] * Math.pow(basisSeq[s], k);
            }
            correlationFunctions[col] = phiK;
        }
        
        // Empty cluster is always 1.0
        correlationFunctions[ncf + numComp - 1] = 1.0;

        // 3. Update CVCF cache
        this.cvcfCorrelationFunctions = geo.getCvcfCorrelationFunctions(correlationFunctions, composition);
    }
    
    public double[] getCvcfCorrelationFunctions() {
        return cvcfCorrelationFunctions;
    }

    /**
     * Calculates the change in energy for swapping occupations at sites i and j.
     * Uses a high-performance incremental update logic (O(1) complexity).
     */
    private double calculateDeltaE(int i, int j, int oldOccI, int oldOccJ) {
        if (eciCvcf == null || eciCvcf.length == 0) return 0.0;
        if (oldOccI == oldOccJ) return 0.0;

        int[] occ = config.getRawOcc();
        Embeddings.FlatEmbData flat = geo.getFlatEmbData();
        double[] flatBasisMatrix = geo.getFlatBasisMatrix();
        Embeddings.CsrSiteToCfIndex siteToCfIndex = geo.getSiteToCfIndex();
        double[][] Tinv = geo.getBasis().Tinv;

        int ac = 0;
        int startI = siteToCfIndex.offsets[i], endI = siteToCfIndex.offsets[i + 1];
        for (int idx = startI; idx < endI; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        int startJ = siteToCfIndex.offsets[j], endJ = siteToCfIndex.offsets[j + 1];
        for (int idx = startJ; idx < endJ; idx++) {
            int l  = siteToCfIndex.dataL[idx];
            int ei = siteToCfIndex.dataEI[idx];
            int key = l * maxEmbPerCol + ei;
            if (!scratch.seen[key]) {
                scratch.seen[key] = true;
                scratch.affectedL[ac]  = l;
                scratch.affectedEI[ac] = ei;
                ac++;
            }
        }
        scratch.affectedCount = ac;

        // 1. Calculate old local products
        for (int a = 0; a < ac; a++) {
            int l   = scratch.affectedL[a];
            int ei  = scratch.affectedEI[a];
            int flatEmbIdx = flat.cfOffsets[l] + ei;
            int sStart = flat.embSiteStart[flatEmbIdx];
            int sEnd   = flat.embSiteStart[flatEmbIdx + 1];
            double prod = 1.0;
            for (int k = sStart; k < sEnd; k++)
                prod *= flatBasisMatrix[occ[flat.siteData[k]] * numComp + flat.alphaData[k]];
            scratch.oldSumDelta[l] += prod;
        }

        // 2. Simulate Swap
        occ[i] = oldOccJ;
        occ[j] = oldOccI;

        // 3. Calculate new local products
        for (int a = 0; a < ac; a++) {
            int l   = scratch.affectedL[a];
            int ei  = scratch.affectedEI[a];
            int flatEmbIdx = flat.cfOffsets[l] + ei;
            int sStart = flat.embSiteStart[flatEmbIdx];
            int sEnd   = flat.embSiteStart[flatEmbIdx + 1];
            double prod = 1.0;
            for (int k = sStart; k < sEnd; k++)
                prod *= flatBasisMatrix[occ[flat.siteData[k]] * numComp + flat.alphaData[k]];
            scratch.newSumDelta[l] += prod;
        }

        // 4. Revert Swap (Metropolis will do it properly if accepted)
        occ[i] = oldOccI;
        occ[j] = oldOccJ;

        // 5. Compute Orthogonal Change (du)
        int cc = 0;
        for (int a = 0; a < ac; a++) {
            int l = scratch.affectedL[a];
            if (!scratch.seenCol[l]) {
                scratch.seenCol[l] = true;
                scratch.affectedCols[cc++] = l;
                double diff = scratch.newSumDelta[l] - scratch.oldSumDelta[l];
                if (diff != 0.0)
                    scratch.deltaUOrth[l] = diff / flat.cfEmbCount[l];
            }
        }
        scratch.affectedColCount = cc;

        // 6. Transform to CVCF basis and compute dE
        double dE = 0.0;
        if (Tinv != null) {
            int tCols = Tinv[0].length;
            for (int l = 0; l < ncf && l < Tinv.length; l++) {
                double sum = 0.0;
                for (int m = 0; m < ncf && m < tCols; m++)
                    sum += Tinv[l][m] * scratch.deltaUOrth[m];
                scratch.deltaVCvcf[l] = sum;
            }
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaVCvcf[l];
        } else {
            // Binary Fallback
            for (int l = 0; l < ncf && l < eciCvcf.length; l++)
                dE += eciCvcf[l] * scratch.deltaUOrth[l];
        }
        
        // Final energy per supercell
        dE *= nSites;

        // 7. Cleanup scratch 'seen' flags
        for (int a = 0; a < ac; a++)
            scratch.seen[scratch.affectedL[a] * maxEmbPerCol + scratch.affectedEI[a]] = false;
        for (int a = 0; a < cc; a++) {
            int l = scratch.affectedCols[a];
            scratch.seenCol[l]     = false;
            scratch.oldSumDelta[l] = 0.0;
            scratch.newSumDelta[l] = 0.0;
            scratch.deltaUOrth[l]  = 0.0;
            scratch.deltaVCvcf[l]  = 0.0;
        }
        scratch.affectedColCount = 0;

        return dE;
    }

    /**
     * Samples thermodynamic properties during the averaging phase.
     */
    private void sampleProperties(double energy) {
        if (sumCvcf == null) {
            sumCvcf = new double[cvcfCorrelationFunctions.length];
            sumCvcfSq = new double[cvcfCorrelationFunctions.length];
        }
        sumEnergy += energy;
        sumEnergySq += energy * energy;
        
        for (int i = 0; i < cvcfCorrelationFunctions.length; i++) {
            double v = cvcfCorrelationFunctions[i];
            sumCvcf[i] += v;
            sumCvcfSq[i] += v * v;
        }
        nSamples++;
    }

    public double getAverageEnergyPerSite() {
        return nSamples == 0 ? 0 : (sumEnergy / nSamples) / nSites;
    }

    public double getStdDevEnergyPerSite() {
        if (nSamples <= 1) return 0;
        double avg = sumEnergy / nSamples;
        double var = (sumEnergySq / nSamples) - (avg * avg);
        return Math.sqrt(Math.max(0, var)) / nSites;
    }

    public double[] getAverageCvcf() {
        if (nSamples == 0) return cvcfCorrelationFunctions.clone();
        double[] avg = new double[sumCvcf.length];
        for (int i = 0; i < sumCvcf.length; i++) {
            avg[i] = sumCvcf[i] / nSamples;
        }
        return avg;
    }

    public double[] getStdDevCvcf() {
        if (nSamples <= 1 || sumCvcfSq == null) return new double[cvcfCorrelationFunctions.length];
        double[] std = new double[sumCvcf.length];
        for (int i = 0; i < sumCvcf.length; i++) {
            double avg = sumCvcf[i] / nSamples;
            double var = (sumCvcfSq[i] / nSamples) - (avg * avg);
            std[i] = Math.sqrt(Math.max(0, var));
        }
        return std;
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
