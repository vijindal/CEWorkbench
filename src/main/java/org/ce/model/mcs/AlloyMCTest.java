package org.ce.model.mcs;

import org.ce.model.ModelSession;
import org.ce.model.storage.DataStore;
import org.ce.model.storage.Workspace;

import java.util.Random;
import java.util.Arrays;

/**
 * Quick test to verify AlloyMC constructor and state management.
 * Specifically targets the Nb-Ti BCC_A2 system.
 */
public class AlloyMCTest {
    public static void main(String[] args) {
        try {
            runTest();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void runTest() throws Exception {
        System.out.println("=== AlloyMC Constructor & State Test ===");
        
        // 1. Setup Data Infrastructure
        Workspace ws = new Workspace();
        // The HamiltonianStore constructor takes the Workspace instance
        DataStore.HamiltonianStore hStore = new DataStore.HamiltonianStore(ws);
        ModelSession.Builder builder = new ModelSession.Builder(hStore);
        
        // 2. Define System Identity (Nb-Ti-V BCC_A2 T)
        Workspace.SystemId id = new Workspace.SystemId("Nb-Ti-V", "BCC_A2", "T");
        
        // 3. Build ModelSession
        System.out.println("\n[1/3] Building ModelSession for Nb-Ti-V...");
        ModelSession session = builder.build(id, ModelSession.EngineConfig.MCS, msg -> {
            if (msg.startsWith("  [Session] ✓")) System.out.println(msg);
        });
        
        // 4. Instantiate AlloyMC
        System.out.println("\n[2/3] Instantiating AlloyMC (L=4)...");
        AlloyMC engine = new AlloyMC(session, 4, msg -> {
            if (msg.contains("built") || msg.contains("building")) System.out.println(msg);
        });
        
        // 5. Verify Fields & State
        System.out.println("\n[3/3] Verifying AlloyMC Internal State (Direct Getters):");
        System.out.println("--------------------------------------------------");
        System.out.println("Supercell Dim (L)      : " + engine.getL());
        System.out.println("Total Sites (N)        : " + engine.getNSites());
        System.out.println("Components (numComp)   : " + engine.getNumComp());
        System.out.println("Geometry CFs (ncf)     : " + engine.getNcf());
        System.out.println("Database Terms (CEC)   : " + (engine.getCecEntry() != null ? engine.getCecEntry().ncf : "NULL"));
        System.out.println("CF Vector length       : " + (engine.getCorrelationFunctions() != null ? engine.getCorrelationFunctions().length : "0"));
        System.out.println("Gas Constant (R)       : " + engine.getR_Gas() + " J/(mol*K)");
        
        System.out.println("\nTesting State Setters:");
        System.out.println("Setting T=1000K...");
        engine.setTemperature(1000, msg -> System.out.println("  " + msg));
        
        System.out.println("Setting x=[1.0, 0.0, 0.0] (Pure Nb)...");
        engine.setComposition(new double[]{1.0, 0.0, 0.0});
        engine.getConfig().randomise(engine.getComposition(), engine.getRng());
        System.out.println("Energy (Nb) : " + engine.calculateTotalEnergy() + " J");
        System.out.println("CFs (CVCF)  : " + Arrays.toString(engine.getCvcfCorrelationFunctions()));
        
        System.out.println("\nSetting x=[0.0, 1.0, 0.0] (Pure Ti)...");
        engine.setComposition(new double[]{0.0, 1.0, 0.0});
        engine.getConfig().randomise(engine.getComposition(), engine.getRng());
        System.out.println("Energy (Ti) : " + engine.calculateTotalEnergy() + " J");
        System.out.println("CFs (CVCF)  : " + Arrays.toString(engine.getCvcfCorrelationFunctions()));
        
        System.out.println("\nSetting x=[0.0, 0.0, 1.0] (Pure V)...");
        engine.setComposition(new double[]{0.0, 0.0, 1.0});
        engine.getConfig().randomise(engine.getComposition(), engine.getRng());
        System.out.println("Energy (V)  : " + engine.calculateTotalEnergy() + " J");
        System.out.println("CFs (CVCF)  : " + Arrays.toString(engine.getCvcfCorrelationFunctions()));
        
        System.out.println("\nFinal Verification:");
        System.out.println("Temperature       : " + engine.getTemperature() + " K");
        System.out.println("Composition       : " + Arrays.toString(engine.getComposition()));
        System.out.println("ECI Vector Size   : " + (engine.getEciCvcf() != null ? engine.getEciCvcf().length : 0));
        
        System.out.println("\n[4/4] Randomizing and Computing Energy...");
        engine.getConfig().randomise(engine.getComposition(), engine.getRng());
        double eInitial = engine.calculateTotalEnergy();
        System.out.println("Initial Energy    : " + eInitial + " J");
        System.out.println("CFs (Orthogonal)  : " + Arrays.toString(engine.getCorrelationFunctions()));
        System.out.println("CFs (CVCF Basis)  : " + Arrays.toString(engine.getCvcfCorrelationFunctions()));
        
        if (engine.getEciCvcf() != null && engine.getEciCvcf().length > 0) {
            System.out.println("First few ECIs    : " + 
                engine.getEciCvcf()[0] + ", " + 
                (engine.getEciCvcf().length > 1 ? engine.getEciCvcf()[1] : "N/A"));
        }
        
        System.out.println("\n=== Test Complete: PASSED ===");
    }
}
