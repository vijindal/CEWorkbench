package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.ModelSession;
import org.ce.model.ThermodynamicResult;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.calculation.CalculationDescriptor.*;
import org.ce.calculation.Conditions;
import org.ce.calculation.workflow.CalculationService;

import java.util.Map;

public class ThermodynamicComparisonTest {
    public static void main(String[] args) {
        try {
            CEWorkbenchContext appCtx = new CEWorkbenchContext();
            CalculationService service = appCtx.getCalculationService();

            String elements = "Nb-Ti-V";
            String structure = "BCC_A2";
            String model = "T";
            double temp = 1000.0;

            System.out.println("=== Thermodynamic Comparison: Nb-Ti-V @ 1000K (Equiatomic) ===");
            System.out.println("System: " + elements + " / " + structure + " / " + model);

            // 1. Run CVM — via the named calculate() entry point (Nb derived).
            System.out.println("\n[1] Running CVM...");
            ModelSpecifications cvmSpecs = new ModelSpecifications(elements, structure, model, EngineConfig.CVM);
            ModelSession cvmSession = service.getOrBuildSession(cvmSpecs, System.out::println);
            Conditions conditions = new Conditions(temp, Map.of("Ti", 1.0/3.0, "V", 1.0/3.0));

            ThermodynamicResult cvmResult =
                    service.calculate(cvmSession, conditions, Property.GIBBS_ENERGY, System.out::println, null);

            // 2. Run AlloyMC (Direct Execution)
            int L = 8;
            int nEquil = 100;
            int nAvg = 500;
            System.out.println("\n[2] Running AlloyMC (Direct execution, L=" + L + ", " + nEquil + "+" + nAvg + " sweeps)...");
            org.ce.model.ModelSession mcsSession = new org.ce.model.ModelSession.Builder(appCtx.getHamiltonianStore())
                .build(new org.ce.model.storage.Workspace.SystemId(elements, structure, model), EngineConfig.MCS, null);
            
            org.ce.model.mcs.AlloyMC alloyMC = new org.ce.model.mcs.AlloyMC(mcsSession, L, null);
            alloyMC.setTemperature(temp);
            alloyMC.setComposition(new double[]{1.0/3.0, 1.0/3.0, 1.0/3.0});
            alloyMC.run(nEquil, nAvg); 
            
            ThermodynamicResult mcsResult = new ThermodynamicResult(
                temp,
                new double[]{1.0/3.0, 1.0/3.0, 1.0/3.0},
                Double.NaN,                         // gibbsEnergy
                alloyMC.getAverageEnergyPerSite(),  // enthalpy
                Double.NaN,                         // entropy
                alloyMC.getStdDevEnergyPerSite(),   // stdEnthalpy
                Double.NaN,                         // heatCapacity
                null,                               // optimizedCFs
                alloyMC.getAverageCvcf(),           // avgCFs
                alloyMC.getStdDevCvcf()             // stdCFs
            );

            // 3. Compare
            printComparison(cvmResult, mcsResult);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void printComparison(ThermodynamicResult cvm, ThermodynamicResult mcs) {
        System.out.println("\n=== Comparison of Results ===");
        System.out.format("Quantity        | CVM                  | MCS (avg \u00B1 \u03C3)         \n");
        System.out.println("-".repeat(60));
        System.out.format("G (J/mol)       | %20.4f | %20s\n", cvm.gibbsEnergy, formatVal(mcs.gibbsEnergy));
        System.out.format("H (J/mol)       | %20.4f | %10.4f \u00B1 %7.4f\n", cvm.enthalpy, mcs.enthalpy, mcs.stdEnthalpy);
        System.out.format("S (J/mol.K)     | %20.4f | %20s\n", cvm.entropy, formatVal(mcs.entropy));

        System.out.println("\n=== Correlation Function Comparison (CVCF Basis) ===");
        double[] cvmCFs = cvm.optimizedCFs;
        double[] mcsCFs = mcs.avgCFs;
        double[] stdCFs = mcs.stdCFs;

        if (cvmCFs == null || mcsCFs == null) {
            System.out.println("Error: CFs not available for one or both engines.");
            return;
        }

        System.out.println("CVM CF Array Length: " + cvmCFs.length);
        System.out.println("MCS CF Array Length: " + mcsCFs.length);

        System.out.format("%-4s | %-12s | %-20s | %-12s\n", "Idx", "CVM", "MCS (avg \u00B1 \u03C3)", "Diff");
        System.out.println("-".repeat(60));
        int n = Math.max(cvmCFs.length, mcsCFs.length);
        for (int i = 0; i < n; i++) {
            double v1 = (i < cvmCFs.length) ? cvmCFs[i] : Double.NaN;
            double v2 = (i < mcsCFs.length) ? mcsCFs[i] : Double.NaN;
            double s2 = (stdCFs != null && i < stdCFs.length) ? stdCFs[i] : 0.0;
            
            String cvmStr = (i < cvmCFs.length) ? String.format("%12.6f", v1) : "      N/A   ";
            String mcsStr = (i < mcsCFs.length) ? String.format("%8.6f \u00B1 %7.6f", v2, s2) : "      N/A   ";
            String diffStr = (i < cvmCFs.length && i < mcsCFs.length) 
                             ? String.format("%12.6f", v2 - v1) 
                             : "      N/A   ";
            
            System.out.println(String.format("[%02d] | %s | %s | %s", i, cvmStr, mcsStr, diffStr));
        }
    }

    private static String formatVal(double val) {
        return Double.isNaN(val) ? "N/A" : String.format("%.4f", val);
    }
}
