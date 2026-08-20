package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.CalculationService;
import org.ce.calculation.workflow.QuaternarySquareScan;
import org.ce.calculation.workflow.QuaternarySquareScan.PropertyQuantity;
import org.ce.calculation.workflow.QuaternarySquareScan.Result;
import org.ce.calculation.workflow.QuaternarySquareScan.Variant;
import org.ce.calculation.workflow.SquarePlotRenderer;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * End-to-end smoke test: computes a coarse quaternary square scan and
 * renders it through the real Python pipeline, saving the PNG next to
 * scripts/ for visual inspection.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.SquarePlotRenderSmokeTest
 * </pre>
 */
public final class SquarePlotRenderSmokeTest {

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        Workspace workspace = new Workspace();
        CEWorkbenchContext context = new CEWorkbenchContext(workspace);
        ModelSession session = new ModelSession.Builder(context.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
        CalculationService service = context.getCalculationService();

        List<String> slotOrder = List.of("Nb", "Ti", "V", "Zr");
        double temperature = 1273.0;
        int n = 20;

        for (Variant variant : Variant.values()) {
            System.out.println("Computing " + variant + " ...");
            Result result = QuaternarySquareScan.run(service, session, slotOrder, variant,
                    temperature, new PropertyQuantity(Property.GIBBS_ENERGY), n, null);
            System.out.printf("  points=%d skipped=%d%n", result.points().size(), result.skipped());

            File png = SquarePlotRenderer.render(slotOrder, "BCC_A2", "T", temperature, result);
            Path dest = Path.of("scripts", "Nb-Ti-V-Zr_" + variant + "_1273K_GIBBS_ENERGY_smoketest.png");
            Files.copy(png.toPath(), dest, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            System.out.println("  Saved: " + dest.toAbsolutePath());
        }

        int nEntropy = 50;
        System.out.println("Computing STANDARD entropy at n=" + nEntropy + " ...");
        Result entropyResult = QuaternarySquareScan.run(service, session, slotOrder, Variant.STANDARD,
                temperature, new PropertyQuantity(Property.ENTROPY), nEntropy, null);
        System.out.printf("  points=%d skipped=%d%n", entropyResult.points().size(), entropyResult.skipped());
        File entropyPng = SquarePlotRenderer.render(slotOrder, "BCC_A2", "T", temperature, entropyResult);
        Path entropyDest = Path.of("scripts", "Nb-Ti-V-Zr_STANDARD_1273K_ENTROPY_n50_smoketest.png");
        Files.copy(entropyPng.toPath(), entropyDest, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        System.out.println("  Saved: " + entropyDest.toAbsolutePath());

        System.out.println("\nRESULT: PASS");
    }

    private SquarePlotRenderSmokeTest() {}
}
