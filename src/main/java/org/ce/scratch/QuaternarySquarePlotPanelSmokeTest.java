package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.model.storage.Workspace;
import org.ce.ui.gui.ActivityBar;
import org.ce.ui.gui.ExplorerPanel;
import org.ce.ui.gui.MainWindow;
import org.ce.ui.gui.OutputPanel;
import org.ce.ui.gui.QuaternarySquarePlotPanel;
import org.ce.ui.gui.WorkbenchContext;

import javax.swing.*;
import java.awt.image.BufferedImage;

/**
 * Constructs the real Swing wiring for the quaternary square-plot panel
 * (without showing a window) to catch construction-time and array-index
 * bugs that {@code SquarePlotRenderSmokeTest} can't see, since that test
 * calls {@link org.ce.calculation.workflow.QuaternarySquareScan} and
 * {@link org.ce.calculation.workflow.SquarePlotRenderer} directly and never
 * touches Swing.
 *
 * <p>Covers: {@link QuaternarySquarePlotPanel} constructs without throwing;
 * {@link ExplorerPanel#addCard} accepts slot index 4;
 * {@code ActivityBar}'s 5-item {@code ITEMS} array and the 5-callback
 * {@code navCallbacks} array used in {@link MainWindow} stay in sync;
 * {@link OutputPanel#showSquarePlot}/{@code showSquareError} route to the
 * correct card without throwing. Does not cover layout/rendering/Nimbus
 * theming, which require an actual visible window.</p>
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.QuaternarySquarePlotPanelSmokeTest
 * </pre>
 */
public final class QuaternarySquarePlotPanelSmokeTest {

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(72));
        System.out.println("  QuaternarySquarePlotPanel Swing-wiring smoke test");
        System.out.println("=".repeat(72));

        SwingUtilities.invokeAndWait(QuaternarySquarePlotPanelSmokeTest::runOnEdt);

        System.out.println("\n" + "=".repeat(72));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(72));
        if (failures > 0) {
            throw new AssertionError(failures + " Swing-wiring checks failed");
        }
    }

    private static void runOnEdt() {
        try {
            CEWorkbenchContext appCtx = new CEWorkbenchContext(new Workspace());
            WorkbenchContext context = new WorkbenchContext();

            OutputPanel outputPanel = new OutputPanel(context, new org.ce.calculation.QuantityDescriptor.SelectionModel());
            check("OutputPanel constructs", outputPanel != null);

            QuaternarySquarePlotPanel panel = new QuaternarySquarePlotPanel(
                    appCtx, context, s -> {}, s -> {},
                    outputPanel::showSquarePlot, outputPanel::showSquareError);
            check("QuaternarySquarePlotPanel constructs without throwing", panel != null);

            // showSquarePlot / showSquareError must not throw when routed through
            // OutputPanel's CardLayout (catches a stale/missing CARD_SQUARE key).
            // showSquarePlot now takes two images (the fixed A-B-C-D / A-B-D-C pair).
            BufferedImage dummyA = new BufferedImage(4, 4, BufferedImage.TYPE_INT_RGB);
            BufferedImage dummyB = new BufferedImage(4, 4, BufferedImage.TYPE_INT_RGB);
            boolean threw = false;
            try {
                outputPanel.showSquarePlot(dummyA, dummyB);
                outputPanel.showSquareError("test error");
            } catch (Exception e) {
                threw = true;
                System.out.println("    [!] showSquarePlot/showSquareError threw: " + e);
            }
            check("showSquarePlot/showSquareError do not throw", !threw);

            // ExplorerPanel must accept the panel at slot 4 (the new "Square" card)
            // without throwing ArrayIndexOutOfBounds -- this is exactly the bug class
            // a mismatched TITLES/CARDS array length would produce.
            ExplorerPanel explorer = new ExplorerPanel();
            boolean explorerThrew = false;
            try {
                explorer.addCard(panel, 4);
                explorer.showCard(4);
            } catch (Exception e) {
                explorerThrew = true;
                System.out.println("    [!] ExplorerPanel slot 4 threw: " + e);
            }
            check("ExplorerPanel accepts card at slot 4", !explorerThrew);

            // ActivityBar must accept a 5-callback array matching its 5-item ITEMS
            // table -- this is exactly the bug class a length mismatch would produce
            // (ArrayIndexOutOfBoundsException inside the constructor's for-loop).
            Runnable[] callbacks = { () -> {}, () -> {}, () -> {}, () -> {}, () -> {} };
            boolean activityBarThrew = false;
            try {
                ActivityBar bar = new ActivityBar(callbacks);
                bar.setActive(4);
            } catch (Exception e) {
                activityBarThrew = true;
                System.out.println("    [!] ActivityBar 5-item construction threw: " + e);
            }
            check("ActivityBar accepts 5 callbacks matching 5-item ITEMS table", !activityBarThrew);

        } catch (Exception e) {
            failures++;
            System.out.println("    [!] FAIL  unexpected exception: " + e);
            e.printStackTrace();
        }
    }

    private static void check(String what, boolean ok) {
        System.out.println((ok ? "    [ok] " : "    [!] FAIL  ") + what);
        if (!ok) failures++;
    }

    private QuaternarySquarePlotPanelSmokeTest() {}
}
