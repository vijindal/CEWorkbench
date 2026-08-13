package org.ce.calculation.workflow;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

import java.io.File;
import java.nio.file.Files;
import java.util.List;

/**
 * Shells out to {@code scripts/isothermal_section.py} (mpltern) to render a
 * PNG from a {@link TernaryGridScan.Result}.
 *
 * <p>No maintained Java library offers ternary contour plotting (see
 * CLAUDE.md), so rendering is delegated to Python; Java remains the sole
 * source of truth for the physics — this class only turns Java's computed
 * numbers into a JSON file and shells out for pixels. Shared by the GUI's
 * {@code TernaryPlotPanel} and the CLI's {@code TernaryGridCommand} so the
 * subprocess invocation exists in exactly one place.</p>
 */
public final class TernaryPlotRenderer {

    private TernaryPlotRenderer() {}

    /**
     * Writes {@code result} as the JSON shape {@code isothermal_section.py
     * --from-json} expects, invokes the renderer, and returns the produced
     * PNG file. Throws with the renderer's stderr on failure.
     */
    public static File render(List<String> elements, String structure, String model,
            double temperature, TernaryGridScan.Result result) throws Exception {
        File jsonFile = File.createTempFile("ternary_grid_", ".json");
        jsonFile.deleteOnExit();
        writeGridJson(jsonFile, elements, structure, model, temperature, result);

        File pngFile = File.createTempFile("ternary_plot_", ".png");
        pngFile.deleteOnExit();
        renderWithPython(jsonFile, pngFile);
        return pngFile;
    }

    private static void writeGridJson(File file, List<String> elements, String structure, String model,
            double temperature, TernaryGridScan.Result result) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        ObjectNode root = mapper.createObjectNode();
        root.put("ok", true);
        ArrayNode els = root.putArray("elements");
        elements.forEach(els::add);
        root.put("structure", structure);
        root.put("model", model);
        root.put("temperature", temperature);
        root.put("calculation", result.quantity().label());
        root.put("skipped", result.skipped());

        ArrayNode points = root.putArray("points");
        for (TernaryGridScan.Point p : result.points()) {
            ObjectNode pn = mapper.createObjectNode();
            pn.put(elements.get(0), p.fa());
            pn.put(elements.get(1), p.fb());
            pn.put(elements.get(2), p.fc());
            pn.put("value", p.value());
            pn.put("interpolated", p.interpolated());
            points.add(pn);
        }
        Files.writeString(file.toPath(), mapper.writerWithDefaultPrettyPrinter().writeValueAsString(root));
    }

    private static void renderWithPython(File jsonFile, File pngFile) throws Exception {
        File repoRoot = findRepoRoot();
        File script = new File(repoRoot, "scripts/isothermal_section.py");
        if (!script.exists())
            throw new IllegalStateException("Renderer script not found: " + script);

        ProcessBuilder pb = new ProcessBuilder("python", script.getAbsolutePath(),
                "--from-json", jsonFile.getAbsolutePath(),
                "--out", pngFile.getAbsolutePath());
        pb.redirectErrorStream(false);
        Process proc = pb.start();
        String stderr = new String(proc.getErrorStream().readAllBytes());
        int exit = proc.waitFor();
        if (exit != 0 || !pngFile.exists()) {
            throw new IllegalStateException("Renderer failed (exit " + exit + "): "
                    + (stderr.isBlank() ? "no output" : stderr.strip()));
        }
    }

    private static File findRepoRoot() {
        File dir = new File(System.getProperty("user.dir"));
        while (dir != null) {
            if (new File(dir, "scripts/isothermal_section.py").exists()) return dir;
            dir = dir.getParentFile();
        }
        return new File(System.getProperty("user.dir"));
    }
}
