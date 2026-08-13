package org.ce.ui.cli;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.ModelSpecifications;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.workflow.CalculationService;
import org.ce.calculation.workflow.TernaryGridScan;
import org.ce.calculation.workflow.TernaryPlotRenderer;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;

import java.io.File;
import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Base64;
import java.util.List;

/**
 * JSON request/response entry point for a ternary composition-grid scan at
 * fixed temperature, used to render isothermal sections.
 *
 * <p>Same stdin/stdout/stderr discipline as {@link ApiCommand}: one JSON
 * request in, one JSON response out on the real stdout, all physics/pipeline
 * diagnostics redirected to stderr so the payload stays parseable. Exit code
 * 0 on success, 1 on any error.</p>
 *
 * <p>Reuses {@link TernaryGridScan}, the same in-process engine the GUI's
 * ternary panel calls directly — this command exists so the same scan can be
 * driven from the CLI (e.g. by an external renderer) without paying
 * per-point subprocess/JSON overhead the way repeated {@code api} calls
 * would.</p>
 *
 * <p>Request shape:</p>
 * <pre>
 * {"system": {"elements":"Nb-Ti-V","structure":"BCC_A2","model":"T","engine":"CVM"},
 *  "calculation": "GIBBS_ENERGY",
 *  "temperature": 1273,
 *  "n": 20,
 *  "render": true}
 * </pre>
 *
 * <p>{@code "calculation"} is either a standard {@code Property} name
 * ({@code GIBBS_ENERGY}, {@code ENTHALPY}, {@code ENTROPY}) or the literal
 * {@code "SRO"}, in which case a {@code "pair"} field is required — a
 * 2-element array naming the unlike species pair, e.g.
 * {@code {"calculation":"SRO","pair":["Nb","Ti"]}}. Only 1st-neighbour pair
 * SRO is exposed for now (see {@link TernaryGridScan.Quantity}).</p>
 *
 * <p>{@code "render"} is optional (default {@code false}). When {@code true},
 * the same {@code scripts/isothermal_section.py} (mpltern) renderer the GUI
 * uses is invoked on the computed grid, and the response gains an
 * {@code "image"} object: {@code {"format":"png","base64":"..."}}. Points are
 * always included regardless of {@code render}, so callers that only want
 * raw numbers pay no rendering cost. Rendering requires a working
 * {@code python} on PATH with matplotlib/mpltern installed; a rendering
 * failure is reported via {@code "renderError"} without failing the whole
 * request (points are still returned).</p>
 */
public final class TernaryGridCommand {

    private static final ObjectMapper IN = new ObjectMapper()
            .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);
    private static final ObjectMapper OUT = new ObjectMapper();

    private TernaryGridCommand() {}

    public static int run(CEWorkbenchContext appCtx, InputStream stdin) {
        return run(appCtx, stdin, false);
    }

    /**
     * Runs the ternary_grid subcommand. Returns the process exit code.
     *
     * @param helpRequested if {@code true} (caller passed {@code --help}/{@code -h}),
     *        stdin is not read and a self-documenting usage payload is printed
     *        instead. An empty/missing stdin request also gets this payload
     *        (under the {@code "usage"} key) alongside its error.
     */
    public static int run(CEWorkbenchContext appCtx, InputStream stdin, boolean helpRequested) {
        final PrintStream realOut = System.out;
        final PrintStream errOut = new PrintStream(new FileOutputStream(FileDescriptor.err), true);
        System.setOut(errOut);

        ObjectNode response;
        int exit = 0;
        try {
            if (helpRequested) {
                response = usage();
                System.setOut(realOut);
                realOut.println(OUT.writerWithDefaultPrettyPrinter().writeValueAsString(response));
                realOut.flush();
                return 0;
            }

            JsonNode req = IN.readTree(new String(stdin.readAllBytes(), StandardCharsets.UTF_8));
            if (req == null || req.isNull() || req.isMissingNode()) {
                ObjectNode errJson = OUT.createObjectNode();
                errJson.put("ok", false);
                errJson.put("error", "EMPTY_REQUEST");
                errJson.put("message", "No JSON request received on stdin. "
                        + "Run 'ternary_grid --help' for usage, or see 'usage' below.");
                errJson.set("usage", usage());
                System.setOut(realOut);
                realOut.println(OUT.writerWithDefaultPrettyPrinter().writeValueAsString(errJson));
                realOut.flush();
                return 1;
            }
            response = handle(appCtx, req);
        } catch (IllegalArgumentException e) {
            response = OUT.createObjectNode();
            response.put("ok", false);
            response.put("error", "INVALID_REQUEST");
            response.put("message", String.valueOf(e.getMessage()));
            exit = 1;
        } catch (Exception e) {
            response = OUT.createObjectNode();
            response.put("ok", false);
            response.put("error", "INTERNAL_ERROR");
            response.put("message", String.valueOf(e.getMessage()));
            exit = 1;
        } finally {
            System.setOut(realOut);
        }

        try {
            realOut.println(OUT.writerWithDefaultPrettyPrinter().writeValueAsString(response));
        } catch (Exception e) {
            realOut.println("{\"ok\":false,\"error\":\"SERIALIZATION_FAILED\"}");
            exit = 1;
        }
        realOut.flush();
        return exit;
    }

    /**
     * Self-documenting usage payload — request shape, fields, error codes —
     * so an external caller can discover this API in-band ({@code
     * ternary_grid --help} or an empty/missing stdin request) without needing
     * to find API.md first.
     */
    private static ObjectNode usage() {
        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        out.put("command", "ternary_grid");
        out.put("summary", "JSON request on stdin -> JSON response on stdout. "
                + "Computes a ternary composition-grid scan (isothermal section) at fixed temperature "
                + "for a 3-component system. Full docs: API.md ('Ternary isothermal sections') in the repo root.");

        ObjectNode example = out.putObject("example");
        ObjectNode req = OUT.createObjectNode();
        ObjectNode sys = req.putObject("system");
        sys.put("elements", "Nb-Ti-V");
        sys.put("structure", "BCC_A2");
        sys.put("model", "T");
        sys.put("engine", "CVM");
        req.put("calculation", "GIBBS_ENERGY");
        req.put("temperature", 1273);
        req.put("n", 20);
        req.put("render", false);
        example.set("request", req);

        ObjectNode sroExample = out.putObject("sroExample");
        ObjectNode sroReq = OUT.createObjectNode();
        ObjectNode sroSys = sroReq.putObject("system");
        sroSys.put("elements", "Nb-Ti-V");
        sroSys.put("structure", "BCC_A2");
        sroSys.put("model", "T");
        sroSys.put("engine", "CVM");
        sroReq.put("calculation", "SRO");
        ArrayNode sroPair = sroReq.putArray("pair");
        sroPair.add("Nb");
        sroPair.add("Ti");
        sroReq.put("temperature", 1273);
        sroReq.put("n", 20);
        sroExample.set("request", sroReq);
        sroExample.put("note", "1st-neighbour Cowley-Warren pair SRO only for now "
                + "(2NN and triangle/tetrahedron multi-site SRO not yet exposed).");

        ArrayNode fields = out.putArray("fields");
        fields.add("system.elements: exactly 3 elements, e.g. 'Nb-Ti-V' (required)");
        fields.add("system.structure, system.model: as in 'api' (required)");
        fields.add("system.engine: 'CVM' or 'MCS', default 'CVM'");
        fields.add("calculation: GIBBS_ENERGY | ENTHALPY | ENTROPY | SRO, default GIBBS_ENERGY");
        fields.add("pair: required when calculation is 'SRO' - a 2-element array of element names "
                + "from system.elements, e.g. [\"Nb\",\"Ti\"] (1st-neighbour pair SRO only)");
        fields.add("temperature: fixed temperature in K (required)");
        fields.add("n: grid subdivisions per triangle edge, default 20 "
                + "(point count = (n+1)(n+2)/2)");
        fields.add("render: if true, also return a rendered PNG as base64 under 'image' "
                + "(requires python + matplotlib/mpltern on PATH); failures reported via "
                + "'renderError' without failing the request");

        ArrayNode errorCodes = out.putArray("errorCodes");
        for (String code : new String[] {
                "EMPTY_REQUEST", "INVALID_REQUEST", "INTERNAL_ERROR", "SERIALIZATION_FAILED" }) {
            errorCodes.add(code);
        }

        out.put("relatedCommand", "api --help  (single-point / 1-D scan calculations)");
        return out;
    }

    private static ObjectNode handle(CEWorkbenchContext appCtx, JsonNode req) throws Exception {
        JsonNode sys = req.get("system");
        if (sys == null || sys.isNull())
            throw new IllegalArgumentException("Request has no 'system' object.");

        String elementsRaw = text(sys, "elements");
        String structure = text(sys, "structure");
        String model = text(sys, "model");
        EngineConfig engine = sys.hasNonNull("engine")
                ? EngineConfig.valueOf(sys.get("engine").asText().trim().toUpperCase())
                : EngineConfig.CVM;

        List<String> elements = java.util.Arrays.asList(elementsRaw.split("-"));
        if (elements.size() != 3)
            throw new IllegalArgumentException(
                    "Ternary grid scan requires exactly 3 elements in system.elements, got '" + elementsRaw + "'.");

        String calcName = req.hasNonNull("calculation")
                ? req.get("calculation").asText().trim().toUpperCase()
                : Property.GIBBS_ENERGY.name();

        TernaryGridScan.Quantity quantity;
        if ("SRO".equals(calcName)) {
            JsonNode pairNode = req.get("pair");
            if (pairNode == null || !pairNode.isArray() || pairNode.size() != 2)
                throw new IllegalArgumentException(
                        "calculation 'SRO' requires a 'pair' field: a 2-element array of element names, e.g. [\"Nb\",\"Ti\"].");
            String a = pairNode.get(0).asText().trim();
            String b = pairNode.get(1).asText().trim();
            if (!elements.contains(a) || !elements.contains(b))
                throw new IllegalArgumentException(
                        "'pair' elements must both be in system.elements (" + elementsRaw + "), got [" + a + "," + b + "].");
            quantity = new TernaryGridScan.PairSroQuantity(a, b);
        } else {
            quantity = new TernaryGridScan.PropertyQuantity(Property.valueOf(calcName));
        }

        if (!req.hasNonNull("temperature"))
            throw new IllegalArgumentException("Request must specify 'temperature'.");
        double temperature = req.get("temperature").asDouble();

        int n = req.hasNonNull("n") ? req.get("n").asInt() : 20;
        if (n < 1) throw new IllegalArgumentException("'n' must be >= 1.");

        CalculationService service = appCtx.getCalculationService();
        ModelSpecifications specs = new ModelSpecifications(elementsRaw, structure, model, engine);
        ModelSession session = service.getOrBuildSession(specs, null);

        TernaryGridScan.Result result = TernaryGridScan.run(service, session, elements, temperature, quantity, n, null);

        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        ArrayNode els = out.putArray("elements");
        elements.forEach(els::add);
        out.put("structure", structure);
        out.put("model", model);
        out.put("engine", engine.name());
        out.put("temperature", temperature);
        out.put("calculation", quantity.label());
        out.put("skipped", result.skipped());

        ArrayNode points = out.putArray("points");
        for (TernaryGridScan.Point p : result.points()) {
            ObjectNode pn = OUT.createObjectNode();
            pn.put(elements.get(0), p.fa());
            pn.put(elements.get(1), p.fb());
            pn.put(elements.get(2), p.fc());
            pn.put("value", p.value());
            pn.put("interpolated", p.interpolated());
            points.add(pn);
        }

        boolean render = req.hasNonNull("render") && req.get("render").asBoolean();
        if (render) {
            try {
                File pngFile = TernaryPlotRenderer.render(elements, structure, model, temperature, result);
                byte[] bytes = Files.readAllBytes(pngFile.toPath());
                ObjectNode image = out.putObject("image");
                image.put("format", "png");
                image.put("base64", Base64.getEncoder().encodeToString(bytes));
            } catch (Exception e) {
                out.put("renderError", String.valueOf(e.getMessage()));
            }
        }

        return out;
    }

    private static String text(JsonNode s, String field) {
        if (!s.hasNonNull(field) || s.get(field).asText().isBlank())
            throw new IllegalArgumentException("system." + field + " is required.");
        return s.get(field).asText().trim();
    }
}
