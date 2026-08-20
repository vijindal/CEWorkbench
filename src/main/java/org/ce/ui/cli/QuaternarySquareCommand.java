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
import org.ce.calculation.workflow.QuaternarySquareScan;
import org.ce.calculation.workflow.SquarePlotRenderer;
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
 * JSON request/response entry point for a quaternary (X,Y)-square
 * composition scan at fixed temperature, used to render "square plots"
 * (Fig. 20, Jindal &amp; Lele 2025, CALPHAD 89, 102825).
 *
 * <p>Same stdin/stdout/stderr discipline as {@link TernaryGridCommand}: one
 * JSON request in, one JSON response out on the real stdout, physics/pipeline
 * diagnostics redirected to stderr. Exit code 0 on success, 1 on any error.</p>
 *
 * <p>Reuses {@link QuaternarySquareScan}, the same in-process engine the
 * GUI's square-plot panel calls directly.</p>
 *
 * <p><b>Always computes exactly two fixed slot orderings</b> for the given
 * elements A-B-C-D (in the order given in {@code system.elements}): A-B-C-D
 * and A-B-D-C (last two swapped). There is no {@code slotOrder} or
 * {@code variant} request field — together, these two square
 * parametrizations reach all six binary edges of the quaternary composition
 * tetrahedron (see {@link QuaternarySquareScan.Variant}'s doc for why a
 * single square cannot), so both are produced on every call.</p>
 *
 * <p>Request shape:</p>
 * <pre>
 * {"system": {"elements":"Nb-Ti-V-Zr","structure":"BCC_A2","model":"T","engine":"CVM"},
 *  "calculation": "GIBBS_ENERGY",
 *  "temperature": 1273,
 *  "n": 50,
 *  "render": true}
 * </pre>
 *
 * <p>{@code "calculation"} is either a standard {@code Property} name
 * ({@code GIBBS_ENERGY}, {@code ENTHALPY}, {@code ENTROPY}) or the literal
 * {@code "SRO"}, in which case a {@code "pair"} field is required — a
 * 2-element array naming the unlike species pair. Only 1st-neighbour pair
 * SRO is exposed (see {@link QuaternarySquareScan.Quantity}).</p>
 *
 * <p>{@code "render"} is optional (default {@code false}). When {@code true},
 * {@code scripts/square_section.py} is invoked on each computed grid, and
 * each entry of the response's {@code "results"} array gains an
 * {@code "image"} object: {@code {"format":"png","base64":"..."}}. Points
 * are always included regardless of {@code render}. Rendering requires a
 * working {@code python} on PATH with matplotlib installed; a rendering
 * failure is reported per-entry via {@code "renderError"} without failing
 * the whole request.</p>
 *
 * <p>Response shape: a top-level {@code "results"} array with exactly two
 * entries, tagged {@code "slotOrder":"A-B-C-D"} and
 * {@code "slotOrder":"A-B-D-C"} respectively (using the caller's actual
 * element symbols), each carrying its own {@code points}/{@code skipped}/
 * {@code image}.</p>
 */
public final class QuaternarySquareCommand {

    private static final ObjectMapper IN = new ObjectMapper()
            .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);
    private static final ObjectMapper OUT = new ObjectMapper();

    private QuaternarySquareCommand() {}

    public static int run(CEWorkbenchContext appCtx, InputStream stdin) {
        return run(appCtx, stdin, false);
    }

    /**
     * Runs the quaternary_square subcommand. Returns the process exit code.
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
                        + "Run 'quaternary_square --help' for usage, or see 'usage' below.");
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
     * quaternary_square --help} or an empty/missing stdin request) without
     * needing to find API.md first.
     */
    private static ObjectNode usage() {
        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        out.put("command", "quaternary_square");
        out.put("summary", "JSON request on stdin -> JSON response on stdout. "
                + "Computes a quaternary (X,Y)-square composition scan (Fig. 20 'square plot') "
                + "at fixed temperature for a 4-component system. Always returns TWO results, "
                + "for slot orders A-B-C-D and A-B-D-C (elements as given in system.elements), "
                + "since together they cover all 6 binary edges of the composition tetrahedron.");

        ObjectNode example = out.putObject("example");
        ObjectNode req = OUT.createObjectNode();
        ObjectNode sys = req.putObject("system");
        sys.put("elements", "Nb-Ti-V-Zr");
        sys.put("structure", "BCC_A2");
        sys.put("model", "T");
        sys.put("engine", "CVM");
        req.put("calculation", "GIBBS_ENERGY");
        req.put("temperature", 1273);
        req.put("n", 50);
        req.put("render", false);
        example.set("request", req);

        ArrayNode fields = out.putArray("fields");
        fields.add("system.elements: exactly 4 elements, e.g. 'Nb-Ti-V-Zr' (required)");
        fields.add("system.structure, system.model: as in 'api' (required)");
        fields.add("system.engine: 'CVM' or 'MCS', default 'CVM'");
        fields.add("calculation: GIBBS_ENERGY | ENTHALPY | ENTROPY | SRO, default GIBBS_ENERGY");
        fields.add("pair: required when calculation is 'SRO' - a 2-element array of element names "
                + "from system.elements, e.g. [\"Nb\",\"Ti\"] (1st-neighbour pair SRO only)");
        fields.add("temperature: fixed temperature in K (required)");
        fields.add("n: grid subdivisions per square axis, default 50 (point count = (n+1)^2)");
        fields.add("render: if true, also return a rendered PNG as base64 under each result's 'image' "
                + "(requires python + matplotlib on PATH); failures reported per-result via "
                + "'renderError' without failing the request");

        ArrayNode errorCodes = out.putArray("errorCodes");
        for (String code : new String[] {
                "EMPTY_REQUEST", "INVALID_REQUEST", "INTERNAL_ERROR", "SERIALIZATION_FAILED" }) {
            errorCodes.add(code);
        }

        out.put("relatedCommand", "ternary_grid --help  (3-component isothermal sections)");
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
        if (elements.size() != 4)
            throw new IllegalArgumentException(
                    "Quaternary square scan requires exactly 4 elements in system.elements, got '" + elementsRaw + "'.");

        String calcName = req.hasNonNull("calculation")
                ? req.get("calculation").asText().trim().toUpperCase()
                : Property.GIBBS_ENERGY.name();

        QuaternarySquareScan.Quantity quantity;
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
            quantity = new QuaternarySquareScan.PairSroQuantity(a, b);
        } else {
            quantity = new QuaternarySquareScan.PropertyQuantity(Property.valueOf(calcName));
        }

        if (!req.hasNonNull("temperature"))
            throw new IllegalArgumentException("Request must specify 'temperature'.");
        double temperature = req.get("temperature").asDouble();

        int n = req.hasNonNull("n") ? req.get("n").asInt() : 50;
        if (n < 1) throw new IllegalArgumentException("'n' must be >= 1.");

        boolean render = req.hasNonNull("render") && req.get("render").asBoolean();

        CalculationService service = appCtx.getCalculationService();
        ModelSpecifications specs = new ModelSpecifications(elementsRaw, structure, model, engine);
        ModelSession session = service.getOrBuildSession(specs, null);

        // Fixed pair: A-B-C-D as given, and A-B-D-C (last two swapped) — together
        // these reach all six binary edges of the composition tetrahedron.
        List<String> slotOrderAbcd = elements;
        List<String> slotOrderAbdc = List.of(elements.get(0), elements.get(1), elements.get(3), elements.get(2));

        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        out.put("structure", structure);
        out.put("model", model);
        out.put("engine", engine.name());
        out.put("temperature", temperature);
        out.put("calculation", quantity.label());

        ArrayNode results = out.putArray("results");
        results.add(runOne(service, session, slotOrderAbcd, temperature, quantity, n, structure, model, render));
        results.add(runOne(service, session, slotOrderAbdc, temperature, quantity, n, structure, model, render));

        return out;
    }

    private static ObjectNode runOne(CalculationService service, ModelSession session, List<String> slotOrder,
            double temperature, QuaternarySquareScan.Quantity quantity, int n,
            String structure, String model, boolean render) throws Exception {

        QuaternarySquareScan.Result result = QuaternarySquareScan.run(
                service, session, slotOrder, QuaternarySquareScan.Variant.STANDARD, temperature, quantity, n, null);

        ObjectNode out = OUT.createObjectNode();
        out.put("slotOrder", String.join("-", slotOrder));
        ArrayNode els = out.putArray("elements");
        slotOrder.forEach(els::add);
        out.put("skipped", result.skipped());

        ArrayNode points = out.putArray("points");
        for (QuaternarySquareScan.Point p : result.points()) {
            ObjectNode pn = OUT.createObjectNode();
            pn.put("x", p.x());
            pn.put("y", p.y());
            pn.put(slotOrder.get(0), p.fSlot0());
            pn.put(slotOrder.get(1), p.fSlot1());
            pn.put(slotOrder.get(2), p.fSlot2());
            pn.put(slotOrder.get(3), p.fSlot3());
            pn.put("value", p.value());
            pn.put("region", p.region().name());
            points.add(pn);
        }

        if (render) {
            try {
                File pngFile = SquarePlotRenderer.render(slotOrder, structure, model, temperature, result);
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
