package org.ce.ui.cli;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.CalculationDescriptor.ModelSpecifications;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.calculation.ConditionsScan;
import org.ce.calculation.EciValidator;
import org.ce.calculation.QuantityDescriptor;
import org.ce.calculation.Range;
import org.ce.calculation.workflow.CalculationService;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.ThermodynamicResult;
import org.ce.model.cvm.CvCfBasis;
import org.ce.model.cvm.SroCalculator;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace.SystemId;

import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * JSON request/response entry point for external (non-JVM) callers.
 *
 * <p>Reads one JSON request from stdin, writes one JSON response to stdout, and
 * sends all diagnostics to stderr so the payload is never polluted. Exit code is
 * 0 on success, 1 on any error (with a JSON error object still written to stdout
 * so callers can parse the reason).</p>
 *
 * <p>Two Hamiltonian modes:</p>
 * <ul>
 *   <li>{@code "hamiltonian"} present — caller supplies ECIs (must already be in
 *       the CVCF basis, J/mol). Validated strictly before use.</li>
 *   <li>{@code "hamiltonian"} omitted — the stored CEC database is used.</li>
 * </ul>
 */
public final class ApiCommand {

    /** Lenient mapper: unknown request keys (e.g. "$schema", "comment") are ignored. */
    private static final ObjectMapper IN = new ObjectMapper()
            .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false);

    private static final ObjectMapper OUT = new ObjectMapper();

    private ApiCommand() {}

    /** Runs the api subcommand. Returns the process exit code. */
    public static int run(CEWorkbenchContext appCtx, InputStream stdin) {
        // Physics and pipeline code writes progress directly to System.out (e.g.
        // "[Stage 3b] Testing..."). That would corrupt the JSON payload, so stdout is
        // redirected to stderr for the duration and the response is written to the
        // real stdout captured here.
        final PrintStream realOut = System.out;
        final PrintStream errOut = new PrintStream(new FileOutputStream(FileDescriptor.err), true);
        System.setOut(errOut);

        ObjectNode response;
        int exit = 0;
        try {
            JsonNode req = IN.readTree(new String(stdin.readAllBytes(), StandardCharsets.UTF_8));
            if (req == null || req.isNull() || req.isMissingNode())
                throw new ApiError("EMPTY_REQUEST", "No JSON request received on stdin.");

            response = req.has("describe")
                    ? handleDescribe(req.get("describe"))
                    : handleCalculate(appCtx, req);

        } catch (ApiError e) {
            response = e.toJson();
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

    // =========================================================================
    // describe — capability + expected-ECI-name discovery
    // =========================================================================

    private static ObjectNode handleDescribe(JsonNode d) throws ApiError {
        SystemSpec sys = SystemSpec.parse(d);

        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        out.set("system", sys.toJson());

        boolean supported = CvCfBasis.isSupported(sys.structure, sys.model, sys.elements.size());
        out.put("supported", supported);

        if (supported) {
            List<String> names = CvCfBasis.getNonPointCfNames(sys.structure, sys.model, sys.elements.size());
            out.put("ncf", names.size());
            ArrayNode arr = out.putArray("expectedEciNames");
            names.forEach(arr::add);
        } else {
            ArrayNode keys = out.putArray("supportedSystems");
            CvCfBasis.supportedKeys().forEach(keys::add);
        }

        ArrayNode calcs = out.putArray("calculations");
        for (Property p : implementedProperties(sys.engine)) calcs.add(p.name());

        ArrayNode notImpl = out.putArray("notImplemented");
        notImpl.add("PHASE_EQUILIBRIUM");

        return out;
    }

    /** Properties this build can actually compute for the given engine. */
    private static List<Property> implementedProperties(EngineConfig engine) {
        return engine.isCvm()
                ? List.of(Property.GIBBS_ENERGY, Property.ENTHALPY, Property.ENTROPY)
                : List.of(Property.ENTHALPY, Property.CORRELATION_FUNCTIONS);
    }

    // =========================================================================
    // calculate
    // =========================================================================

    private static ObjectNode handleCalculate(CEWorkbenchContext appCtx, JsonNode req) throws Exception {
        SystemSpec sys = SystemSpec.parse(req.get("system"));

        // ── Property ──────────────────────────────────────────────────────────
        String calcName = req.hasNonNull("calculation")
                ? req.get("calculation").asText().trim().toUpperCase()
                : Property.GIBBS_ENERGY.name();

        Property property;
        try {
            property = Property.valueOf(calcName);
        } catch (IllegalArgumentException e) {
            throw new ApiError("NOT_IMPLEMENTED",
                    "Unknown calculation '" + calcName + "'.")
                    .with("calculation", calcName);
        }
        if (!implementedProperties(sys.engine).contains(property)) {
            throw new ApiError("NOT_IMPLEMENTED",
                    "Calculation '" + calcName + "' is not implemented for engine " + sys.engine + ".")
                    .with("calculation", calcName);
        }

        // ── Conditions ────────────────────────────────────────────────────────
        ConditionsScan scan = parseConditions(req.get("conditions"));

        // ── Session (inline ECIs or stored database) ──────────────────────────
        CalculationService service = appCtx.getCalculationService();
        ModelSpecifications specs =
                new ModelSpecifications(sys.elementsRaw, sys.structure, sys.model, sys.engine);

        boolean inline = req.hasNonNull("hamiltonian");
        ModelSession session;
        if (inline) {
            CECEntry entry = parseHamiltonian(req.get("hamiltonian"), sys);
            validateEcis(entry, sys);
            // Bypass getOrBuildSession: its cache is keyed only on (systemId, engine)
            // and would collide with a stored-Hamiltonian session for the same system.
            session = new ModelSession.Builder(appCtx.getHamiltonianStore())
                    .build(new SystemId(sys.elementsRaw, sys.structure, sys.model), sys.engine, entry, null);
        } else {
            session = service.getOrBuildSession(specs, null);
        }

        // ── Run ───────────────────────────────────────────────────────────────
        List<ThermodynamicResult> points =
                service.calculateScan(session, scan, property, null, null);

        // ── Response ──────────────────────────────────────────────────────────
        ObjectNode out = OUT.createObjectNode();
        out.put("ok", true);
        out.set("system", sys.toJson());
        out.put("hamiltonianSource", inline ? "inline" : "store");
        if (!inline) out.put("hamiltonianId", session.resolvedHamiltonianId);
        out.put("calculation", property.name());

        ArrayNode arr = out.putArray("points");
        for (ThermodynamicResult r : points) arr.add(pointToJson(r, sys, session));
        return out;
    }

    /** Serializes one result; absent (NaN) quantities are omitted, never emitted as NaN. */
    private static ObjectNode pointToJson(ThermodynamicResult r, SystemSpec sys, ModelSession session) {
        ObjectNode p = OUT.createObjectNode();
        p.put("temperature", r.temperature);

        // Minimizer status. A false value means the numbers below came from a run
        // that hit its iteration limit — callers must not treat them as equilibrium.
        if (r.converged != null) {
            p.put("converged", r.converged);
            if (!r.converged) {
                p.put("iterations", r.iterations);
                if (!Double.isNaN(r.finalGradientNorm))
                    p.put("finalGradientNorm", r.finalGradientNorm);
            }
        }

        ObjectNode comp = p.putObject("composition");
        List<String> order = session.elements();
        for (int i = 0; i < order.size() && i < r.composition.length; i++) {
            comp.put(order.get(i), r.composition[i]);
        }

        // Drive quantity output from QuantityDescriptor so this stays in sync as
        // quantities are added, and so NaN never reaches the JSON (invalid JSON).
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            if (q.available(r)) p.put(jsonNameOf(q), q.extract(r));
        }

        double[] cfs = (r.optimizedCFs != null) ? r.optimizedCFs : r.avgCFs;
        if (cfs != null) {
            ObjectNode cfNode = p.putObject("correlationFunctions");
            List<String> names = CvCfBasis.isSupported(sys.structure, sys.model, sys.elements.size())
                    ? CvCfBasis.getNonPointCfNames(sys.structure, sys.model, sys.elements.size())
                    : List.of();
            // The CF vector trails point CFs (mole fractions, already reported under
            // "composition"); emit only the named non-point entries.
            int n = names.isEmpty() ? cfs.length : Math.min(cfs.length, names.size());
            for (int i = 0; i < n; i++) {
                cfNode.put(names.isEmpty() ? ("cf_" + i) : names.get(i), cfs[i]);
            }
        }

        // Cowley-Warren SRO (Jindal & Lele 2025, Eq. 40), keyed by neighbour shell.
        // alpha < 0 ordering, 0 random, > 0 clustering.
        if (r.sro != null && !r.sro.isEmpty()) {
            ObjectNode sroNode = p.putObject("sro");
            r.sro.forEach((shell, pairs) -> {
                ObjectNode shellNode = sroNode.putObject(shell);
                for (SroCalculator.PairSro s : pairs) {
                    if (s.i >= order.size() || s.j >= order.size()) continue;
                    ObjectNode e = shellNode.putObject(order.get(s.i) + "-" + order.get(s.j));
                    e.put("alpha", s.alpha);
                    e.put("probability", s.probability);
                    e.put("random", s.reference);
                }
            });
        }
        return p;
    }

    private static String jsonNameOf(QuantityDescriptor q) {
        return switch (q) {
            case GIBBS_ENERGY -> "gibbsEnergy";
            case ENTHALPY     -> "enthalpy";
            case ENTROPY      -> "entropy";
            case HEAT_CAPACITY-> "heatCapacity";
            case STD_ENTHALPY -> "stdEnthalpy";
        };
    }

    // =========================================================================
    // Parsing helpers
    // =========================================================================

    private static ConditionsScan parseConditions(JsonNode c) throws ApiError {
        if (c == null || c.isNull())
            throw new ApiError("MISSING_CONDITIONS", "Request has no 'conditions' object.");

        Range t = parseRange(c.get("temperature"), "temperature");
        if (t == null)
            throw new ApiError("MISSING_CONDITIONS", "conditions.temperature is required.");

        JsonNode compNode = c.get("composition");
        if (compNode == null || !compNode.isObject() || compNode.isEmpty())
            throw new ApiError("MISSING_CONDITIONS",
                    "conditions.composition must be a non-empty object of element -> fraction (or range).");

        Map<String, Range> comp = new LinkedHashMap<>();
        for (Iterator<String> it = compNode.fieldNames(); it.hasNext(); ) {
            String el = it.next();
            Range r = parseRange(compNode.get(el), "composition." + el);
            if (r != null) comp.put(el, r);
        }

        try {
            return new ConditionsScan(t, comp);
        } catch (IllegalArgumentException e) {
            throw new ApiError("INVALID_CONDITIONS", e.getMessage());
        }
    }

    /** Accepts a scalar (fixed) or {"start":..,"end":..,"step":..} object. */
    private static Range parseRange(JsonNode n, String where) throws ApiError {
        if (n == null || n.isNull()) return null;
        if (n.isNumber()) return Range.fixed(n.asDouble());
        if (n.isObject()) {
            if (!n.hasNonNull("start") || !n.hasNonNull("end"))
                throw new ApiError("INVALID_CONDITIONS",
                        where + " range requires 'start' and 'end' (and 'step' when they differ).");
            double start = n.get("start").asDouble();
            double end   = n.get("end").asDouble();
            double step  = n.hasNonNull("step") ? n.get("step").asDouble() : 0.0;
            if (Math.abs(start - end) > 1e-9 && Math.abs(step) < 1e-12)
                throw new ApiError("INVALID_CONDITIONS",
                        where + " range has start != end but step is 0 — this would not terminate.");
            return new Range(start, end, step);
        }
        throw new ApiError("INVALID_CONDITIONS", where + " must be a number or a {start,end,step} object.");
    }

    private static CECEntry parseHamiltonian(JsonNode h, SystemSpec sys) throws ApiError {
        // Basis and units must be declared explicitly — CEWorkbench applies no
        // transformation, so a caller sending orthogonal/eV values would otherwise
        // silently get wrong physics.
        String basis = h.hasNonNull("basis") ? h.get("basis").asText().trim() : "CVCF";
        if (!basis.equalsIgnoreCase("CVCF"))
            throw new ApiError("UNSUPPORTED_BASIS",
                    "hamiltonian.basis must be 'CVCF' (got '" + basis + "'). "
                    + "CEWorkbench applies no basis transformation; convert before sending.");

        String units = h.hasNonNull("units") ? h.get("units").asText().trim() : "J/mol";
        if (!units.equalsIgnoreCase("J/mol"))
            throw new ApiError("UNSUPPORTED_UNITS",
                    "hamiltonian.units must be 'J/mol' (got '" + units + "').");

        JsonNode terms = h.get("cecTerms");
        if (terms == null || !terms.isArray() || terms.isEmpty())
            throw new ApiError("INVALID_HAMILTONIAN", "hamiltonian.cecTerms must be a non-empty array.");

        List<CECEntry.CECTerm> list = new ArrayList<>();
        for (JsonNode t : terms) {
            if (!t.hasNonNull("name"))
                throw new ApiError("INVALID_HAMILTONIAN", "Every cecTerms entry requires a 'name'.");
            CECEntry.CECTerm term = new CECEntry.CECTerm();
            term.name = t.get("name").asText();
            term.a = t.hasNonNull("a") ? t.get("a").asDouble() : 0.0;
            term.b = t.hasNonNull("b") ? t.get("b").asDouble() : 0.0;
            term.validate();
            list.add(term);
        }

        CECEntry entry = new CECEntry();
        entry.elements = sys.elementsRaw;
        entry.structurePhase = sys.structure;
        entry.model = sys.model;
        entry.cecUnits = "J/mol";
        entry.cecTerms = list.toArray(new CECEntry.CECTerm[0]);
        entry.ncf = list.size();
        entry.notes = "Supplied inline via api";
        return entry;
    }

    /**
     * Strict ECI check. Rejects both unmatched names and incomplete coverage —
     * either would silently leave interactions at 0.0 in {@code CECEvaluator}.
     */
    private static void validateEcis(CECEntry entry, SystemSpec sys) throws Exception {
        int K = sys.elements.size();
        if (!CvCfBasis.isSupported(sys.structure, sys.model, K))
            throw new ApiError("UNSUPPORTED_SYSTEM",
                    "No CVCF basis registered for " + sys.structure + "/" + sys.model + "/" + K
                    + " components. " + CvCfBasis.supportedSummary());

        List<String> expected = CvCfBasis.getNonPointCfNames(sys.structure, sys.model, K);

        // Delegates to the same alias-resolution the evaluator uses, so a name that
        // passes validation here is guaranteed to map at evaluation time.
        EciValidator.Result result = EciValidator.validate(entry, expected);
        if (!result.isValid()) {
            ApiError err = new ApiError("ECI_VALIDATION_FAILED",
                    "Supplied ECIs do not map cleanly onto the "
                    + sys.structure + "/" + sys.model + " basis. Unmapped ECIs would "
                    + "silently default to 0.0, producing wrong energies.");
            err.arrays.put("unmatched", result.unmatched);
            err.arrays.put("unmapped", result.unmapped);
            err.arrays.put("expected", result.expected);
            throw err;
        }
    }

    // =========================================================================
    // Value types
    // =========================================================================

    /** Parsed and validated {@code "system"} block. */
    private static final class SystemSpec {
        final String elementsRaw;
        final List<String> elements;
        final String structure;
        final String model;
        final EngineConfig engine;

        private SystemSpec(String elementsRaw, List<String> elements,
                           String structure, String model, EngineConfig engine) {
            this.elementsRaw = elementsRaw;
            this.elements = elements;
            this.structure = structure;
            this.model = model;
            this.engine = engine;
        }

        static SystemSpec parse(JsonNode s) throws ApiError {
            if (s == null || s.isNull())
                throw new ApiError("MISSING_SYSTEM", "Request has no 'system' object.");

            String elements  = text(s, "elements");
            String structure = text(s, "structure");
            String model     = text(s, "model");

            // The _CVCF suffix belongs only to the derived storage ID; CvCfBasis.generate
            // rejects anything but the bare model name.
            if (model.toUpperCase().endsWith("_CVCF"))
                throw new ApiError("INVALID_SYSTEM",
                        "system.model must be the bare model name (e.g. 'T'), not '" + model
                        + "'. The _CVCF suffix is applied internally to the stored Hamiltonian ID.");

            String engineStr = s.hasNonNull("engine") ? s.get("engine").asText().trim() : "CVM";
            EngineConfig engine;
            try {
                engine = EngineConfig.valueOf(engineStr.toUpperCase());
            } catch (IllegalArgumentException e) {
                throw new ApiError("INVALID_SYSTEM", "system.engine must be 'CVM' or 'MCS', got '" + engineStr + "'.");
            }

            SystemId id = new SystemId(elements, structure, model);
            List<String> list;
            try {
                list = id.elementList();
            } catch (IllegalArgumentException e) {
                throw new ApiError("INVALID_SYSTEM", e.getMessage());
            }
            return new SystemSpec(elements, list, structure, model, engine);
        }

        private static String text(JsonNode s, String field) throws ApiError {
            if (!s.hasNonNull(field) || s.get(field).asText().isBlank())
                throw new ApiError("INVALID_SYSTEM", "system." + field + " is required.");
            return s.get(field).asText().trim();
        }

        ObjectNode toJson() {
            ObjectNode n = OUT.createObjectNode();
            ArrayNode els = n.putArray("elements");
            elements.forEach(els::add);
            n.put("structure", structure);
            n.put("model", model);
            n.put("engine", engine.name());
            return n;
        }
    }

    /** Error carrying a machine-readable code plus optional detail arrays. */
    private static final class ApiError extends Exception {
        final String code;
        final Map<String, String> extras = new LinkedHashMap<>();
        final Map<String, List<String>> arrays = new LinkedHashMap<>();

        ApiError(String code, String message) {
            super(message);
            this.code = code;
        }

        ApiError with(String k, String v) {
            extras.put(k, v);
            return this;
        }

        ObjectNode toJson() {
            ObjectNode n = OUT.createObjectNode();
            n.put("ok", false);
            n.put("error", code);
            n.put("message", getMessage());
            extras.forEach(n::put);
            arrays.forEach((k, vals) -> {
                ArrayNode a = n.putArray(k);
                vals.forEach(a::add);
            });
            return n;
        }
    }
}
