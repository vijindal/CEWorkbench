package org.ce.model.equilibrium;

import org.ce.model.storage.Workspace;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Reader for the SGTE Unary database ({@code inputs/unary.dat}), supplying the
 * pure-element reference energies {@code G0(element, phase, T)} and their
 * temperature derivatives.
 *
 * <p>This is the database-backed counterpart to {@link LatticeStability},
 * which hardcodes the same polynomials for eight elements. Both are
 * transcriptions of the same SGTE source, so they can be cross-checked against
 * each other -- see {@code org.ce.scratch.LatticeStabilityVerification}. The
 * database version is preferred going forward because it covers every element
 * and phase in the file rather than a hand-picked subset, and because
 * transcription errors become impossible: a hardcoded coefficient can silently
 * lose a digit, and nothing in the program can detect it.</p>
 *
 * <h2>File format</h2>
 *
 * <p>Each element's block is delimited by a {@code "$ <El>"} header and
 * {@code "$ ----"} rules. Within it, one {@code PARAMETER G(<PHASE>,<EL>:VA;0)}
 * entry carries the piecewise polynomial, as a sequence of
 * {@code <expression>; <upperT> Y} segments terminated by {@code N !}:</p>
 *
 * <pre>
 *   PARAMETER G(BCC_A2,TI:VA;0)   298.15
 *    -1272.064+134.71418*T-25.5768*T*LN(T)-0.663845E-3*T**2-0.278803E-6*T**3
 *   +7208*T**(-1);  1155.00 Y
 *    6667.385+105.366379*T-22.3771*T*LN(T)+...;  1941.00 Y
 *    ...;  4000.00 N !
 * </pre>
 *
 * <p>Each expression is reduced to a coefficient vector over the fixed SGTE
 * term set, so evaluation and differentiation are both closed-form:</p>
 *
 * <pre>
 *   G0(T) = c0 + c1 T + c2 T ln T + c3 T^2 + c4 T^3 + c5 T^4
 *         + c6 T^7 + c7 T^-1 + c8 T^-2 + c9 T^-3 + c10 T^-9
 * </pre>
 *
 * <p>Adapted from a standalone {@code database.sgte} implementation. The
 * parsing logic is preserved; the surrounding behaviour is not. That original
 * called {@code System.exit(0)} on an unknown element or phase, resolved the
 * data file by walking up from the working directory, logged through a global
 * printer, and threw {@code IOException} from every accessor. Here lookups
 * throw {@link IllegalArgumentException}, the path comes from {@link Workspace},
 * and the parse cache is thread-safe -- the GUI evaluates on a SwingWorker
 * while scans may run concurrently.</p>
 */
public final class SgteDatabase {

    /** Number of SGTE polynomial terms; see the class doc for the ordering. */
    private static final int NUM_TERMS = 11;

    private static final String DATA_FILE = "unary.dat";

    /** Parsed piecewise polynomial for one (element, phase) pair. */
    private record Entry(String element, String phase, double[] breakpoints, String[] expressions, int numRanges) {

        /**
         * Coefficients valid at {@code t}. Ranges are half-open
         * {@code (lower, upper]} as the database defines them, except that the
         * first range includes its lower bound so the tabulated start
         * temperature (usually 298.15 K) is usable.
         */
        double[] coefficientsAt(double t) {
            for (int j = 0; j < numRanges - 1; j++) {
                boolean lowerOk = (j == 0) ? t >= breakpoints[j] : t > breakpoints[j];
                if (lowerOk && t <= breakpoints[j + 1]) {
                    return parseExpression(expressions[j]);
                }
            }
            throw new IllegalArgumentException(String.format(
                    "SGTE: %s %s has no data at T = %.2f K (valid range %.2f to %.2f K)",
                    element, phase, t, breakpoints[0], breakpoints[numRanges - 1]));
        }
    }

    /** Cache keyed on "EL-PHASE"; the file parse is far costlier than the lookup. */
    private static final Map<String, Entry> CACHE = new ConcurrentHashMap<>();

    private SgteDatabase() {
    }

    // =========================================================================
    // Public API -- mirrors LatticeStability so callers need not change
    // =========================================================================

    /**
     * Pure-element reference energy {@code G0(element, phase, T)} in J/mol.
     *
     * @param element element symbol, e.g. {@code "Ti"} (case-insensitive)
     * @param phase   phase name, e.g. {@code "BCC_A2"}; an ordered phase is
     *                mapped to its disordered parent by {@link #sgtePhaseName}
     * @param t       temperature in K
     */
    public static double g0(String element, String phase, double t) {
        return evaluate(coefficients(element, phase, t), t);
    }

    /** Analytic {@code dG0/dT} in J/mol/K. */
    public static double dG0Dt(String element, String phase, double t) {
        return evaluateFirstDerivative(coefficients(element, phase, t), t);
    }

    /** Analytic {@code d2G0/dT2} in J/mol/K^2. */
    public static double d2G0Dt2(String element, String phase, double t) {
        return evaluateSecondDerivative(coefficients(element, phase, t), t);
    }

    /**
     * Composition-weighted reference energy
     * {@code G0m = sum_i x_i * g0(element_i, phase, t)} -- the mechanical
     * mixture of pure elements in {@code G = G0m + Gm}.
     */
    public static double g0m(List<String> elements, String phase, double[] composition, double t) {
        if (elements.size() != composition.length) {
            throw new IllegalArgumentException(
                    "elements.size()=" + elements.size() + " != composition.length=" + composition.length);
        }
        double sum = 0.0;
        for (int i = 0; i < elements.size(); i++) {
            sum += composition[i] * g0(elements.get(i), phase, t);
        }
        return sum;
    }

    /** Composition-weighted {@code d(G0m)/dT}. */
    public static double dG0mDt(List<String> elements, String phase, double[] composition, double t) {
        if (elements.size() != composition.length) {
            throw new IllegalArgumentException(
                    "elements.size()=" + elements.size() + " != composition.length=" + composition.length);
        }
        double sum = 0.0;
        for (int i = 0; i < elements.size(); i++) {
            sum += composition[i] * dG0Dt(elements.get(i), phase, t);
        }
        return sum;
    }

    /**
     * Maps a structure name to the phase whose pure-element reference applies.
     *
     * <p>An ordered phase takes the reference of its disordered parent -- a
     * pure element has no ordering, so {@code BCC_B2} and {@code BCC_A2} share
     * one {@code G0}. Names already in SGTE form pass through unchanged.</p>
     */
    public static String sgtePhaseName(String phase) {
        return switch (phase.toUpperCase()) {
            case "BCC_A2", "A2", "BCC_B2", "B2", "B32" -> "BCC_A2";
            case "FCC_A1", "A1", "L10", "L12" -> "FCC_A1";
            case "HCP_A3", "A3", "B19", "D019" -> "HCP_A3";
            case "DIAMOND_A4", "A4" -> "DIAMOND_A4";
            case "LIQUID", "L" -> "LIQUID";
            default -> phase.toUpperCase();
        };
    }

    // =========================================================================
    // Lookup and caching
    // =========================================================================

    private static double[] coefficients(String element, String phase, double t) {
        String el = normaliseElement(element);
        String ph = sgtePhaseName(phase);
        Entry entry = CACHE.computeIfAbsent(el + "-" + ph, key -> read(el, ph));
        return entry.coefficientsAt(t);
    }

    /** SGTE writes element symbols in upper case ({@code TI}, not {@code Ti}). */
    private static String normaliseElement(String element) {
        return element.toUpperCase();
    }

    private static Path dataFile() {
        return new Workspace().inputsDir().resolve(DATA_FILE);
    }

    /**
     * Scans the database for one (element, phase) block.
     *
     * <p>Element blocks are located by their {@code "$ <EL>"} header and bounded
     * by {@code "$ ----"} rules, so a {@code PARAMETER} line for a different
     * element cannot be picked up by accident.</p>
     */
    private static Entry read(String element, String phase) {
        Path path = dataFile();
        if (!Files.exists(path)) {
            throw new IllegalStateException(
                    "SGTE database not found at " + path + ". Expected " + DATA_FILE
                            + " under the workspace inputs directory.");
        }

        try (BufferedReader br = Files.newBufferedReader(path, StandardCharsets.ISO_8859_1)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.trim().equalsIgnoreCase("$ " + element)) {
                    continue;
                }
                Entry entry = readElementBlock(br, element, phase);
                if (entry != null) {
                    return entry;
                }
            }
        } catch (IOException e) {
            throw new IllegalStateException("Failed reading SGTE database at " + path, e);
        }

        throw new IllegalArgumentException(
                "SGTE: no data for " + element + " " + phase + " in " + path.getFileName()
                        + " (check the element symbol and phase name)");
    }

    /** Reads one element's block, returning the requested phase or null. */
    private static Entry readElementBlock(BufferedReader br, String element, String phase)
            throws IOException {

        // Solid phases carry the vacancy sublattice, "G(BCC_A2,TI:VA;0)";
        // LIQUID does not, "G(LIQUID,TI;0)". Accept either terminator so the
        // element name cannot match a prefix of a longer symbol.
        String stem = "PARAMETER G(" + phase + "," + element;
        boolean started = false;
        String line;

        while ((line = br.readLine()) != null) {
            if (line.startsWith("$ ------")) {
                if (!started) {
                    started = true;
                    continue;
                }
                return null; // end of this element's block
            }
            String upper = line.toUpperCase();
            if (!upper.startsWith(stem + ":") && !upper.startsWith(stem + ";")) {
                continue;
            }

            // The parameter's segments may wrap across lines; accumulate until
            // the terminating "N !".
            StringBuilder body = new StringBuilder(line.substring(line.indexOf(')') + 1));
            while (!body.toString().trim().endsWith("!")) {
                String next = br.readLine();
                if (next == null) {
                    break;
                }
                body.append(' ').append(next);
            }
            return parseParameter(element, phase, body.toString());
        }
        return null;
    }

    /**
     * Parses a parameter body into breakpoints and per-range expressions.
     *
     * <p>Shape: {@code <T0> <expr>; <T1> Y <expr>; <T2> N !} -- a start
     * temperature, then repeating (expression, upper bound) pairs. Expressions
     * end at {@code ';'}; the token after each bound is a continuation flag
     * ({@code Y}/{@code N}) that carries no coefficient information.</p>
     */
    private static Entry parseParameter(String element, String phase, String body) {
        // Segments are separated by ';'. The text before the first ';' is
        // "<T_start> <expr>"; each later segment is "<T_upper> <Y|N> <expr>",
        // with the final one carrying no expression.
        String text = body.replace("!", " ").trim();
        String[] segments = text.split(";", -1);

        List<Double> bounds = new ArrayList<>();
        List<String> exprs = new ArrayList<>();

        String head = segments[0].trim();
        int split = indexOfWhitespace(head);
        if (split < 0) {
            throw new IllegalStateException(
                    "SGTE: malformed parameter for " + element + " " + phase + ": " + body);
        }
        bounds.add(Double.parseDouble(head.substring(0, split).trim()));
        exprs.add(head.substring(split).trim().replaceAll("\\s+", ""));

        for (int i = 1; i < segments.length; i++) {
            String seg = segments[i].trim();
            if (seg.isEmpty()) {
                continue;
            }
            // "<T_upper> <flag> [<expr>]"
            String[] parts = seg.split("\\s+", 3);
            bounds.add(Double.parseDouble(parts[0]));
            boolean continues = parts.length > 1 && parts[1].equalsIgnoreCase("Y");
            if (!continues) {
                break;
            }
            if (parts.length < 3) {
                throw new IllegalStateException(
                        "SGTE: range continues but no expression follows for "
                                + element + " " + phase + ": " + seg);
            }
            exprs.add(parts[2].replaceAll("\\s+", ""));
        }

        if (bounds.size() != exprs.size() + 1) {
            throw new IllegalStateException(String.format(
                    "SGTE: %s %s has %d breakpoints but %d expressions -- expected one more"
                            + " breakpoint than expression",
                    element, phase, bounds.size(), exprs.size()));
        }

        double[] bp = new double[bounds.size()];
        for (int i = 0; i < bp.length; i++) {
            bp[i] = bounds.get(i);
        }
        return new Entry(element, phase, bp, exprs.toArray(new String[0]), bp.length);
    }

    private static int indexOfWhitespace(String s) {
        for (int i = 0; i < s.length(); i++) {
            if (Character.isWhitespace(s.charAt(i))) {
                return i;
            }
        }
        return -1;
    }

    // =========================================================================
    // Expression -> coefficient vector
    // =========================================================================

    /**
     * Reduces one SGTE expression to its coefficient vector.
     *
     * <p>Term splitting is on {@code +}/{@code -} at the top level, taking care
     * not to break an exponent ({@code 0.663845E-3}) or a negative power
     * ({@code T**(-1)}). Each term is then classified by its {@code T} factor.
     * Unrecognised powers throw rather than being dropped: a silently ignored
     * term would shift the energy by a plausible-looking amount.</p>
     */
    static double[] parseExpression(String expression) {
        double[] c = new double[NUM_TERMS];
        String s = expression.replace(" ", "");
        if (s.endsWith(";")) {
            s = s.substring(0, s.length() - 1);
        }

        for (String term : splitTerms(s)) {
            if (term.isEmpty()) {
                continue;
            }
            applyTerm(c, term);
        }
        return c;
    }

    /** Splits on top-level +/-, protecting exponent signs and bracketed powers. */
    private static List<String> splitTerms(String s) {
        List<String> terms = new ArrayList<>();
        int depth = 0;
        int start = 0;
        for (int i = 0; i < s.length(); i++) {
            char ch = s.charAt(i);
            if (ch == '(') {
                depth++;
            } else if (ch == ')') {
                depth--;
            } else if ((ch == '+' || ch == '-') && depth == 0 && i > start) {
                char prev = s.charAt(i - 1);
                // "1.5E-3": the sign belongs to the exponent, not a new term.
                // "T**-1": likewise belongs to the power.
                if (prev == 'E' || prev == 'e' || prev == '*') {
                    continue;
                }
                terms.add(s.substring(start, i));
                start = i;
            }
        }
        terms.add(s.substring(start));
        return terms;
    }

    /** Classifies one signed term and accumulates it into the coefficient vector. */
    private static void applyTerm(double[] c, String term) {
        String[] factors = term.split("\\*(?!\\*)"); // split on single '*' only

        if (factors.length == 1) {
            c[0] += Double.parseDouble(stripSign(factors[0]));
            return;
        }

        double value = Double.parseDouble(stripSign(factors[0]));
        String rest = term.substring(factors[0].length() + 1);

        if (rest.equalsIgnoreCase("T")) {
            c[1] += value;
        } else if (rest.toUpperCase().startsWith("T*LN(T)") || rest.toUpperCase().equals("T*LN(T)")) {
            c[2] += value;
        } else if (rest.toUpperCase().startsWith("T**")) {
            int power = parsePower(rest.substring(3));
            c[termIndexForPower(power, term)] += value;
        } else {
            throw new IllegalStateException("SGTE: unrecognised term '" + term + "'");
        }
    }

    private static String stripSign(String s) {
        return s.isEmpty() ? "0" : s;
    }

    private static int parsePower(String raw) {
        String p = raw.replace("(", "").replace(")", "").trim();
        return Integer.parseInt(p);
    }

    /** Maps a power of T onto its slot in the coefficient vector. */
    private static int termIndexForPower(int power, String term) {
        return switch (power) {
            case 2 -> 3;
            case 3 -> 4;
            case 4 -> 5;
            case 7 -> 6;
            case -1 -> 7;
            case -2 -> 8;
            case -3 -> 9;
            case -9 -> 10;
            default -> throw new IllegalStateException(
                    "SGTE: unsupported power T**" + power + " in term '" + term
                            + "'. Add a coefficient slot rather than dropping the term.");
        };
    }

    // =========================================================================
    // Evaluation
    // =========================================================================

    private static double evaluate(double[] c, double t) {
        return c[0]
                + c[1] * t
                + c[2] * t * Math.log(t)
                + c[3] * t * t
                + c[4] * Math.pow(t, 3)
                + c[5] * Math.pow(t, 4)
                + c[6] * Math.pow(t, 7)
                + c[7] / t
                + c[8] * Math.pow(t, -2)
                + c[9] * Math.pow(t, -3)
                + c[10] * Math.pow(t, -9);
    }

    private static double evaluateFirstDerivative(double[] c, double t) {
        return c[1]
                + c[2] * (1.0 + Math.log(t))
                + c[3] * 2 * t
                + c[4] * 3 * t * t
                + c[5] * 4 * Math.pow(t, 3)
                + c[6] * 7 * Math.pow(t, 6)
                - c[7] * Math.pow(t, -2)
                - c[8] * 2 * Math.pow(t, -3)
                - c[9] * 3 * Math.pow(t, -4)
                - c[10] * 9 * Math.pow(t, -10);
    }

    private static double evaluateSecondDerivative(double[] c, double t) {
        return c[2] / t
                + c[3] * 2
                + c[4] * 6 * t
                + c[5] * 12 * t * t
                + c[6] * 42 * Math.pow(t, 5)
                + c[7] * 2 * Math.pow(t, -3)
                + c[8] * 6 * Math.pow(t, -4)
                + c[9] * 12 * Math.pow(t, -5)
                + c[10] * 90 * Math.pow(t, -11);
    }

    /** Clears the parse cache; for tests that swap the database file. */
    static void clearCache() {
        CACHE.clear();
    }

    /** Diagnostic: coefficient vector for one (element, phase, T). */
    static double[] coefficientsFor(String element, String phase, double t) {
        return coefficients(element, phase, t);
    }

    /** Diagnostic: all cached keys. */
    static Map<String, String> cacheSnapshot() {
        Map<String, String> out = new HashMap<>();
        CACHE.forEach((k, v) -> out.put(k, v.numRanges() + " ranges"));
        return out;
    }
}
