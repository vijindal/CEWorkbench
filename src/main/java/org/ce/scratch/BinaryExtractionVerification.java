package org.ce.scratch;

import org.ce.CEWorkbenchContext;
import org.ce.calculation.workflow.BinarySubsystemExtractor;
import org.ce.model.ModelSession;
import org.ce.model.ModelSession.EngineConfig;
import org.ce.model.hamiltonian.CECEntry;
import org.ce.model.storage.Workspace;
import org.ce.model.storage.Workspace.SystemId;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Verifies {@link BinarySubsystemExtractor} pulls every binary subsystem out of
 * the quaternary Hamiltonian, with the values the paper publishes.
 *
 * <pre>
 *   ./gradlew runScratch -PscratchClass=org.ce.scratch.BinaryExtractionVerification
 * </pre>
 *
 * <p>The quaternary edge scan needs all six binaries, so this gate is a
 * prerequisite for it: at a square-plot edge one element is zero and the point
 * is solved as a true binary rather than on the quaternary Hamiltonian.</p>
 *
 * <p><b>Why this needed a fix.</b> The extractor matched only the internal
 * letter form ({@code v21AB}), but {@code Nb-Ti-V-Zr_BCC_A2_T_CVCF} stores
 * positional {@code CF_<index>} names, so all six pairs failed. Names are now
 * canonicalised first, which also admits the paper's own {@code e2NbTi1}
 * spelling.</p>
 *
 * <p><b>Expected values are the paper's, not our file's.</b> They are typed in
 * from Table 17 ("CECs for the Nb-Ti-V-Zr system in the CVCF basis", Jindal
 * &amp; Lele 2025, CALPHAD 89, 102825) rather than read back from the
 * Hamiltonian, so this checks the stored data against the publication instead
 * of checking it against itself.</p>
 */
public final class BinaryExtractionVerification {

    /** Table 17: {pair, e2..1 (1NN), e2..2 (2NN), e3.., e4..}. */
    private static final String[][] TABLE_17 = {
            { "Nb", "Ti",  "6240",    "3120",     "0", "0" },
            { "Nb", "V",  "14080",    "7040",     "0", "0" },
            { "Nb", "Zr",  "7401.6",  "3700.8",  "-4224", "0" },
            { "Ti", "V",   "8160",    "4080",     "0", "0" },
            { "Ti", "Zr", "-9264",   "-4632",     "0", "0" },
            { "V",  "Zr", "17920",    "8960",  "-5760", "0" },
    };

    private static int failures = 0;

    public static void main(String[] args) throws Exception {
        java.util.logging.LogManager.getLogManager().reset();

        System.out.println("=".repeat(74));
        System.out.println("  Binary extraction from Nb-Ti-V-Zr vs the paper's Table 17");
        System.out.println("=".repeat(74));

        CEWorkbenchContext ctx = new CEWorkbenchContext(new Workspace());
        ModelSession session = new ModelSession.Builder(ctx.getHamiltonianStore())
                .build(new SystemId("Nb-Ti-V-Zr", "BCC_A2", "T"), EngineConfig.CVM, null);
        List<String> elements = List.of("Nb", "Ti", "V", "Zr");

        System.out.printf("%n  parent stores its terms as: %s ...%n",
                session.cecEntry.cecTerms[0].name);

        for (String[] row : TABLE_17) {
            String a = row[0], b = row[1];
            String pair = a + "-" + b;
            System.out.printf("%n  --- %s ---%n", pair);

            CECEntry binary;
            try {
                binary = BinarySubsystemExtractor.extractBinary(
                        session.cecEntry, elements, a, b);
            } catch (RuntimeException e) {
                fail(pair + " extraction: " + e.getMessage());
                continue;
            }

            // In the extracted binary the pair is always A,B.
            Map<String, Double> expected = new LinkedHashMap<>();
            expected.put("v21AB", Double.parseDouble(row[2]));
            expected.put("v22AB", Double.parseDouble(row[3]));
            expected.put("v3AB",  Double.parseDouble(row[4]));
            expected.put("v4AB",  Double.parseDouble(row[5]));

            Map<String, Double> got = new LinkedHashMap<>();
            for (CECEntry.CECTerm t : binary.cecTerms) {
                got.put(t.name, t.a);
            }

            for (Map.Entry<String, Double> e : expected.entrySet()) {
                Double v = got.get(e.getKey());
                boolean ok = v != null && Math.abs(v - e.getValue()) < 1e-6;
                System.out.printf("    %-6s %12s   paper %12.1f   %s%n",
                        e.getKey(), v == null ? "(absent)" : String.format("%.1f", v),
                        e.getValue(), ok ? "ok" : "MISMATCH");
                if (!ok) {
                    failures++;
                }
            }

            // A binary Hamiltonian must not carry higher-order interactions.
            for (String name : got.keySet()) {
                if (!expected.containsKey(name)) {
                    fail(pair + " leaked a non-binary term: " + name);
                }
            }
            check(pair + " elements are " + a + "-" + b,
                    binary.elements.equalsIgnoreCase(a + "-" + b));
        }

        System.out.println("\n" + "=".repeat(74));
        System.out.printf("RESULT: %s   (%d failures)%n", failures == 0 ? "PASS" : "FAIL", failures);
        System.out.println("=".repeat(74));
        if (failures > 0) {
            throw new AssertionError(failures + " binary extraction checks failed");
        }
    }

    private static void check(String what, boolean ok) {
        if (!ok) {
            fail(what);
        }
    }

    private static void fail(String what) {
        failures++;
        System.out.println("    [!] FAIL  " + what);
    }

    private BinaryExtractionVerification() {
    }
}
