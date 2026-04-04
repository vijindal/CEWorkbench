import java.util.Arrays;

public class TernaryDisorderedVerification {
    public static void main(String[] args) {
        double xA = 0.33, xB = 0.33, xC = 0.34;
        
        // CFs (VJ fallback rules)
        double v4AB = xA*xA*xB*xB;
        double v4AC = xA*xA*xC*xC;
        double v4BC = xB*xB*xC*xC;
        double v4ABC1 = xA*xA*xB*xC;
        double v4ABC2 = xA*xB*xB*xC;
        double v4ABC3 = xA*xB*xC*xC;
        double v3AB = -xA*xA*xB + xA*xB*xB;
        double v3AC = -xA*xA*xC + xA*xC*xC;
        double v3BC = -xB*xB*xC + xB*xC*xC;
        double v3ABC1 = xA*xB*xC;
        double v3ABC2 = xA*xB*xC;
        double v3ABC3 = xA*xB*xC;
        double v22AB = xA*xB;
        double v22AC = xA*xC;
        double v22BC = xB*xC;
        double v21AB = xA*xB;
        double v21AC = xA*xC;
        double v21BC = xB*xC;

        System.out.println("--- Java Code Result (CFs) ---");
        System.out.printf("v4AB: %.7f, v4AC: %.7f, v4BC: %.7f%n", v4AB, v4AC, v4BC);
        System.out.printf("v4ABC1: %.7f, v4ABC2: %.7f, v4ABC3: %.7f%n", v4ABC1, v4ABC2, v4ABC3);
        System.out.printf("v3AB: %.7f, v3AC: %.7f, v3BC: %.7f%n", v3AB, v3AC, v3BC);
        System.out.printf("v3ABC1: %.7f, v3ABC2: %.7f, v3ABC3: %.7f%n", v3ABC1, v3ABC2, v3ABC3);
        System.out.printf("v22AB: %.7f, v22AC: %.7f, v22BC: %.7f%n", v22AB, v22AC, v22BC);
        System.out.printf("v21AB: %.7f, v21AC: %.7f, v21BC: %.7f%n", v21AB, v21AC, v21BC);
        
        System.out.println("\n--- Java Code Result (CV Probabilities) ---");
        // Tetra (P(A,A,A,A), P(A,A,A,B), P(A,A,B,B), P(A,A,C,C), P(A,B,C,D) types)
        System.out.println("Tetra sample: P(A,A,A,A)=" + (xA*xA*xA*xA));
        System.out.println("Tetra sample: P(A,A,B,B)=" + (xA*xA*xB*xB));
        System.out.println("Tetra sample: P(A,A,C,C)=" + (xA*xA*xC*xC));
        System.out.println("Tetra sample: P(A,A,B,C)=" + (xA*xA*xB*xC));
        System.out.println("Tetra sample: P(C,C,C,C)=" + (xC*xC*xC*xC));
        
        System.out.println("Tri sample: P(A,A,A)=" + (xA*xA*xA));
        System.out.println("Tri sample: P(A,A,C)=" + (xA*xA*xC));
        System.out.println("Tri sample: P(C,C,C)=" + (xC*xC*xC));

        System.out.println("Pair sample: P(A,A)=" + (xA*xA));
        System.out.println("Pair sample: P(A,C)=" + (xA*xC));
        System.out.println("Pair sample: P(C,C)=" + (xC*xC));
    }
}
