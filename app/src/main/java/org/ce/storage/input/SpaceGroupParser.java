package org.ce.storage.input;

import org.ce.domain.cluster.SpaceGroup;
import org.ce.domain.cluster.SymmetryOperation;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class SpaceGroupParser {

    public static SpaceGroup parseFromResources(String baseName)
            throws Exception {

        // ==========================================
        // 1ï¸âƒ£ Read symmetry operations
        // ==========================================
        InputStream is = SpaceGroupParser.class
                .getClassLoader()
                .getResourceAsStream("sym/" + baseName + ".txt");

        if (is == null)
            throw new RuntimeException("File not found: " + baseName + ".txt");

        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        StringBuilder sb = new StringBuilder();
        String line;

        while ((line = br.readLine()) != null) {
            sb.append(line.trim());
        }

        br.close();

        String content = sb.toString()
                .replaceAll("\\{", "")
                .replaceAll("}", "");

        String[] tokens = content.split(",");

        List<Double> numbers = new ArrayList<>();
        for (String t : tokens) {
            if (!t.trim().isEmpty())
                numbers.add(Double.parseDouble(t.trim()));
        }

        int matrixSize = 12;
        int numOps = numbers.size() / matrixSize;

        List<SymmetryOperation> ops = new ArrayList<>();

        for (int i = 0; i < numOps; i++) {

            double[][] rot = new double[3][3];
            double[] trans = new double[3];

            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 3; c++) {
                    rot[r][c] = numbers.get(i * matrixSize + r * 4 + c);
                }
                trans[r] = numbers.get(i * matrixSize + r * 4 + 3);
            }

            ops.add(new SymmetryOperation(rot, trans));
        }

        // ==========================================
        // 2ï¸âƒ£ Read rotateMat + translateMat
        // ==========================================
        InputStream isMat = SpaceGroupParser.class
                .getClassLoader()
                .getResourceAsStream("sym/" + baseName + "_mat.txt");

        if (isMat == null)
            throw new RuntimeException("File not found: " + baseName + "_mat.txt");

        BufferedReader brMat = new BufferedReader(new InputStreamReader(isMat));

        StringBuilder sbMat = new StringBuilder();

        while ((line = brMat.readLine()) != null) {
            sbMat.append(line.trim());
        }

        brMat.close();

        String matContent = sbMat.toString()
                .replaceAll("\\{", "")
                .replaceAll("}", "");

        String[] matTokens = matContent.split(",");

        List<Double> matNumbers = new ArrayList<>();
        for (String t : matTokens) {
            if (!t.trim().isEmpty())
                matNumbers.add(Double.parseDouble(t.trim()));
        }

        double[][] rotateMat = new double[3][3];
        double[] translateMat = new double[3];

        for (int i = 0; i < 9; i++) {
            rotateMat[i / 3][i % 3] = matNumbers.get(i);
        }

        for (int i = 0; i < 3; i++) {
            translateMat[i] = matNumbers.get(9 + i);
        }

        return new SpaceGroup(
                baseName,
                ops,
                rotateMat,
                translateMat
        );
    }
}



