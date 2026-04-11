package org.ce.ui.gui;

import org.ce.calculation.CalculationDescriptor.Parameter;

import javax.swing.*;
import java.awt.*;

/**
 * Factory for creating Swing input components for specific Calculation Parameters.
 */
public final class ParameterFieldFactory {

    private ParameterFieldFactory() {}

    /**
     * Creates a component for the given parameter.
     * @param parameter The parameter metadata.
     * @param initialValue The initial value to set.
     * @return A JComponent representing the input field.
     */
    public static JComponent createEditor(Parameter parameter, Object initialValue) {
        Object value = (initialValue != null) ? initialValue : parameter.defaultValue;

        if (parameter.type == Double.class) {
            double v = (value != null) ? (Double) value : 0.0;
            return new JSpinner(new SpinnerNumberModel(v, -1e10, 1e10, 0.1));
        }

        if (parameter.type == Integer.class) {
            int v = (value != null) ? (Integer) value : 0;
            return new JSpinner(new SpinnerNumberModel(v, 0, 1000000, 1));
        }

        if (parameter.type == double[].class) {
            JTextField field = new JTextField();
            if (value != null) {
                field.setText(formatArray((double[]) value));
            }
            return field;
        }

        return new JTextField(value != null ? value.toString() : "");
    }

    private static String formatArray(double[] arr) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            sb.append(String.format("%.4f", arr[i]));
            if (i < arr.length - 1) sb.append(", ");
        }
        return sb.toString();
    }

    /**
     * Extracts the current value from a component created by this factory.
     */
    public static Object getValue(Parameter parameter, JComponent component) {
        if (component instanceof JSpinner) {
            return ((JSpinner) component).getValue();
        }

        if (component instanceof JTextField) {
            String text = ((JTextField) component).getText().trim();
            if (parameter.type == double[].class) {
                return parseArray(text);
            }
            return text;
        }

        return null;
    }

    private static double[] parseArray(String text) {
        if (text.isEmpty()) return new double[0];
        String[] parts = text.split("[,\\s]+");
        double[] vals = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            vals[i] = Double.parseDouble(parts[i]);
        }
        return vals;
    }
}
