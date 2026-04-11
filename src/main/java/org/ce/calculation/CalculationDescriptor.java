package org.ce.calculation;

import org.ce.model.ModelSession.EngineConfig;

import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * Defines the core vocabulary and schemas for the calculation discovery system.
 *
 * <p>This class contains the enums and types that the UI uses to describe what
 * it wants to calculate (Property, Mode) and what parameters it must provide.</p>
 */
public final class CalculationDescriptor {

    private CalculationDescriptor() {}

    /**
     * The thermodynamic quantity to be calculated.
     */
    public enum Property {
        GIBBS_ENERGY("G", "Gibbs Energy", "J/mol"),
        ENTHALPY("H", "Enthalpy", "J/mol"),
        ENTROPY("S", "Entropy", "J/mol\u00B7K"),
        HEAT_CAPACITY("Cp", "Heat Capacity", "J/mol\u00B7K"),
        CORRELATION_FUNCTIONS("CF", "Correlation Functions", "");

        public final String symbol;
        public final String displayName;
        public final String unit;

        Property(String symbol, String displayName, String unit) {
            this.symbol = symbol;
            this.displayName = displayName;
            this.unit = unit;
        }
    }

    /**
     * The dimensionality or "shape" of the calculation.
     */
    public enum Mode {
        SINGLE_POINT("Single Point"),
        LINE_SCAN("Line Scan"),
        GRID_SCAN("Grid Scan / Map"),
        FINITE_SIZE_SCALING("Finite-Size Scaling (FSS)");

        public final String displayName;

        Mode(String displayName) {
            this.displayName = displayName;
        }
    }

    /**
     * Defines a single input requirement for a calculation.
     */
    public static final class Parameter {
        public static final Parameter TEMPERATURE = new Parameter("Temperature", Double.class, 1000.0);
        public static final Parameter COMPOSITION = new Parameter("Composition", double[].class, null);
        public static final Parameter T_START      = new Parameter("T Start", Double.class, 500.0);
        public static final Parameter T_END        = new Parameter("T End", Double.class, 1500.0);
        public static final Parameter T_STEP       = new Parameter("T Step", Double.class, 50.0);
        public static final Parameter X_START      = new Parameter("X Start", Double.class, 0.05);
        public static final Parameter X_END        = new Parameter("X End", Double.class, 0.95);
        public static final Parameter X_STEP       = new Parameter("X Step", Double.class, 0.05);
        public static final Parameter MCS_L        = new Parameter("Lattice Size L", Integer.class, 4);
        public static final Parameter MCS_NEQUIL   = new Parameter("Equil. Sweeps", Integer.class, 1000);
        public static final Parameter MCS_NAVG     = new Parameter("Avg. Sweeps", Integer.class, 2000);
        public static final Parameter FIXED_CORRELATIONS = new Parameter("Fixed Correlations", double[].class, null);

        public final String name;
        public final Class<?> type;
        public final Object defaultValue;

        public Parameter(String name, Class<?> type, Object defaultValue) {
            this.name = name;
            this.type = type;
            this.defaultValue = defaultValue;
        }

        @Override
        public String toString() { return name; }
    }

    /**
     * Carries global system identity from the UI to the calculation layer.
     * Use this to construct or retrieve a ModelSession.
     */
    public static final class ModelSpecifications {
        public final String elements;
        public final String structure;
        public final String modelName;
        public final EngineConfig engineConfig;

        public ModelSpecifications(String elements, String structure, String modelName, EngineConfig engineConfig) {
            this.elements = Objects.requireNonNull(elements);
            this.structure = Objects.requireNonNull(structure);
            this.modelName = Objects.requireNonNull(modelName);
            this.engineConfig = Objects.requireNonNull(engineConfig);
        }

        @Override
        public String toString() {
            return String.format("%s / %s / %s [%s]", elements, structure, modelName, engineConfig);
        }
    }
}
