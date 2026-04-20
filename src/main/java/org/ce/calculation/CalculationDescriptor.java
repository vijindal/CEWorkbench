package org.ce.calculation;

import org.ce.model.ModelSession.EngineConfig;
import java.util.*;

/**
 * Defines the core vocabulary and schemas for the calculation discovery system.
 */
public final class CalculationDescriptor {

    private CalculationDescriptor() {}

    /** The thermodynamic quantity to be calculated. */
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

    /** The dimensionality or "shape" of the calculation. */
    public enum Mode {
        ANALYSIS("Analysis"),
        FINITE_SIZE_SCALING("Finite-Size Scaling (FSS)");

        public final String displayName;
        Mode(String displayName) { this.displayName = displayName; }
    }

    /** Defines a single input requirement for a calculation. */
    public static final class Parameter {
        public static final Parameter TEMPERATURE = new Parameter("Temperature", Double.class, 1000.0);
        public static final Parameter COMPOSITION = new Parameter("Composition", double[].class, null);
        public static final Parameter T_START      = new Parameter("T Start", Double.class, 1000.0);
        public static final Parameter T_END        = new Parameter("T End", Double.class, 1000.0);
        public static final Parameter T_STEP       = new Parameter("T Step", Double.class, 100.0);
        public static final Parameter X_START      = new Parameter("X Start", Double.class, 0.5);
        public static final Parameter X_END        = new Parameter("X End", Double.class, 0.5);
        public static final Parameter X_STEP       = new Parameter("X Step", Double.class, 0.1);

        public static final Parameter X_STARTS     = new Parameter("X Starts", double[].class, null);
        public static final Parameter X_ENDS       = new Parameter("X Ends", double[].class, null);
        public static final Parameter X_STEPS      = new Parameter("X Steps", double[].class, null);

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
        @Override public String toString() { return name; }
    }

    /** Carries global system identity (elements, structure, model). */
    public record ModelSpecifications(String elements, String structure, String modelName, EngineConfig engineConfig) {
        @Override public String toString() {
            return String.format("%s / %s / %s [%s]", elements, structure, modelName, engineConfig);
        }
    }

    /** Value object representing the specifications for a single calculation job. */
    public static final class JobSpecifications {
        private final Property property;
        private final Mode mode;
        private final Map<Parameter, Object> parameters = new HashMap<>();

        public JobSpecifications(Property property, Mode mode) {
            this.property = Objects.requireNonNull(property);
            this.mode = Objects.requireNonNull(mode);
        }

        public void set(Parameter param, Object value) { parameters.put(param, value); }

        @SuppressWarnings("unchecked")
        public <T> Optional<T> get(Parameter param) { return Optional.ofNullable((T) parameters.get(param)); }

        @SuppressWarnings("unchecked")
        public <T> T getOrDefault(Parameter param) { return (T) parameters.getOrDefault(param, param.defaultValue); }

        public Property getProperty() { return property; }
        public Mode getMode() { return mode; }

        @Override public String toString() {
            return String.format("Request[%s in %s mode | %d params]", property, mode, parameters.size());
        }
    }

    /** Metadata provider for discoverable properties and requirements. */
    public static final class Registry {
        public static List<Property> getAvailableProperties(EngineConfig engine) {
            if (engine.isCvm()) return Arrays.asList(Property.GIBBS_ENERGY, Property.ENTHALPY, Property.ENTROPY);
            return Arrays.asList(Property.ENTHALPY, Property.HEAT_CAPACITY, Property.CORRELATION_FUNCTIONS);
        }

        public static List<Mode> getAvailableModes(Property property, EngineConfig engine) {
            if (property == Property.HEAT_CAPACITY) return Arrays.asList(Mode.FINITE_SIZE_SCALING);
            return Arrays.asList(Mode.ANALYSIS);
        }

        public static List<Parameter> getRequirements(Property property, Mode mode, EngineConfig engine) {
            List<Parameter> requirements = new ArrayList<>();
            switch (mode) {
                case ANALYSIS:
                    requirements.addAll(Arrays.asList(Parameter.T_START, Parameter.T_END, Parameter.T_STEP,
                                       Parameter.X_STARTS, Parameter.X_ENDS, Parameter.X_STEPS));
                    break;
                case FINITE_SIZE_SCALING:
                    requirements.addAll(Arrays.asList(Parameter.TEMPERATURE, Parameter.COMPOSITION));
                    break;
            }
            if (engine.isMcs()) {
                if (mode != Mode.FINITE_SIZE_SCALING) requirements.add(Parameter.MCS_L);
                requirements.addAll(Arrays.asList(Parameter.MCS_NEQUIL, Parameter.MCS_NAVG));
            }
            return requirements;
        }
    }
}
