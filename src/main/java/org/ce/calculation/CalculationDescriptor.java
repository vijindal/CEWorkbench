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
        ANALYSIS("Analysis");

        public final String displayName;
        Mode(String displayName) { this.displayName = displayName; }
    }

    /** Defines a single input requirement for a calculation. */
    public static final class Parameter {
        /** Single-point conditions (temperature + composition). See {@link org.ce.calculation.Conditions}. */
        public static final Parameter COMPOSITION = new Parameter("Composition", org.ce.calculation.Conditions.class, null);
        /** Conditions scan (temperature range and/or one composition axis). See {@link org.ce.calculation.ConditionsScan}. */
        public static final Parameter CONDITIONS_SCAN = new Parameter("Conditions Scan", org.ce.calculation.ConditionsScan.class, null);

        // GUI spinner-editor templates only — never read as job parameters directly.
        public static final Parameter T_START      = new Parameter("T Start", Double.class, 1000.0);
        public static final Parameter T_END        = new Parameter("T End", Double.class, 1000.0);
        public static final Parameter T_STEP       = new Parameter("T Step", Double.class, 100.0);
        public static final Parameter X_START      = new Parameter("X Start", Double.class, 0.5);
        public static final Parameter X_END        = new Parameter("X End", Double.class, 0.5);
        public static final Parameter X_STEP       = new Parameter("X Step", Double.class, 0.1);

        public static final Parameter MCS_L        = new Parameter("Lattice Size L", Integer.class, 4);
        public static final Parameter MCS_NEQUIL   = new Parameter("Equil. Sweeps", Integer.class, 100);
        public static final Parameter MCS_NAVG     = new Parameter("Avg. Sweeps", Integer.class, 500);
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

        public void set(Parameter param, Object value) {
            if (value != null && !param.type.isInstance(value))
                throw new IllegalArgumentException(
                        "Parameter '" + param.name + "' expects " + param.type.getSimpleName()
                        + ", got " + value.getClass().getSimpleName());
            parameters.put(param, value);
        }

        @SuppressWarnings("unchecked")
        public <T> Optional<T> get(Parameter param) { return Optional.ofNullable((T) parameters.get(param)); }

        @SuppressWarnings("unchecked")
        public <T> T getOrDefault(Parameter param) { return (T) parameters.getOrDefault(param, param.defaultValue); }

        /** Like {@link #getOrDefault} but throws if neither an explicit value nor a default is present. */
        @SuppressWarnings("unchecked")
        public <T> T require(Parameter param, Class<T> type) {
            Object v = parameters.getOrDefault(param, param.defaultValue);
            if (v == null) throw new IllegalStateException("Required parameter '" + param.name + "' not set");
            return type.cast(v);
        }

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
            return Arrays.asList(Mode.ANALYSIS);
        }

        public static List<Parameter> getRequirements(Property property, Mode mode, EngineConfig engine) {
            List<Parameter> requirements = new ArrayList<>();
            switch (mode) {
                case ANALYSIS:
                    requirements.add(Parameter.CONDITIONS_SCAN);
                    break;
            }
            if (engine.isMcs()) {
                requirements.add(Parameter.MCS_L);
                requirements.addAll(Arrays.asList(Parameter.MCS_NEQUIL, Parameter.MCS_NAVG));
            }
            return requirements;
        }
    }
}
