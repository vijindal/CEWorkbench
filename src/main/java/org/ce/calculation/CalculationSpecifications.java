package org.ce.calculation;

import org.ce.calculation.CalculationDescriptor.Mode;
import org.ce.calculation.CalculationDescriptor.Parameter;
import org.ce.calculation.CalculationDescriptor.Property;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

/**
 * A unified value object representing the specifications for a single calculation job.
 *
 * <p>Carries the target property, the calculation mode, and a map of parameter
 * values. The UI constructs this based on the metadata provided by the
 * {@link CalculationRegistry}.</p>
 */
public final class CalculationSpecifications {

    private final Property property;
    private final Mode mode;
    private final Map<Parameter, Object> parameters = new HashMap<>();

    public CalculationSpecifications(Property property, Mode mode) {
        this.property = Objects.requireNonNull(property);
        this.mode = Objects.requireNonNull(mode);
    }

    /**
     * Sets a parameter value. This should be used by the UI layer after
     * collecting user input.
     */
    public void set(Parameter param, Object value) {
        parameters.put(param, value);
    }

    /**
     * Gets a parameter value, or empty if not set.
     */
    @SuppressWarnings("unchecked")
    public <T> Optional<T> get(Parameter param) {
        return Optional.ofNullable((T) parameters.get(param));
    }

    /**
     * Gets a parameter value, falling back to the parameter's default as defined
     * in {@link CalculationDescriptor}.
     */
    @SuppressWarnings("unchecked")
    public <T> T getOrDefault(Parameter param) {
        return (T) parameters.getOrDefault(param, param.defaultValue);
    }

    public Property getProperty() { return property; }
    public Mode getMode() { return mode; }

    @Override
    public String toString() {
        return String.format("Request[%s in %s mode | %d params]", property, mode, parameters.size());
    }
}
