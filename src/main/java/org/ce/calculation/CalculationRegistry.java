package org.ce.calculation;

import org.ce.calculation.CalculationDescriptor.Mode;
import org.ce.calculation.CalculationDescriptor.Parameter;
import org.ce.calculation.CalculationDescriptor.Property;
import org.ce.model.ModelSession.EngineConfig;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Metadata provider that UI uses to discover available properties and modes.
 *
 * <p>Acts as the "Source of Truth" for what can be calculated for a given engine
 * and what parameters are required. This ensures the UI remains decoupled from
 * thermodynamic logic.</p>
 */
public final class CalculationRegistry {

    public static final CalculationRegistry INSTANCE = new CalculationRegistry();

    private CalculationRegistry() {}

    /**
     * Returns properties supported by the given engine.
     */
    public List<Property> getAvailableProperties(EngineConfig engine) {
        if (engine.isCvm()) {
            return Arrays.asList(Property.GIBBS_ENERGY, Property.ENTHALPY, Property.ENTROPY);
        } else {
            // MCS supports H, CF, and Heat Capacity (via FSS)
            return Arrays.asList(Property.ENTHALPY, Property.HEAT_CAPACITY, Property.CORRELATION_FUNCTIONS);
        }
    }

    /**
     * Returns modes supported by the property and engine.
     */
    public List<Mode> getAvailableModes(Property property, EngineConfig engine) {
        if (property == Property.HEAT_CAPACITY) {
            // Heat capacity currently only via FSS in MCS
            return Arrays.asList(Mode.FINITE_SIZE_SCALING);
        }
        
        List<Mode> modes = new ArrayList<>(Arrays.asList(Mode.SINGLE_POINT, Mode.LINE_SCAN));
        if (engine.isCvm()) {
            modes.add(Mode.GRID_SCAN); // Maps currently CVM only
        }
        return modes;
    }

    /**
     * Returns the list of parameters required for a specific combination.
     */
    public List<Parameter> getRequirements(Property property, Mode mode, EngineConfig engine) {
        List<Parameter> requirements = new ArrayList<>();

        // Basic conditions
        switch (mode) {
            case SINGLE_POINT:
            case FINITE_SIZE_SCALING:
                requirements.add(Parameter.TEMPERATURE);
                requirements.add(Parameter.COMPOSITION);
                break;
            case LINE_SCAN:
                // UI handles selection between T-scan and x-scan internally for now,
                // but let's provide all common scan params
                requirements.add(Parameter.T_START);
                requirements.add(Parameter.T_END);
                requirements.add(Parameter.T_STEP);
                requirements.add(Parameter.X_START);
                requirements.add(Parameter.X_END);
                requirements.add(Parameter.X_STEP);
                break;
            case GRID_SCAN:
                requirements.add(Parameter.T_START);
                requirements.add(Parameter.T_END);
                requirements.add(Parameter.T_STEP);
                requirements.add(Parameter.X_START);
                requirements.add(Parameter.X_END);
                requirements.add(Parameter.X_STEP);
                break;
        }

        // Engine-specific additions (MCS params)
        if (engine.isMcs()) {
            if (mode == Mode.FINITE_SIZE_SCALING) {
                // FSS uses hardcoded L values (12, 16, 24) but needs sweeps
                requirements.add(Parameter.MCS_NEQUIL);
                requirements.add(Parameter.MCS_NAVG);
            } else {
                requirements.add(Parameter.MCS_L);
                requirements.add(Parameter.MCS_NEQUIL);
                requirements.add(Parameter.MCS_NAVG);
            }
        }

        return requirements;
    }
}
