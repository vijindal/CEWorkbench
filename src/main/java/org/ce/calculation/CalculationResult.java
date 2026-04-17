package org.ce.calculation;

import org.ce.model.ThermodynamicResult;
import java.util.List;

public sealed interface CalculationResult {
    record Grid  (List<List<ThermodynamicResult>> values) implements CalculationResult {}
    record Single(ThermodynamicResult value)               implements CalculationResult {}
}
