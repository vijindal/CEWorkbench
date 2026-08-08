package org.ce.calculation;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A conditions scan: temperature range + composition ranges by element symbol.
 * At most one axis total (temperature, or exactly one element) may vary —
 * multi-axis scanning is unsupported (would require {@code CalculationResult.Grid}
 * to become N-D).
 */
public record ConditionsScan(Range temperature, Map<String, Range> composition) {

    public ConditionsScan {
        composition = new LinkedHashMap<>(composition); // preserve insertion order for reproducible grids
        long varyingCount = composition.values().stream().filter(Range::varies).count();
        if (temperature.varies()) varyingCount++;
        if (varyingCount > 1) {
            throw new IllegalArgumentException(
                    "Only one condition axis may vary per scan; found " + varyingCount
                    + " varying (temperature and/or composition). Fix all but one to a single value.");
        }
        composition = Map.copyOf(composition);
    }

    public static ConditionsScan fixedAt(Conditions c) {
        Map<String, Range> r = new LinkedHashMap<>();
        for (var e : c.composition().entrySet()) r.put(e.getKey(), Range.fixed(e.getValue()));
        return new ConditionsScan(Range.fixed(c.temperature()), r);
    }

    public int pointCount() {
        if (temperature.varies()) return temperature.pointCount();
        return composition.values().stream().filter(Range::varies)
                .findFirst().map(Range::pointCount).orElse(1);
    }

    /** Returns the {@link Conditions} for scan index {@code i} along the single varying axis. */
    public Conditions pointAt(int i) {
        double T = temperature.valueAt(temperature.varies() ? i : 0);
        Map<String, Double> comp = new LinkedHashMap<>();
        for (var e : composition.entrySet()) {
            Range r = e.getValue();
            comp.put(e.getKey(), r.valueAt(r.varies() ? i : 0));
        }
        return new Conditions(T, comp);
    }

    /** Fail-fast: resolve the first and last grid points before running anything expensive. */
    public void validateAgainst(List<String> canonicalOrder) {
        pointAt(0).resolveComposition(canonicalOrder);
        pointAt(pointCount() - 1).resolveComposition(canonicalOrder);
    }
}
