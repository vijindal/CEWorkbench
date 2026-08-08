package org.ce.calculation;

/** A scan range over one scalar quantity (temperature or one element's mole fraction). */
public record Range(double start, double end, double step) {

    public static Range fixed(double v) {
        return new Range(v, v, 0.0);
    }

    public boolean varies() {
        return Math.abs(start - end) > 1e-9;
    }

    public int pointCount() {
        return varies() ? (int) Math.round((end - start) / step) + 1 : 1;
    }

    public double valueAt(int i) {
        return varies() ? start + i * step : start;
    }
}
