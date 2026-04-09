package org.ce.domain.engine.mcs;

import java.util.ArrayDeque;
import java.util.Deque;

/** Efficient rolling window statistics calculator. */
final class RollingWindow {

    private final Deque<Double> window;
    private final int maxSize;
    private double sum = 0.0;
    private double sumSquares = 0.0;

    RollingWindow(int maxSize) {
        if (maxSize <= 0) throw new IllegalArgumentException("maxSize must be positive");
        this.maxSize = maxSize;
        this.window  = new ArrayDeque<>();
    }

    void add(double value) {
        if (window.size() >= maxSize) {
            double removed = window.removeFirst();
            sum -= removed;
            sumSquares -= removed * removed;
        }
        window.addLast(value);
        sum += value;
        sumSquares += value * value;
    }

    double getMean() {
        if (window.isEmpty()) return 0.0;
        return sum / window.size();
    }

    double getStdDev() {
        if (window.size() < 2) return 0.0;
        double mean     = getMean();
        double variance = (sumSquares / window.size()) - (mean * mean);
        if (variance < 0) variance = 0;
        return Math.sqrt(variance);
    }

    void clear() {
        window.clear();
        sum = 0.0;
        sumSquares = 0.0;
    }
}
