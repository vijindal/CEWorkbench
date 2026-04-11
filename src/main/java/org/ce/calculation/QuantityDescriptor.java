package org.ce.calculation;

import org.ce.model.ThermodynamicResult;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

/**
 * Enumeration of all available thermodynamic quantities that can be extracted
 * from a {@link ThermodynamicResult}.
 *
 * <p>Each constant carries metadata (symbol, full name, unit, defaults) and two
 * abstract methods that both the GUI and the CLI use to extract values and test
 * availability without any per-quantity {@code if/else} branches:</p>
 * <ul>
 *   <li>{@link #extract(ThermodynamicResult)} — returns the numeric value</li>
 *   <li>{@link #available(ThermodynamicResult)} — returns {@code false} when the
 *       value is {@code NaN} or otherwise not meaningful for this result</li>
 * </ul>
 *
 * <p><b>Extending the list:</b> add a new enum constant with its metadata and
 * implement the two abstract methods. No other class requires changes — the GUI
 * ({@code QuantitiesPanel}, {@code ResultChartPanel}) and the CLI
 * ({@code ResultFormatter}) both iterate {@link #values()}.</p>
 *
 * <p>The inner {@link SelectionModel} is an observable {@link EnumSet} that is
 * created once in {@code MainWindow} and shared by all calculation-mode panels
 * and the output chart.</p>
 */
public enum QuantityDescriptor {

    GIBBS_ENERGY("G", "Gibbs Energy", "J/mol", true, true) {
        @Override public double  extract(ThermodynamicResult r) { return r.gibbsEnergy; }
        @Override public boolean available(ThermodynamicResult r) { return r.isFreeEnergyValid(); }
    },

    ENTHALPY("H", "Enthalpy", "J/mol", true, true) {
        @Override public double  extract(ThermodynamicResult r) { return r.enthalpy; }
        @Override public boolean available(ThermodynamicResult r) { return !Double.isNaN(r.enthalpy); }
    },

    ENTROPY("S", "Entropy", "J/mol\u00B7K", false, true) {
        // CVM produces S; MCS result has NaN entropy
        @Override public double  extract(ThermodynamicResult r) { return r.entropy; }
        @Override public boolean available(ThermodynamicResult r) { return !Double.isNaN(r.entropy); }
    },

    HEAT_CAPACITY("Cp", "Heat Capacity", "J/mol\u00B7K", false, false) {
        // MCS only; CVM result has NaN heatCapacity
        @Override public double  extract(ThermodynamicResult r) { return r.heatCapacity; }
        @Override public boolean available(ThermodynamicResult r) { return !Double.isNaN(r.heatCapacity); }
    },

    STD_ENTHALPY("\u03C3H", "Std. Enthalpy", "J/mol", false, false) {
        // MCS only; CVM result has NaN stdEnthalpy
        @Override public double  extract(ThermodynamicResult r) { return r.stdEnthalpy; }
        @Override public boolean available(ThermodynamicResult r) { return !Double.isNaN(r.stdEnthalpy); }
    };

    // ── metadata ──────────────────────────────────────────────────────────────

    /** Short symbol used in axis labels, legend, table headers (e.g. "G", "σH"). */
    public final String symbol;

    /** Full human-readable name (e.g. "Gibbs Energy"). */
    public final String displayName;

    /** SI unit string (e.g. "J/mol", "J/mol·K"). */
    public final String unit;

    /**
     * Whether this quantity is selected by default in a fresh {@link SelectionModel}.
     * G and H are true; S, Cp, σH are false.
     */
    public final boolean selectedByDefault;

    /**
     * Whether CVM produces this quantity ({@code true}) or only MCS does
     * ({@code false}).  Used by {@code QuantitiesPanel} to gray out checkboxes
     * that are incompatible with the active engine.
     */
    public final boolean cvmProduced;

    // ── constructor ───────────────────────────────────────────────────────────

    QuantityDescriptor(String symbol, String displayName, String unit,
                       boolean selectedByDefault, boolean cvmProduced) {
        this.symbol           = symbol;
        this.displayName      = displayName;
        this.unit             = unit;
        this.selectedByDefault = selectedByDefault;
        this.cvmProduced      = cvmProduced;
    }

    // ── abstract accessors ────────────────────────────────────────────────────

    /** Extracts the numeric value from {@code r}. May return {@code NaN}. */
    public abstract double  extract(ThermodynamicResult r);

    /**
     * Returns {@code true} when this quantity has a meaningful value in {@code r}
     * (i.e. not {@code NaN} and not otherwise invalid).
     */
    public abstract boolean available(ThermodynamicResult r);

    // =========================================================================
    // SelectionModel — observable EnumSet shared by all GUI panels + chart
    // =========================================================================

    /**
     * Observable selection model that tracks which {@link QuantityDescriptor}s
     * the user has chosen to display.
     *
     * <p>Create one instance in {@code MainWindow} and pass it to every panel
     * and to {@code OutputPanel}.  Registered {@link Runnable} listeners fire
     * whenever the selection changes (e.g. a checkbox is toggled), allowing the
     * chart to repaint immediately.</p>
     */
    public static class SelectionModel {

        private final EnumSet<QuantityDescriptor> selected = EnumSet.noneOf(QuantityDescriptor.class);
        private final List<Runnable> listeners = new ArrayList<>();

        /** Initialises selection to all quantities where {@link #selectedByDefault} is true. */
        public SelectionModel() {
            for (QuantityDescriptor q : values()) {
                if (q.selectedByDefault) selected.add(q);
            }
        }

        /** Sets whether {@code q} is selected, then notifies all listeners. */
        public void set(QuantityDescriptor q, boolean on) {
            if (on) selected.add(q); else selected.remove(q);
            listeners.forEach(Runnable::run);
        }

        /** Clears existing selection and sets {@code q} exclusively. */
        public void setExclusive(QuantityDescriptor q) {
            selected.clear();
            selected.add(q);
            listeners.forEach(Runnable::run);
        }

        /** Returns {@code true} when {@code q} is currently selected. */
        public boolean isSelected(QuantityDescriptor q) {
            return selected.contains(q);
        }

        /**
         * Returns a snapshot copy of the currently selected quantities in
         * declaration order.
         */
        public EnumSet<QuantityDescriptor> selected() {
            return selected.isEmpty() ? EnumSet.noneOf(QuantityDescriptor.class)
                                      : EnumSet.copyOf(selected);
        }

        /** Registers a listener that is called whenever the selection changes. */
        public void addChangeListener(Runnable r) {
            listeners.add(r);
        }
    }
}
