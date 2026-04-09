package org.ce.ui.gui;

import org.ce.calculation.QuantityDescriptor;

import javax.swing.*;
import java.awt.*;
import java.util.EnumMap;
import java.util.Map;

/**
 * Compact checkbox strip for thermodynamic quantity selection.
 *
 * <p>Embedded at the top of each calculation-mode parameter card (Single Point,
 * Line Scan, Map Scan). All panels share the same {@link QuantityDescriptor.SelectionModel}
 * instance so that quantity choices are preserved when switching modes.</p>
 *
 * <p>Engine-aware disabling: S (entropy) is CVM-only and is disabled when the
 * active session uses MCS. Cp and σH are MCS-only and are disabled when the
 * active session uses CVM (or when no session is active).</p>
 */
public class QuantitiesPanel extends JPanel {

    private final QuantityDescriptor.SelectionModel selectionModel;
    private final Map<QuantityDescriptor, JCheckBox> checkboxes = new EnumMap<>(QuantityDescriptor.class);

    public QuantitiesPanel(QuantityDescriptor.SelectionModel selectionModel,
                           WorkbenchContext context) {
        this.selectionModel = selectionModel;

        setOpaque(false);
        setLayout(new BorderLayout(0, 4));

        // Section title
        JLabel title = new JLabel("QUANTITIES");
        title.setFont(title.getFont().deriveFont(Font.BOLD, 9f));
        title.setForeground(new Color(0x858585));
        add(title, BorderLayout.NORTH);

        // Checkbox grid — 3 columns
        JPanel grid = new JPanel(new GridLayout(0, 3, 4, 2));
        grid.setOpaque(false);

        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            JCheckBox cb = new JCheckBox(q.symbol);
            cb.setSelected(selectionModel.isSelected(q));
            cb.setOpaque(false);
            cb.setForeground(new Color(0xD4D4D4));
            cb.setFont(cb.getFont().deriveFont(11f));
            cb.setToolTipText(q.displayName + "  (" + q.unit + ")");
            cb.addActionListener(e -> selectionModel.set(q, cb.isSelected()));
            checkboxes.put(q, cb);
            grid.add(cb);
        }
        add(grid, BorderLayout.CENTER);

        // Wire session listener for engine-aware disabling
        context.addSessionListener(session -> {
            boolean isCvm = (session == null) || session.engineConfig.isCvm();
            updateCheckboxStates(isCvm);
        });

        // Initial state — assume CVM (no session yet)
        updateCheckboxStates(true);

        // Sync checkboxes when selection model changes externally
        selectionModel.addChangeListener(this::syncFromModel);
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    private void updateCheckboxStates(boolean isCvm) {
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            JCheckBox cb = checkboxes.get(q);
            boolean enabled = isCompatible(q, isCvm);
            cb.setEnabled(enabled);
            if (!enabled) {
                selectionModel.set(q, false);
                cb.setSelected(false);
            }
        }
    }

    /** Returns true when quantity {@code q} is compatible with the current engine. */
    private static boolean isCompatible(QuantityDescriptor q, boolean isCvm) {
        if (q == QuantityDescriptor.ENTROPY && !isCvm) return false;      // S: CVM only
        if (q == QuantityDescriptor.HEAT_CAPACITY && isCvm) return false; // Cp: MCS only
        if (q == QuantityDescriptor.STD_ENTHALPY  && isCvm) return false; // σH: MCS only
        return true;
    }

    /** Syncs checkbox ticks from the model (e.g. after external programmatic change). */
    private void syncFromModel() {
        for (QuantityDescriptor q : QuantityDescriptor.values()) {
            JCheckBox cb = checkboxes.get(q);
            boolean sel = selectionModel.isSelected(q);
            if (cb.isSelected() != sel) cb.setSelected(sel);
        }
    }

    // =========================================================================
    // Preferred size hint — compact strip
    // =========================================================================

    @Override
    public Dimension getPreferredSize() {
        Dimension d = super.getPreferredSize();
        d.width = Math.max(d.width, 200);
        return d;
    }
}
