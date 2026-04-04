package org.ce.ui.gui;

import org.ce.storage.SystemId;

import java.util.ArrayList;
import java.util.List;

/**
 * Shared, observable session context for the CE Workbench GUI.
 *
 * <p>Holds the currently selected system identity (elements, structure, model).
 * All panels that need this information subscribe via {@link #addChangeListener(Runnable)}
 * and write changes via {@link #setSystem(String, String, String)}.</p>
 *
 * <p>This is the single source of truth for system identity within the GUI session —
 * the GUI equivalent of the CLI's positional arguments.</p>
 */
public class WorkbenchContext {

    private SystemId system;
    private final List<Runnable> listeners = new ArrayList<>();

    /**
     * Updates the current system and notifies all listeners.
     * Ignores the call if any field is blank.
     */
    public void setSystem(String elements, String structure, String model) {
        if (elements == null || structure == null || model == null) return;
        if (elements.isBlank() || structure.isBlank() || model.isBlank()) return;
        try {
            this.system = new SystemId(elements, structure, model);
            listeners.forEach(Runnable::run);
        } catch (IllegalArgumentException ignored) {
            // e.g. unsupported ncomp — don't update
        }
    }

    /** Returns the current system, or {@code null} if none has been set. */
    public SystemId getSystem() {
        return system;
    }

    /** Returns {@code true} if a valid system has been set. */
    public boolean hasSystem() {
        return system != null;
    }

    /**
     * Registers a listener that fires (on the calling thread) whenever
     * the system identity changes.
     */
    public void addChangeListener(Runnable r) {
        listeners.add(r);
    }
}
