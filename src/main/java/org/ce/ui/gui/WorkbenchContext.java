package org.ce.ui.gui;

import org.ce.model.ModelSession;
import org.ce.model.storage.Workspace.SystemId;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Shared, observable session context for the CE Workbench GUI.
 *
 * <p>Holds the currently selected system identity (elements, structure, model)
 * and the active {@link ModelSession} derived from it. All panels subscribe
 * via {@link #addChangeListener(Runnable)} for identity changes and via
 * {@link #addSessionListener(Consumer)} for session availability changes.</p>
 *
 * <p>This is the single source of truth for both system identity and model
 * state within the GUI session.</p>
 */
public class WorkbenchContext {

    private SystemId system;
    private ModelSession activeSession;
    private boolean buildingSession = false;

    private final List<Runnable>                systemListeners  = new ArrayList<>();
    private final List<Consumer<ModelSession>>  sessionListeners = new ArrayList<>();

    // =========================================================================
    // System identity
    // =========================================================================

    /**
     * Updates the current system identity and notifies all listeners.
     * Ignores the call if any field is blank. Also invalidates the active session
     * because the model state is no longer valid for the new identity.
     */
    public void setSystem(String elements, String structure, String model) {
        if (elements == null || structure == null || model == null) return;
        if (elements.isBlank() || structure.isBlank() || model.isBlank()) return;
        try {
            this.system = new SystemId(elements, structure, model);
            invalidateSession();
            systemListeners.forEach(Runnable::run);
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
        systemListeners.add(r);
    }

    // =========================================================================
    // Model session lifecycle
    // =========================================================================

    /**
     * Sets the active session (called after a successful async build).
     * Fires all session listeners with the new session.
     */
    public void setActiveSession(ModelSession session) {
        this.activeSession = session;
        sessionListeners.forEach(l -> l.accept(session));
    }

    /** Marks whether a session build is currently in progress. */
    public void setBuildingSession(boolean b) { this.buildingSession = b; }

    /** Returns {@code true} if a session build is currently running. */
    public boolean isBuildingSession() { return buildingSession; }

    /**
     * Invalidates the current session (sets it to {@code null}) and notifies
     * session listeners. Called automatically when system identity changes.
     */
    public void invalidateSession() {
        this.activeSession = null;
        sessionListeners.forEach(l -> l.accept(null));
    }

    /** Returns the active session, or {@code null} if none has been built yet. */
    public ModelSession getActiveSession() {
        return activeSession;
    }

    /** Returns {@code true} if a session is ready for calculations. */
    public boolean hasActiveSession() {
        return activeSession != null;
    }

    /**
     * Registers a listener that fires whenever the session changes.
     * The argument is the new session, or {@code null} if it was invalidated.
     */
    public void addSessionListener(Consumer<ModelSession> listener) {
        sessionListeners.add(listener);
    }
}
