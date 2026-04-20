package org.ce.model;

/** Physical constants shared across all thermodynamic engines. */
public final class PhysicsConstants {
    private PhysicsConstants() {}

    /** Gas constant (J/mol·K) — CODATA 2014 value. */
    public static final double R_GAS = 8.3144598;

    /** eV per atom to J/mol conversion factor (1 eV * N_A). */
    public static final double EV_TO_J_MOL = 96485.3321;
}
