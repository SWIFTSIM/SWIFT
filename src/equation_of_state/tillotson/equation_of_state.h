/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_TILLOTSON_EQUATION_OF_STATE_H
#define SWIFT_TILLOTSON_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/tillotson/equation_of_state.h
 *
 * Only P(rho, u) and c_s(rho, u) are implemented for now!
 * So, must be used with the Minimal SPH formulation.
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"

extern struct eos_parameters eos;
/* ------------------------------------------------------------------------- */

// Tillotson parameters
struct Til_params {
    float rho_0, a, b, A, B, E_0, E_iv, E_cv, alpha, beta, eta_min;
};

struct eos_parameters {
    struct Til_params Til_iron, Til_granite, Til_water;
};

// Material identifier flags
enum material_id {
    // Tillotson
    Til_iron = 10,
    Til_granite = 11,
    Til_water = 12
};

// Tillotson parameter values for each material
void set_Til_iron(struct Til_params *mat) {
    mat->rho_0 = 7.800;
    mat->a = 0.5;
    mat->b = 1.5;
    mat->A = 1.28e12;
    mat->B = 1.05e12;
    mat->E_0 = 9.5e10;
    mat->E_iv = 2.4e10;
    mat->E_cv = 8.67e10;
    mat->alpha = 5.0;
    mat->beta = 5.0;
    mat->eta_min = 0.0;
}
void set_Til_granite(struct Til_params *mat) {
    mat->rho_0 = 2.680;
    mat->a = 0.5;
    mat->b = 1.3;
    mat->A = 1.8e11;
    mat->B = 1.8e11;
    mat->E_0 = 1.6e11;
    mat->E_iv = 3.5e10;
    mat->E_cv = 1.8e11;
    mat->alpha = 5.0;
    mat->beta = 5.0;
    mat->eta_min = 0.0;
}
void set_Til_water(struct Til_params *mat) {
    mat->rho_0 = 0.998;
    mat->a = 0.7;
    mat->b = 0.15;
    mat->A = 2.18e10;
    mat->B = 1.325e11;
    mat->E_0 = 7.0e10;
    mat->E_iv = 4.19e9;
    mat->E_cv = 2.69e10;
    mat->alpha = 10.0;
    mat->beta = 5.0;
    mat->eta_min = 0.915;
}


/**
 * @brief Returns the internal energy given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float pressure, int mat_id) {

  return 0;
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy, int mat_id) {

  return 0;
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u, int mat_id) {

  return 0;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u, int mat_id) {
    struct Til_params *mat;
    // Select the material parameters
    switch(mat_id) {
        case Til_iron:
            mat = &eos.Til_iron;

        case Til_granite:
            mat = &eos.Til_granite;

        case Til_water:
            mat = &eos.Til_water;
    };

    const float eta = density / params->rho_0;
    const float mu = eta - 1.f;
    const float nu = 1.f/eta - 1.f;
    float P_c, P_e, P;

    // Condensed or cold
    if (eta < eta_min) {
        P_c = 0.f;
    }
    else {
        P_c = (mat->a + mat->b / (u / (mat->E_0 * eta*eta) + 1.f)) * density * u
            + mat->A * mu + mat->B * mu*mu;
    }
    // Expanded and hot
    P_e = mat->a*density*u + (
        mat->b * density * u / (u / (mat->E_0 * eta*eta) + 1.f)
        + mat->A*mu * np.exp(-mat->beta * nu)
        ) * np.exp(-mat->alpha * nu*nu);

    // Condensed or cold state
    if ((1.f < eta) || (u < mat->E_iv)) {
        P = P_c;
    }
    // Expanded and hot state
    else if ((eta < 1.f) && (mat->E_cv < u)) {
        P = P_e;
    }
    // Hybrid state
    else {
        P = ((u - mat->E_iv)*P_e + (mat->E_cv - u)*P_c) / (mat->E_cv - mat->E_iv);
    }

    if (P < 0.f) {
        P = 0.f;
    }

    return P;
}

/**
 * @brief Returns the internal energy given density and pressure.
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The internal energy \f$u\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_pressure(float density, float pressure, int mat_id) {

  return 0;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u, int mat_id) {

  return ; ///wilo
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * NOT IMPLEMENTED!
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P, int mat_id) {

  return 0;
}


/**
 * @brief Initialize the eos parameters
 *
 * @param e The #eos_parameters
 * @param params The parsed parameters
 */
__attribute__((always_inline)) INLINE static void eos_init(
    struct eos_parameters *e, const struct swift_params *params) {

    // Set the Tillotson parameters for each material
    set_Til_iron(e->Til_iron);
    set_Til_granite(e->Til_granite);
    set_Til_water(e->Til_water);
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {

  message("Equation of state: Tillotson.");
}

#if defined(HAVE_HDF5)
/**
 * @brief Write equation of state information to the snapshot
 *
 * @param h_grpsph The HDF5 group in which to write
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print_snapshot(
    hid_t h_grpsph, const struct eos_parameters *e) {

  io_write_attribute_s(h_grpsph, "Equation of state", "Tillotson");
}
#endif

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */
