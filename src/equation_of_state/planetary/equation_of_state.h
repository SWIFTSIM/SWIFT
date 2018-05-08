/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_EQUATION_OF_STATE_H
#define SWIFT_PLANETARY_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/equation_of_state.h
 *
 *              WORK IN PROGRESS!
 *
 * To contain all of the planetary EOS.
 * See material_id for the available choices.
 *
 * Only P(rho, u), c(rho, u), and c(rho, P) are implemented for now!
 * So, must be used with the MinimalMultiMat SPH formulation.
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "units.h"
#include "physical_constants.h"
#include <string.h>

extern struct eos_parameters eos;
/* ------------------------------------------------------------------------- */

// Tillotson parameters
struct Til_params {
    int mat_id;
    float rho_0, a, b, A, B, E_0, E_iv, E_cv, alpha, beta, eta_min, P_min;
    float c_TEMPORARY;
};

// Hubbard & MacFarlane (1980) Uranus/Neptune table parameters
struct HM80_params {
    int mat_id;
    int num_rho, num_u;
    float log_rho_min, log_rho_max, log_rho_step, inv_log_rho_step, log_u_min,
        log_u_max, log_u_step, inv_log_u_step, bulk_mod;
    float **table_P_rho_u;
};
// Table file names
/// to be read in from the parameter file instead once finished testing...
#define HM80_HHe_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_HHe.txt"
#define HM80_ice_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_ice.txt"
#define HM80_rock_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_roc.txt"

// ANEOS (WIP)
struct ANEOS_params {
    int mat_id;
    int num_rho, num_u;
};

// SESAME (WIP)
struct SESAME_params {
    int mat_id;
};

struct eos_parameters {
    struct Til_params Til_iron, Til_granite, Til_water;
    struct HM80_params HM80_HHe, HM80_ice, HM80_rock;
    struct ANEOS_params ANEOS_iron, MANEOS_forsterite;
    struct SESAME_params SESAME_iron;
};

// Material identifier flags (material_ID = type_ID * type_factor + unit_ID)
#define type_factor 10
enum type_id {
    type_Till   = 1,
    type_HM80   = 2,
    type_ANEOS  = 3,
    type_SESAME = 4
};
enum material_id {
    // Tillotson
    Til_iron    = type_Till*type_factor,
    Til_granite = type_Till*type_factor + 1,
    Til_water   = type_Till*type_factor + 2,
    // Hubbard & MacFarlane (1980) Uranus/Neptune
    HM80_HHe    = type_HM80*type_factor,        // Hydrogen-helium atmosphere
    HM80_ice    = type_HM80*type_factor + 1,    // H20-CH4-NH3 ice mix
    HM80_rock   = type_HM80*type_factor + 2,    // SiO2-MgO-FeS-FeO rock mix
    // ANEOS
    ANEOS_iron          = type_ANEOS*type_factor,
    MANEOS_forsterite   = type_ANEOS*type_factor + 1,
    // SESAME
    SESAME_iron = type_SESAME*type_factor,
};

// Parameter values for each material (cgs units)
// Tillotson
INLINE static void set_Til_iron(struct Til_params *mat) {
    mat->mat_id = Til_iron;
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
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}
INLINE static void set_Til_granite(struct Til_params *mat) {
    mat->mat_id = Til_granite;
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
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}
INLINE static void set_Til_water(struct Til_params *mat) {
    mat->mat_id = Til_water;
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
    mat->P_min = 0.0;

    mat->c_TEMPORARY = 9.4e-4;
}

// Hubbard & MacFarlane (1980)
INLINE static void set_HM80_HHe(struct HM80_params *mat) {
    mat->mat_id = HM80_HHe;
    mat->num_rho = 100;
    mat->num_u = 100;
    mat->log_rho_min = -9.2103404;
    mat->log_rho_max = 1.6094379;
    mat->log_rho_step = 0.1092907;
    mat->log_u_min = 9.2103404;
    mat->log_u_max = 22.3327037;
    mat->log_u_step = 0.1325491;
    mat->bulk_mod = 0;

    mat->inv_log_rho_step = 1.f / mat->log_rho_step;
    mat->inv_log_u_step = 1.f / mat->log_u_step;
}
INLINE static void set_HM80_ice(struct HM80_params *mat) {
    mat->mat_id = HM80_ice;
    mat->num_rho = 200;
    mat->num_u = 200;
    mat->log_rho_min = -6.9077553;
    mat->log_rho_max = 2.7080502;
    mat->log_rho_step = 0.0483206;
    mat->log_u_min = 6.9077553;
    mat->log_u_max = 22.3327037;
    mat->log_u_step = 0.0775123;
    mat->bulk_mod = 2.0e10;

    mat->inv_log_rho_step = 1.f / mat->log_rho_step;
    mat->inv_log_u_step = 1.f / mat->log_u_step;
}
INLINE static void set_HM80_rock(struct HM80_params *mat) {
    mat->mat_id = HM80_rock;
    mat->num_rho = 100;
    mat->num_u = 100;
    mat->log_rho_min = -6.9077553;
    mat->log_rho_max = 2.9957323;
    mat->log_rho_step = 0.1000352;
    mat->log_u_min = 9.2103404;
    mat->log_u_max = 20.7232658;
    mat->log_u_step = 0.1162922;
    mat->bulk_mod = 3.49e11;

    mat->inv_log_rho_step = 1.f / mat->log_rho_step;
    mat->inv_log_u_step = 1.f / mat->log_u_step;
}

// Read the table from file
INLINE static void load_HM80_table(struct HM80_params *mat, char *table_file) {
    // Allocate table memory
    mat->table_P_rho_u = (float **) malloc(mat->num_rho*sizeof(float *));
    for (int i=0; i<mat->num_rho; i++) {
        mat->table_P_rho_u[i] = (float *) malloc(mat->num_u*sizeof(float));
    }

    // Load table contents from file
    FILE *f = fopen(table_file, "r");
    for (int i=0; i<mat->num_rho; i++) {
        for (int j=0; j<mat->num_u; j++) {
            fscanf(f, "%f", &mat->table_P_rho_u[i][j]);
        }
    }
    fclose(f);
}

// ANEOS
INLINE static void set_ANEOS_iron(struct ANEOS_params *mat) {
    mat->mat_id = ANEOS_iron;
    mat->num_rho = 100;
    mat->num_u = 100;
}
INLINE static void set_MANEOS_forsterite(struct ANEOS_params *mat) {
    mat->mat_id = MANEOS_forsterite;
    mat->num_rho = 100;
    mat->num_u = 100;
}

// SESAME
INLINE static void set_SESAME_iron(struct SESAME_params *mat) {
    mat->mat_id = SESAME_iron;
}

// Convert from cgs to internal units
// Tillotson
INLINE static void convert_units_Til(
    struct Til_params *mat, const struct unit_system* us) {

    mat->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
    mat->A /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat->B /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
    mat->E_0 /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->E_iv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->E_cv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
    mat->P_min /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

    mat->c_TEMPORARY /= units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
}

// Hubbard & MacFarlane (1980)
#define Mbar_to_Ba 1e12     // Convert Megabar to Barye
INLINE static void convert_units_HM80(
    struct HM80_params *mat, const struct unit_system* us) {

    mat->log_rho_min -= log(units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
    mat->log_rho_max -= log(units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
    mat->log_rho_step -= log(units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));

    mat->log_u_min -= log(units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
    mat->log_u_max -= log(units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
    mat->log_u_step -= log(units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));

    for (int i=0; i<mat->num_rho; i++) {
        for (int j=0; j<mat->num_u; j++) {
            mat->table_P_rho_u[i][j] *= Mbar_to_Ba /
                units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
        }
    }

    mat->bulk_mod /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
}

// ANEOS
INLINE static void convert_units_ANEOS(
    struct ANEOS_params *mat, const struct unit_system* us) {

}

// SESAME
INLINE static void convert_units_SESAME(
    struct SESAME_params *mat, const struct unit_system* us) {

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
__attribute__((always_inline)) INLINE static float
gas_pressure_from_entropy(float density, float entropy, int mat_id) {

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
__attribute__((always_inline)) INLINE static float
gas_entropy_from_pressure(float density, float pressure, int mat_id) {

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
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_entropy(float density, float entropy, int mat_id) {

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

// gas_pressure_from_internal_energy
// Tillotson
INLINE static float Til_pressure_from_internal_energy(float density, float u,
                                                      struct Til_params *mat) {
    const float eta = density / mat->rho_0;
    const float mu = eta - 1.f;
    const float nu = 1.f/eta - 1.f;
    float P_c, P_e, P;

    // Condensed or cold
    if (eta < mat->eta_min) {
        P_c = 0.f;
    }
    else {
        P_c = (mat->a + mat->b / (u / (mat->E_0 * eta*eta) + 1.f)) * density * u
            + mat->A * mu + mat->B * mu*mu;
    }
    // Expanded and hot
    P_e = mat->a*density*u + (
        mat->b * density * u / (u / (mat->E_0 * eta*eta) + 1.f)
        + mat->A*mu * exp(-mat->beta * nu)
        ) * exp(-mat->alpha * nu*nu);

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
        P = ((u - mat->E_iv)*P_e + (mat->E_cv - u)*P_c) /
            (mat->E_cv - mat->E_iv);
    }

    // Minimum pressure
    if (P < mat->P_min) {
        P = mat->P_min;
    }

    return P;
}

// Hubbard & MacFarlane (1980)
INLINE static float HM80_pressure_from_internal_energy(float density, float u,
                                                       struct HM80_params *mat) {
    float P;

    if (u <= 0) {
        return 0;
    }

    int rho_idx, u_idx;
    float intp_rho, intp_u;
    const float log_rho = log(density);
    const float log_u = log(u);

    // 2D interpolation (linear in log(rho), log(u)) to find P(rho, u)
    rho_idx = floor((log_rho - mat->log_rho_min) * mat->inv_log_rho_step);
    u_idx = floor((log_u - mat->log_u_min) * mat->inv_log_u_step);

    intp_rho = (log_rho - mat->log_rho_min - rho_idx*mat->log_rho_step) *
        mat->inv_log_rho_step;
    intp_u = (log_u - mat->log_u_min - u_idx*mat->log_u_step) *
        mat->inv_log_u_step;

    // Return zero pressure if below the table minimum/a
    // Extrapolate the pressure for low densities
    if (rho_idx < 0) {                      // Too-low rho
        P = exp(log((1-intp_u)*mat->table_P_rho_u[0][u_idx]
                    + intp_u*mat->table_P_rho_u[0][u_idx+1])
                + log_rho - mat->log_rho_min);
        if (u_idx < 0) {                    // and too-low u
            P = 0;
        }
    }
    else if (u_idx < 0) {                   // Too-low u
        P = 0;
    }
    // Return an edge value if above the table maximum/a
    else if (rho_idx >= mat->num_rho-1) {   // Too-high rho
        if (u_idx >= mat->num_u-1) {        // and too-high u
            P = mat->table_P_rho_u[mat->num_rho-1][mat->num_u-1];
        }
        else {
            P = mat->table_P_rho_u[mat->num_rho-1][u_idx];
        }
    }
    else if (u_idx >= mat->num_u-1) {       // Too-high u
        P = mat->table_P_rho_u[rho_idx][mat->num_u-1];
    }
    // Normal interpolation within the table
    else {
        P = (1-intp_rho) * ((1-intp_u)*mat->table_P_rho_u[rho_idx][u_idx] +
                            intp_u*mat->table_P_rho_u[rho_idx][u_idx+1]) +
            intp_rho * ((1-intp_u)*mat->table_P_rho_u[rho_idx+1][u_idx] +
                        intp_u*mat->table_P_rho_u[rho_idx+1][u_idx+1]);
    }

    return P;
}

// ANEOS
INLINE static float ANEOS_pressure_from_internal_energy(float density, float u,
                                                        struct ANEOS_params *mat) {
    float P;

    /// Placeholder
    P = mat->num_rho;

    return P;
}

// SESAME
INLINE static float SESAME_pressure_from_internal_energy(float density, float u,
                                                         struct SESAME_params *mat) {
    float P;

    /// Placeholder
    P = mat->mat_id;

    return P;
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u, int mat_id) {
  float P;

  // Material base type
  switch((int)(mat_id/type_factor)) {

    // Tillotson
    case type_Till:;
        // Select the material parameters
        struct Til_params *mat_Til;
        switch(mat_id) {
          case Til_iron:
              mat_Til = &eos.Til_iron;
              break;

          case Til_granite:
              mat_Til = &eos.Til_granite;
              break;

          case Til_water:
              mat_Til = &eos.Til_water;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_Til = &eos.Til_iron; // Ignored, just here to keep the compiler happy
        };

        P = Til_pressure_from_internal_energy(density, u, mat_Til);

        break;

    // Hubbard & MacFarlane (1980)
    case type_HM80:;
        // Select the material parameters
        struct HM80_params *mat_HM80;
        switch(mat_id) {
          case HM80_HHe:
              mat_HM80 = &eos.HM80_HHe;
              break;

          case HM80_ice:
              mat_HM80 = &eos.HM80_ice;
              break;

          case HM80_rock:
              mat_HM80 = &eos.HM80_rock;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        P = HM80_pressure_from_internal_energy(density, u, mat_HM80);

        break;

    /// WIP
    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
          case ANEOS_iron:
              mat_ANEOS = &eos.ANEOS_iron;
              break;

          case MANEOS_forsterite:
              mat_ANEOS = &eos.MANEOS_forsterite;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        P = ANEOS_pressure_from_internal_energy(density, u, mat_ANEOS);

        break;

    /// WIP
    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
          case SESAME_iron:
              mat_SESAME = &eos.SESAME_iron;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_SESAME = &eos.SESAME_iron; // Ignored, just here to keep the compiler happy
        };

        P = SESAME_pressure_from_internal_energy(density, u, mat_SESAME);

        break;

    default:
        error("Unknown material type! mat_id = %d", mat_id);
        P = 0; // Ignored, just here to keep the compiler happy
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

// gas_soundspeed_from_internal_energy
// Tillotson
INLINE static float Til_soundspeed_from_internal_energy(float density, float u,
                                                        struct Til_params *mat) {
//    const float eta = density / mat->rho_0;
//    const float mu = eta - 1.f;
//    const float nu = 1.f/eta - 1.f;
//    float P_c, P_e, P, c_c, c_e, c;
//
//    // Condensed or cold
//    if (eta < mat->eta_min) {
//        P_c = 0.f;
//    }
//    else {
//        P_c = (mat->a + mat->b / (u / (mat->E_0 * eta*eta) + 1.f)) * density * u
//            + mat->A * mu + mat->B * mu*mu;
//    }
//    c_c = mat->a*u + mat->b*u / ((u / (mat->E_0*eta*eta)+1.f) *
//        (u / (mat->E_0*eta*eta)+1.f)) *
//        (3.f*(u / (mat->E_0*eta*eta)+1.f) - 2.f) +
//        (mat->A + 2.f*mat->B*mu) / mat->rho_0  +  P_c / (rho*rho) *
//        (mat->a*rho + mat->b*rho / ((u / (mat->E_0*eta*eta)+1.f) *
//        (u / (mat->E_0*eta*eta)+1.f)));
//
//    c_c = max(c_c, mat->A / mat->rho_0);
//
//    // Expanded and hot
//    P_e = mat->a*density*u + (
//        mat->b * density * u / (u / (mat->E_0 * eta*eta) + 1.f)
//        + mat->A*mu * exp(-mat->beta * nu)
//        ) * exp(-mat->alpha * nu*nu);
//
//    c_e = (mat->a + mat->b / (u / (mat->E_0*eta*eta)+1.f) *
//        exp(-mat->beta*((1.f - eta)/eta)*((1.f - eta)/eta))
//        + 1.f)*P_e/rho + mat->A/mat->rho_0
//        *exp(-(mat->alpha*((1.f - eta)/eta)+mat->beta *
//        ((1.f - eta)/eta)*((1.f - eta)/eta)))*(1.f+mu/(eta*eta)
//        *(mat->alpha+2.f*mat->beta*((1.f - eta)/eta)-eta)) +
//        mat->b*rho*u/((u / (mat->E_0*eta*eta)+1.f)*
//        (u / (mat->E_0*eta*eta)+1.f)*eta*eta)*
//        exp(-mat->beta*((1.f - eta)/eta)*((1.f - eta)/eta))*
//        (2.f*mat->beta*((1.f - eta)/eta)*(u / (mat->E_0*eta*eta)+1.f) /
//         mat->rho_0 + 1.f/(mat->E_0*rho)*(2.f*u-P_e/rho));
//
//    // Condensed or cold state
//    if ((1.f < eta) || (u < mat->E_iv)) {
//        c = c_c;
//    }
//    // Expanded and hot state
//    else if ((eta < 1.f) && (mat->E_cv < u)) {
//        c = c_e;
//    }
//    // Hybrid state
//    else {
//		c = ((u - mat->E_iv)*c_e + (mat->E_cv - u)*c_c) /
//            (mat->E_cv - mat->E_iv);
//
//        c = max(c_c, mat->A / mat->rho0);
//    }
    float c;

    c = mat->c_TEMPORARY; /// VERY TEMPORARY!!!

    return c;
}

// Hubbard & MacFarlane (1980)
INLINE static float HM80_soundspeed_from_internal_energy(float density, float u,
                                                         struct HM80_params *mat) {
    float c, P;

    // Bulk modulus
    if (mat->bulk_mod != 0) {
        c = sqrt(mat->bulk_mod / density);
    }
    // Ideal gas
    else {
        P = gas_pressure_from_internal_energy(density, u, mat->mat_id);
        c = sqrt(5.f/3.f * P / density);
    }

    return c;
}

// ANEOS
INLINE static float ANEOS_soundspeed_from_internal_energy(float density, float u,
                                                          struct ANEOS_params *mat) {
    float c;

    /// Placeholder
    c = mat->num_rho;

    return c;
}

// SESAME
INLINE static float SESAME_soundspeed_from_internal_energy(float density, float u,
                                                           struct SESAME_params *mat) {
    float c;

    /// Placeholder
    c = mat->mat_id;

    return c;
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u, int mat_id) {
  float c;

  // Material base type
  switch((int)(mat_id/type_factor)) {

    // Tillotson
    case type_Till:;
        // Select the material parameters
        struct Til_params *mat_Til;
        switch(mat_id) {
          case Til_iron:
              mat_Til = &eos.Til_iron;
              break;

          case Til_granite:
              mat_Til = &eos.Til_granite;
              break;

          case Til_water:
              mat_Til = &eos.Til_water;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_Til = &eos.Til_iron; // Ignored, just here to keep the compiler happy
        };

        c = Til_soundspeed_from_internal_energy(density, u, mat_Til);

        break;

    // Hubbard & MacFarlane (1980)
    case type_HM80:;
        // Select the material parameters
        struct HM80_params *mat_HM80;
        switch(mat_id) {
          case HM80_HHe:
              mat_HM80 = &eos.HM80_HHe;
              break;

          case HM80_ice:
              mat_HM80 = &eos.HM80_ice;
              break;

          case HM80_rock:
              mat_HM80 = &eos.HM80_rock;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        c = HM80_soundspeed_from_internal_energy(density, u, mat_HM80);

        break;

    /// WIP
    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
          case ANEOS_iron:
              mat_ANEOS = &eos.ANEOS_iron;
              break;

          case MANEOS_forsterite:
              mat_ANEOS = &eos.MANEOS_forsterite;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        c = ANEOS_soundspeed_from_internal_energy(density, u, mat_ANEOS);

        break;

    /// WIP
    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
          case SESAME_iron:
              mat_SESAME = &eos.SESAME_iron;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_SESAME = &eos.SESAME_iron; // Ignored, just here to keep the compiler happy
        };

        c = SESAME_soundspeed_from_internal_energy(density, u, mat_SESAME);

        break;

    default:
        error("Unknown material type! mat_id = %d", mat_id);
        c = 0; // Ignored, just here to keep the compiler happy
  }

  return c;
}

// gas_soundspeed_from_pressure
// Tillotson
INLINE static float Til_soundspeed_from_pressure(float density, float P,
                                                 struct Til_params *mat) {
    float c;

    c = mat->c_TEMPORARY; /// VERY TEMPORARY!!!

    return c;
}

// Hubbard & MacFarlane (1980)
INLINE static float HM80_soundspeed_from_pressure(float density, float P,
                                                  struct HM80_params *mat) {
    float c;

    // Bulk modulus
    if (mat->bulk_mod != 0) {
        c = sqrt(mat->bulk_mod / density);
    }
    // Ideal gas
    else {
        c = sqrt(5.f/3.f * P / density);
    }

    return c;
}

// ANEOS
INLINE static float ANEOS_soundspeed_from_pressure(float density, float P,
                                                   struct ANEOS_params *mat) {
    float c;

    /// Placeholder
    c = mat->num_rho;

    return c;
}

// SESAME
INLINE static float SESAME_soundspeed_from_pressure(float density, float P,
                                                    struct SESAME_params *mat) {
    float c;

    /// Placeholder
    c = mat->mat_id;

    return c;
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_pressure(float density, float P, int mat_id) {
  float c;

  // Material base type
  switch((int)(mat_id/type_factor)) {

    // Tillotson
    case type_Till:;
        // Select the material parameters
        struct Til_params *mat_Til;
        switch(mat_id) {
          case Til_iron:
              mat_Til = &eos.Til_iron;
              break;

          case Til_granite:
              mat_Til = &eos.Til_granite;
              break;

          case Til_water:
              mat_Til = &eos.Til_water;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_Til = &eos.Til_iron; // Ignored, just here to keep the compiler happy
        };

        c = Til_soundspeed_from_internal_energy(density, P, mat_Til);

        break;

    // Hubbard & MacFarlane (1980)
    case type_HM80:;
        // Select the material parameters
        struct HM80_params *mat_HM80;
        switch(mat_id) {
          case HM80_HHe:
              mat_HM80 = &eos.HM80_HHe;
              break;

          case HM80_ice:
              mat_HM80 = &eos.HM80_ice;
              break;

          case HM80_rock:
              mat_HM80 = &eos.HM80_rock;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        c = HM80_soundspeed_from_pressure(density, P, mat_HM80);

        break;

    /// WIP
    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
          case ANEOS_iron:
              mat_ANEOS = &eos.ANEOS_iron;
              break;

          case MANEOS_forsterite:
              mat_ANEOS = &eos.MANEOS_forsterite;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        c = ANEOS_soundspeed_from_pressure(density, P, mat_ANEOS);

        break;

    /// WIP
    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
          case SESAME_iron:
              mat_SESAME = &eos.SESAME_iron;
              break;

          default:
              error("Unknown material ID! mat_id = %d", mat_id);
              mat_SESAME = &eos.SESAME_iron; // Ignored, just here to keep the compiler happy
        };

        c = SESAME_soundspeed_from_pressure(density, P, mat_SESAME);

        break;

    default:
        error("Unknown material type! mat_id = %d", mat_id);
        c = 0; // Ignored, just here to keep the compiler happy
  }

  return c;
}

/**
 * @brief Initialize the eos parameters
 *
 * @param e The #eos_parameters
 * @param params The parsed parameters
 */
__attribute__((always_inline)) INLINE static void eos_init(
    struct eos_parameters *e, const struct phys_const *phys_const,
    const struct unit_system *us, const struct swift_params *params) {

    // Set the parameters and load tables etc. for each material
    // Tillotson
    set_Til_iron(&e->Til_iron);
    set_Til_granite(&e->Til_granite);
    set_Til_water(&e->Til_water);

    // Hubbard & MacFarlane (1980)
    set_HM80_HHe(&e->HM80_HHe);
    set_HM80_ice(&e->HM80_ice);
    set_HM80_rock(&e->HM80_rock);

    load_HM80_table(&e->HM80_HHe, HM80_HHe_table_file);
    load_HM80_table(&e->HM80_ice, HM80_ice_table_file);
    load_HM80_table(&e->HM80_rock, HM80_rock_table_file);

    // ANEOS
    set_ANEOS_iron(&e->ANEOS_iron);
    set_MANEOS_forsterite(&e->MANEOS_forsterite);

    // SESAME
    set_SESAME_iron(&e->SESAME_iron);

    // Convert to internal units
    // Tillotson
    convert_units_Til(&e->Til_iron, us);
    convert_units_Til(&e->Til_granite, us);
    convert_units_Til(&e->Til_water, us);

    // Hubbard & MacFarlane (1980)
    convert_units_HM80(&e->HM80_HHe, us);
    convert_units_HM80(&e->HM80_ice, us);
    convert_units_HM80(&e->HM80_rock, us);

    // ANEOS
    convert_units_ANEOS(&e->ANEOS_iron, us);
    convert_units_ANEOS(&e->MANEOS_forsterite, us);

    // SESAME
    convert_units_SESAME(&e->SESAME_iron, us);
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {

  message("Equation of state: Planetary.");
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

  io_write_attribute_s(h_grpsph, "Equation of state", "Planetary");
}
#endif

#endif /* SWIFT_PLANETARY_EQUATION_OF_STATE_H */
