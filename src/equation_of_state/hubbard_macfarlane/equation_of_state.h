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
#ifndef SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H
#define SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/hubbard_macfarlane/equation_of_state.h
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

extern struct eos_parameters eos;
/* ------------------------------------------------------------------------- */

// Hubbard & MacFarlane (1980) Uranus/Neptune table parameters
struct HM80_params {
    int mat_id;
    int num_rho, num_u;
    float log_rho_min, log_rho_max, log_rho_step, inv_log_rho_step, log_u_min,
        log_u_max, log_u_step, inv_log_u_step, bulk_mod;
    float **table_P_rho_u;
};
// Table file names
/// to be read in from the parameter file instead once tested...
#define HM80_HHe_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_HHe.txt"
#define HM80_ice_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_ice.txt"
#define HM80_rock_table_file "/gpfs/data/dc-kege1/gihr_data/P_rho_u_roc.txt"

struct eos_parameters {
    struct HM80_params HM80_HHe, HM80_ice, HM80_rock;
};

// Material identifier flags (material_ID = type_ID * type_factor + unit_ID)
#define type_factor 10
enum type_id {
    type_HM80   = 2
};
enum material_id {
    // Hubbard & MacFarlane (1980) Uranus/Neptune
    HM80_HHe    = type_HM80*type_factor,        // Hydrogen-helium atmosphere
    HM80_ice    = type_HM80*type_factor + 1,    // H20-CH4-NH3 ice mix
    HM80_rock   = type_HM80*type_factor + 2     // SiO2-MgO-FeS-FeO rock mix
};

// Parameter values for each material (cgs units)
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

// Convert from cgs to internal units
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

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u, int mat_id) {
    float P;

    // Select the material parameters
    struct HM80_params *mat;
    switch(mat_id) {
      case HM80_HHe:
          mat = &eos.HM80_HHe;
          break;

      case HM80_ice:
          mat = &eos.HM80_ice;
          break;

      case HM80_rock:
          mat = &eos.HM80_rock;
          break;

      default:
          error("Unknown material ID! mat_id = %d", mat_id);
          mat = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
    };

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
    // Select the material parameters
    struct HM80_params *mat;
    switch(mat_id) {
      case HM80_HHe:
          mat = &eos.HM80_HHe;
          break;

      case HM80_ice:
          mat = &eos.HM80_ice;
          break;

      case HM80_rock:
          mat = &eos.HM80_rock;
          break;

      default:
          error("Unknown material ID! mat_id = %d", mat_id);
          mat = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
    };

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

/**
 * @brief Returns the sound speed given density and pressure
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_pressure(float density, float P, int mat_id) {
    // Select the material parameters
    struct HM80_params *mat;
    switch(mat_id) {
      case HM80_HHe:
          mat = &eos.HM80_HHe;
          break;

      case HM80_ice:
          mat = &eos.HM80_ice;
          break;

      case HM80_rock:
          mat = &eos.HM80_rock;
          break;

      default:
          error("Unknown material ID! mat_id = %d", mat_id);
          mat = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
    };

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
    set_HM80_HHe(&e->HM80_HHe);
    set_HM80_ice(&e->HM80_ice);
    set_HM80_rock(&e->HM80_rock);

    load_HM80_table(&e->HM80_HHe, HM80_HHe_table_file);
    load_HM80_table(&e->HM80_ice, HM80_ice_table_file);
    load_HM80_table(&e->HM80_rock, HM80_rock_table_file);

    // Convert from cgs units to internal units
    convert_units_HM80(&e->HM80_HHe, us);
    convert_units_HM80(&e->HM80_ice, us);
    convert_units_HM80(&e->HM80_rock, us);
}

/**
 * @brief Print the equation of state
 *
 * @param e The #eos_parameters
 */
__attribute__((always_inline)) INLINE static void eos_print(
    const struct eos_parameters *e) {

  message("Equation of state: Hubbard & MacFarlane (1980).");
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

  io_write_attribute_s(h_grpsph, "Equation of state", "Hubbard & MacFarlane (1980)");
}
#endif

#endif /* SWIFT_HUBBARD_MACFARLANE_EQUATION_OF_STATE_H */
