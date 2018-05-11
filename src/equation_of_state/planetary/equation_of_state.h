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
 * For any/all of the planetary EOS. Each EOS type's functions are set in its
 * own header file: equation_of_state/planetary/<eos_type>.h
 * See material_id for the available choices.
 *
 * Not all functions are implemented for all EOS types, so not all can be used
 * with all hydro formulations yet.
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "inline.h"
#include "units.h"
#include "physical_constants.h"

/* EOS function headers. */
#include "tillotson.h"
#include "hm80.h"
#include "aneos.h"
#include "sesame.h"

extern struct eos_parameters eos;
/* ------------------------------------------------------------------------- */

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
    id_Til_iron     = type_Till*type_factor,
    id_Til_granite  = type_Till*type_factor + 1,
    id_Til_water    = type_Till*type_factor + 2,
    // Hubbard & MacFarlane (1980) Uranus/Neptune
    id_HM80_HHe     = type_HM80*type_factor,        // Hydrogen-helium atmosphere
    id_HM80_ice     = type_HM80*type_factor + 1,    // H20-CH4-NH3 ice mix
    id_HM80_rock    = type_HM80*type_factor + 2,    // SiO2-MgO-FeS-FeO rock mix
    // ANEOS
    id_ANEOS_iron           = type_ANEOS*type_factor,
    id_MANEOS_forsterite    = type_ANEOS*type_factor + 1,
    // SESAME
    id_SESAME_iron  = type_SESAME*type_factor,
};

/**
 * @brief Returns the internal energy given density and entropy
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

  // Material base type
  switch((int)(mat_id/type_factor)) {

    // Tillotson
    case type_Till:;
        // Select the material parameters
        struct Til_params *mat_Til;
        switch(mat_id) {
            case id_Til_iron:
                mat_Til = &eos.Til_iron;
                break;

            case id_Til_granite:
                mat_Til = &eos.Til_granite;
                break;

            case id_Til_water:
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
            case id_HM80_HHe:
                mat_HM80 = &eos.HM80_HHe;
                break;

            case id_HM80_ice:
                mat_HM80 = &eos.HM80_ice;
                break;

            case id_HM80_rock:
                mat_HM80 = &eos.HM80_rock;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        P = HM80_pressure_from_internal_energy(density, u, mat_HM80);

        break;

    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
            case id_ANEOS_iron:
                mat_ANEOS = &eos.ANEOS_iron;
                break;

            case id_MANEOS_forsterite:
                mat_ANEOS = &eos.MANEOS_forsterite;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        P = ANEOS_pressure_from_internal_energy(density, u, mat_ANEOS);

        break;

    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
            case id_SESAME_iron:
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
            case id_Til_iron:
                mat_Til = &eos.Til_iron;
                break;

            case id_Til_granite:
                mat_Til = &eos.Til_granite;
                break;

            case id_Til_water:
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
            case id_HM80_HHe:
                mat_HM80 = &eos.HM80_HHe;
                break;

            case id_HM80_ice:
                mat_HM80 = &eos.HM80_ice;
                break;

            case id_HM80_rock:
                mat_HM80 = &eos.HM80_rock;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        c = HM80_soundspeed_from_internal_energy(density, u, mat_HM80);

        break;

    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
            case id_ANEOS_iron:
                mat_ANEOS = &eos.ANEOS_iron;
                break;

            case id_MANEOS_forsterite:
                mat_ANEOS = &eos.MANEOS_forsterite;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        c = ANEOS_soundspeed_from_internal_energy(density, u, mat_ANEOS);

        break;

    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
            case id_SESAME_iron:
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
            case id_Til_iron:
                mat_Til = &eos.Til_iron;
                break;

            case id_Til_granite:
                mat_Til = &eos.Til_granite;
                break;

            case id_Til_water:
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
            case id_HM80_HHe:
                mat_HM80 = &eos.HM80_HHe;
                break;

            case id_HM80_ice:
                mat_HM80 = &eos.HM80_ice;
                break;

            case id_HM80_rock:
                mat_HM80 = &eos.HM80_rock;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_HM80 = &eos.HM80_HHe; // Ignored, just here to keep the compiler happy
        };

        c = HM80_soundspeed_from_pressure(density, P, mat_HM80);

        break;

    // ANEOS
    case type_ANEOS:;
        struct ANEOS_params *mat_ANEOS;
        // Select the material parameters
        switch(mat_id) {
            case id_ANEOS_iron:
                mat_ANEOS = &eos.ANEOS_iron;
                break;

            case id_MANEOS_forsterite:
                mat_ANEOS = &eos.MANEOS_forsterite;
                break;

            default:
                error("Unknown material ID! mat_id = %d", mat_id);
                mat_ANEOS = &eos.ANEOS_iron; // Ignored, just here to keep the compiler happy
        };

        c = ANEOS_soundspeed_from_pressure(density, P, mat_ANEOS);

        break;

    // SESAME
    case type_SESAME:;
        struct SESAME_params *mat_SESAME;
        // Select the material parameters
        switch(mat_id) {
            case id_SESAME_iron:
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

    // Set the parameters and material IDs, load tables, etc. for each material
    // Tillotson
    set_Til_iron(&e->Til_iron, id_Til_iron);
    set_Til_granite(&e->Til_granite, id_Til_granite);
    set_Til_water(&e->Til_water, id_Til_water);

    // Hubbard & MacFarlane (1980)
    set_HM80_HHe(&e->HM80_HHe, id_HM80_HHe);
    set_HM80_ice(&e->HM80_ice, id_HM80_ice);
    set_HM80_rock(&e->HM80_rock, id_HM80_rock);

    load_HM80_table(&e->HM80_HHe, HM80_HHe_table_file);
    load_HM80_table(&e->HM80_ice, HM80_ice_table_file);
    load_HM80_table(&e->HM80_rock, HM80_rock_table_file);

    // ANEOS
    set_ANEOS_iron(&e->ANEOS_iron, id_ANEOS_iron);
    set_MANEOS_forsterite(&e->MANEOS_forsterite, id_MANEOS_forsterite);

    // SESAME
    set_SESAME_iron(&e->SESAME_iron, id_SESAME_iron);

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
