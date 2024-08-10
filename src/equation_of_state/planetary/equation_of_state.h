/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
 * own header file: `equation_of_state/planetary/<eos_type>.h`.
 * See `eos_planetary_material_id` for the available choices.
 *
 * Not all functions are implemented for all EOS types, so not all can be used
 * with all hydro formulations yet.
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "eos_setup.h"
#include "eos_utilities.h"
#include "inline.h"
#include "material_properties.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

/**
 * @brief Returns the internal energy given density and entropy
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float
gas_internal_energy_from_entropy(float density, float entropy,
                                 enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_internal_energy_from_entropy(density, entropy,
                                              &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_internal_energy_from_entropy(density, entropy,
                                              &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_internal_energy_from_entropy(density, entropy,
                                              &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_internal_energy_from_entropy(density, entropy,
                                               &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_internal_energy_from_entropy(density, entropy,
                                                 &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_internal_energy_from_entropy(density, entropy,
                                                 &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_internal_energy_from_entropy(density, entropy,
                                                 &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_internal_energy_from_entropy(density, entropy,
                                                 &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the pressure given density and entropy
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_pressure_from_entropy(
    float density, float entropy, enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_pressure_from_entropy(density, entropy, &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_pressure_from_entropy(density, entropy, &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_pressure_from_entropy(density, entropy,
                                       &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_pressure_from_entropy(density, entropy,
                                        &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_pressure_from_entropy(density, entropy,
                                          &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_pressure_from_entropy(density, entropy,
                                          &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_pressure_from_entropy(density, entropy,
                                          &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_pressure_from_entropy(density, entropy,
                                          &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the entropy given density and pressure.
 *
 * @param density The density \f$\rho\f$.
 * @param pressure The pressure \f$P\f$.
 * @return The entropy \f$A\f$.
 */
__attribute__((always_inline)) INLINE static float gas_entropy_from_pressure(
    float density, float P, enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_entropy_from_pressure(density, P, &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_entropy_from_pressure(density, P, &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_entropy_from_pressure(density, P,
                                       &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_entropy_from_pressure(density, P, &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_entropy_from_pressure(density, P, &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_entropy_from_pressure(density, P, &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_entropy_from_pressure(density, P, &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_entropy_from_pressure(density, P, &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the sound speed given density and entropy
 *
 * @param density The density \f$\rho\f$.
 * @param entropy The entropy \f$S\f$.
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_entropy(
    float density, float entropy, enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_soundspeed_from_entropy(density, entropy,
                                         &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_soundspeed_from_entropy(density, entropy,
                                         &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_soundspeed_from_entropy(density, entropy,
                                         &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_soundspeed_from_entropy(density, entropy,
                                          &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_soundspeed_from_entropy(density, entropy,
                                            &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_soundspeed_from_entropy(density, entropy,
                                            &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_soundspeed_from_entropy(density, entropy,
                                            &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_soundspeed_from_entropy(density, entropy,
                                            &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the entropy given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_entropy_from_internal_energy(float density, float u,
                                 enum eos_planetary_material_id mat_id) {
  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_entropy_from_internal_energy(density, u,
                                              &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_entropy_from_internal_energy(density, u,
                                              &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_entropy_from_internal_energy(density, u,
                                              &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_entropy_from_internal_energy(density, u,
                                               &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_entropy_from_internal_energy(density, u,
                                                 &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_entropy_from_internal_energy(density, u,
                                                 &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_entropy_from_internal_energy(density, u,
                                                 &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_entropy_from_internal_energy(density, u,
                                                 &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the pressure given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_pressure_from_internal_energy(float density, float u,
                                  enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

#ifdef SWIFT_DEBUG_CHECKS
  // Check the corresponding EoS has been initialised for this material ID

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      if (eos.all_idg[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for ideal gas mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /* Tillotson EoS */
    case eos_type_Til:
      if (eos.all_Til[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for Tillotson mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      if (eos.all_Til_custom[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for custom Tillotson mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      if (eos.all_HM80[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for HM80 mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /* SESAME EoS */
    case eos_type_SESAME:
      if (eos.all_SESAME[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for SESAME mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      if (eos.all_ANEOS[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for ANEOS mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      if (eos.all_linear[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for linear mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      if (eos.all_custom[unit_id].mat_id != mat_id)
        error(
            "EoS not enabled for custom mat_id = %d, "
            "please set the corresponding EoS:planetary_use_: 1",
            mat_id);
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      break;
  }
#endif

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_pressure_from_internal_energy(density, u,
                                               &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_pressure_from_internal_energy(density, u,
                                               &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_pressure_from_internal_energy(density, u,
                                               &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_pressure_from_internal_energy(density, u,
                                                &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_pressure_from_internal_energy(density, u,
                                                  &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_pressure_from_internal_energy(density, u,
                                                  &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_pressure_from_internal_energy(density, u,
                                                  &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_pressure_from_internal_energy(density, u,
                                                  &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
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
gas_internal_energy_from_pressure(float density, float P,
                                  enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_internal_energy_from_pressure(density, P,
                                               &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_internal_energy_from_pressure(density, P,
                                               &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_internal_energy_from_pressure(density, P,
                                               &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_internal_energy_from_pressure(density, P,
                                                &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_internal_energy_from_pressure(density, P,
                                                  &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_internal_energy_from_pressure(density, P,
                                                  &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_internal_energy_from_pressure(density, P,
                                                  &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_internal_energy_from_pressure(density, P,
                                                  &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the sound speed given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_soundspeed_from_internal_energy(float density, float u,
                                    enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_soundspeed_from_internal_energy(density, u,
                                                 &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_soundspeed_from_internal_energy(density, u,
                                                 &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_soundspeed_from_internal_energy(density, u,
                                                 &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_soundspeed_from_internal_energy(density, u,
                                                  &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_soundspeed_from_internal_energy(density, u,
                                                    &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_soundspeed_from_internal_energy(density, u,
                                                    &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_soundspeed_from_internal_energy(density, u,
                                                    &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_soundspeed_from_internal_energy(density, u,
                                                    &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the sound speed given density and pressure
 *
 * @param density The density \f$\rho\f$
 * @param P The pressure \f$P\f$
 */
__attribute__((always_inline)) INLINE static float gas_soundspeed_from_pressure(
    float density, float P, enum eos_planetary_material_id mat_id) {

  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_soundspeed_from_pressure(density, P, &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_soundspeed_from_pressure(density, P, &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_soundspeed_from_pressure(density, P,
                                          &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_soundspeed_from_pressure(density, P, &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_soundspeed_from_pressure(density, P,
                                             &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_soundspeed_from_pressure(density, P,
                                             &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_soundspeed_from_pressure(density, P,
                                             &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_soundspeed_from_pressure(density, P,
                                             &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the temperature given density and internal energy
 *
 * @param density The density \f$\rho\f$
 * @param u The internal energy \f$u\f$
 */
__attribute__((always_inline)) INLINE static float
gas_temperature_from_internal_energy(float density, float u,
                                     enum eos_planetary_material_id mat_id) {
  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_temperature_from_internal_energy(density, u,
                                                  &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_temperature_from_internal_energy(density, u,
                                                  &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_temperature_from_internal_energy(density, u,
                                                  &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_temperature_from_internal_energy(density, u,
                                                   &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_temperature_from_internal_energy(density, u,
                                                     &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_temperature_from_internal_energy(density, u,
                                                     &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_temperature_from_internal_energy(density, u,
                                                     &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_temperature_from_internal_energy(density, u,
                                                     &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the density given pressure and temperature
 *
 * @param P The pressure \f$P\f$
 * @param T The temperature \f$T\f$
 */
__attribute__((always_inline)) INLINE static float
gas_density_from_pressure_and_temperature(
    float P, float T, enum eos_planetary_material_id mat_id) {
  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_density_from_pressure_and_temperature(P, T,
                                                       &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_density_from_pressure_and_temperature(P, T,
                                                       &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_density_from_pressure_and_temperature(
          P, T, &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_density_from_pressure_and_temperature(P, T,
                                                        &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_density_from_pressure_and_temperature(
          P, T, &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_density_from_pressure_and_temperature(
          P, T, &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_density_from_pressure_and_temperature(
          P, T, &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_density_from_pressure_and_temperature(
          P, T, &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Returns the density given pressure and internal energy
 *
 * @param P The pressure \f$P\f$
 * @param T The temperature \f$T\f$
 */
__attribute__((always_inline)) INLINE static float
gas_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    enum eos_planetary_material_id mat_id) {
  const enum eos_planetary_type_id type =
      (enum eos_planetary_type_id)(mat_id / eos_type_factor);
  const int unit_id = mat_id % eos_type_factor;

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_type_idg:
      return idg_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_idg[unit_id]);

    /* Tillotson EoS */
    case eos_type_Til:
      return Til_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_Til[unit_id]);

    /* Custom user-provided Tillotson EoS */
    case eos_type_Til_custom:
      return Til_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_Til_custom[unit_id]);

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_type_HM80:
      return HM80_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_HM80[unit_id]);

    /* SESAME EoS */
    case eos_type_SESAME:
      return SESAME_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_SESAME[unit_id]);

    /* ANEOS -- using SESAME-style tables */
    case eos_type_ANEOS:
      return SESAME_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_ANEOS[unit_id]);

    /*! Linear EoS -- user-provided parameters */
    case eos_type_linear:
      return linear_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_linear[unit_id]);

    /*! Generic user-provided custom tables */
    case eos_type_custom:
      return SESAME_density_from_pressure_and_internal_energy(
          P, u, rho_ref, rho_sph, &eos.all_custom[unit_id]);

    default:
      return -1.f;
  }
}

/**
 * @brief Initialize the eos parameters
 *
 * @param e The #eos_parameters
 * @param params The parsed parameters
 */
__attribute__((always_inline)) INLINE static void eos_init(
    struct eos_parameters *e, const struct phys_const *phys_const,
    const struct unit_system *us, struct swift_params *params) {

  // Prepare any/all requested EoS: Set the parameters and material IDs, load
  // tables etc., and convert to internal units

  char eos_file[PARSER_MAX_LINE_SIZE], mat_params_file[PARSER_MAX_LINE_SIZE];
  char param_name[PARSER_MAX_LINE_SIZE];

  // Ideal gas
  if (parser_get_opt_param_int(params, "EoS:planetary_use_idg_def", 0)) {
    set_idg_def(&e->all_idg[eos_unit_id_idg_def], eos_mat_id_idg_def);

    sprintf(param_name, "EoS:idg_def_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_idg_def, mat_params_file,
                        us);
  }

  // Tillotson
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_iron", 0)) {
    set_Til_iron(&e->all_Til[eos_unit_id_Til_iron], eos_mat_id_Til_iron);
    set_Til_u_cold(&e->all_Til[eos_unit_id_Til_iron], eos_mat_id_Til_iron);
    convert_units_Til(&e->all_Til[eos_unit_id_Til_iron], us);

    sprintf(param_name, "EoS:Til_iron_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_Til_iron, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_granite", 0)) {
    set_Til_granite(&e->all_Til[eos_unit_id_Til_granite],
                    eos_mat_id_Til_granite);
    set_Til_u_cold(&e->all_Til[eos_unit_id_Til_granite],
                   eos_mat_id_Til_granite);
    convert_units_Til(&e->all_Til[eos_unit_id_Til_granite], us);

    sprintf(param_name, "EoS:Til_granite_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_Til_granite,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_water", 0)) {
    set_Til_water(&e->all_Til[eos_unit_id_Til_water], eos_mat_id_Til_water);
    set_Til_u_cold(&e->all_Til[eos_unit_id_Til_water], eos_mat_id_Til_water);
    convert_units_Til(&e->all_Til[eos_unit_id_Til_water], us);

    sprintf(param_name, "EoS:Til_water_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_Til_water,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_basalt", 0)) {
    set_Til_basalt(&e->all_Til[eos_unit_id_Til_basalt], eos_mat_id_Til_basalt);
    set_Til_u_cold(&e->all_Til[eos_unit_id_Til_basalt], eos_mat_id_Til_basalt);
    convert_units_Til(&e->all_Til[eos_unit_id_Til_basalt], us);

    sprintf(param_name, "EoS:Til_basalt_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_Til_basalt,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_ice", 0)) {
    set_Til_ice(&e->all_Til[eos_unit_id_Til_ice], eos_mat_id_Til_ice);
    set_Til_u_cold(&e->all_Til[eos_unit_id_Til_ice], eos_mat_id_Til_ice);
    convert_units_Til(&e->all_Til[eos_unit_id_Til_ice], us);

    sprintf(param_name, "EoS:Til_ice_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_Til_ice, mat_params_file,
                        us);
  }

  // Custom user-provided Tillotson
  for (int i_custom = 0; i_custom <= 9; i_custom++) {
    sprintf(param_name, "EoS:planetary_use_Til_custom_%d", i_custom);
    if (parser_get_opt_param_int(params, param_name, 0)) {
      int mat_id = eos_type_Til_custom * eos_type_factor + i_custom;

      sprintf(param_name, "EoS:planetary_Til_custom_%d_param_file", i_custom);
      parser_get_param_string(params, param_name, eos_file);

      set_Til_custom(&e->all_Til_custom[i_custom],
                     (enum eos_planetary_material_id)mat_id, eos_file);
      set_Til_u_cold(&e->all_Til_custom[i_custom],
                     (enum eos_planetary_material_id)mat_id);
      convert_units_Til(&e->all_Til_custom[i_custom], us);

      sprintf(param_name, "EoS:Til_custom_%d_mat_params_file", i_custom);
      parser_get_opt_param_string(params, param_name, mat_params_file,
                                  "NoFile");
      set_material_params(e->all_mat_params, &e->method_params, mat_id, mat_params_file, us);
    }
  }

  // Hubbard & MacFarlane (1980)
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_HHe", 0)) {
    set_HM80_HHe(&e->all_HM80[eos_unit_id_HM80_HHe], eos_mat_id_HM80_HHe);
    parser_get_param_string(params, "EoS:planetary_HM80_HHe_table_file",
                            eos_file);
    load_table_HM80(&e->all_HM80[eos_unit_id_HM80_HHe], eos_file);
    prepare_table_HM80(&e->all_HM80[eos_unit_id_HM80_HHe]);
    convert_units_HM80(&e->all_HM80[eos_unit_id_HM80_HHe], us);

    sprintf(param_name, "EoS:HM80_HHe_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_HM80_HHe, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_ice", 0)) {
    set_HM80_ice(&e->all_HM80[eos_unit_id_HM80_ice], eos_mat_id_HM80_ice);
    parser_get_param_string(params, "EoS:planetary_HM80_ice_table_file",
                            eos_file);
    load_table_HM80(&e->all_HM80[eos_unit_id_HM80_ice], eos_file);
    prepare_table_HM80(&e->all_HM80[eos_unit_id_HM80_ice]);
    convert_units_HM80(&e->all_HM80[eos_unit_id_HM80_ice], us);

    sprintf(param_name, "EoS:HM80_ice_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_HM80_ice, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_rock", 0)) {
    set_HM80_rock(&e->all_HM80[eos_unit_id_HM80_rock], eos_mat_id_HM80_rock);
    parser_get_param_string(params, "EoS:planetary_HM80_rock_table_file",
                            eos_file);
    load_table_HM80(&e->all_HM80[eos_unit_id_HM80_rock], eos_file);
    prepare_table_HM80(&e->all_HM80[eos_unit_id_HM80_rock]);
    convert_units_HM80(&e->all_HM80[eos_unit_id_HM80_rock], us);

    sprintf(param_name, "EoS:HM80_rock_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_HM80_rock,
                        mat_params_file, us);
  }

  // SESAME
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_iron", 0)) {
    set_SESAME_iron(&e->all_SESAME[eos_unit_id_SESAME_iron],
                    eos_mat_id_SESAME_iron);
    parser_get_param_string(params, "EoS:planetary_SESAME_iron_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_iron], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_iron]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_SESAME_iron], us);

    sprintf(param_name, "EoS:SESAME_iron_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_SESAME_iron,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_basalt", 0)) {
    set_SESAME_basalt(&e->all_SESAME[eos_unit_id_SESAME_basalt],
                      eos_mat_id_SESAME_basalt);
    parser_get_param_string(params, "EoS:planetary_SESAME_basalt_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_basalt], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_basalt]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_SESAME_basalt], us);

    sprintf(param_name, "EoS:SESAME_basalt_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_SESAME_basalt,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_water", 0)) {
    set_SESAME_water(&e->all_SESAME[eos_unit_id_SESAME_water],
                     eos_mat_id_SESAME_water);
    parser_get_param_string(params, "EoS:planetary_SESAME_water_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_water], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_SESAME_water]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_SESAME_water], us);

    sprintf(param_name, "EoS:SESAME_water_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_SESAME_water,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SS08_water", 0)) {
    set_SS08_water(&e->all_SESAME[eos_unit_id_SS08_water],
                   eos_mat_id_SS08_water);
    parser_get_param_string(params, "EoS:planetary_SS08_water_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_SS08_water], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_SS08_water]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_SS08_water], us);

    sprintf(param_name, "EoS:SS08_water_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_SS08_water,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_AQUA", 0)) {
    set_AQUA(&e->all_SESAME[eos_unit_id_AQUA], eos_mat_id_AQUA);
    parser_get_param_string(params, "EoS:planetary_AQUA_table_file", eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_AQUA], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_AQUA]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_AQUA], us);

    sprintf(param_name, "EoS:AQUA_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_AQUA, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_CMS19_H", 0)) {
    set_CMS19_H(&e->all_SESAME[eos_unit_id_CMS19_H], eos_mat_id_CMS19_H);
    parser_get_param_string(params, "EoS:planetary_CMS19_H_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_CMS19_H], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_CMS19_H]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_CMS19_H], us);

    sprintf(param_name, "EoS:CMS19_H_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_CMS19_H, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_CMS19_He", 0)) {
    set_CMS19_He(&e->all_SESAME[eos_unit_id_CMS19_He], eos_mat_id_CMS19_He);
    parser_get_param_string(params, "EoS:planetary_CMS19_He_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_CMS19_He], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_CMS19_He]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_CMS19_He], us);

    sprintf(param_name, "EoS:CMS19_He_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_CMS19_He, mat_params_file,
                        us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_CD21_HHe", 0)) {
    set_CD21_HHe(&e->all_SESAME[eos_unit_id_CD21_HHe], eos_mat_id_CD21_HHe);
    parser_get_param_string(params, "EoS:planetary_CD21_HHe_table_file",
                            eos_file);
    load_table_SESAME(&e->all_SESAME[eos_unit_id_CD21_HHe], eos_file);
    prepare_table_SESAME(&e->all_SESAME[eos_unit_id_CD21_HHe]);
    convert_units_SESAME(&e->all_SESAME[eos_unit_id_CD21_HHe], us);

    sprintf(param_name, "EoS:CD21_HHe_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_CD21_HHe, mat_params_file,
                        us);
  }

  // ANEOS -- using SESAME-style tables
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_forsterite",
                               0)) {
    set_ANEOS_forsterite(&e->all_ANEOS[eos_unit_id_ANEOS_forsterite],
                         eos_mat_id_ANEOS_forsterite);
    parser_get_param_string(params, "EoS:planetary_ANEOS_forsterite_table_file",
                            eos_file);
    load_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_forsterite], eos_file);
    prepare_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_forsterite]);
    convert_units_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_forsterite], us);

    sprintf(param_name, "EoS:ANEOS_forsterite_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_ANEOS_forsterite,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_iron", 0)) {
    set_ANEOS_iron(&e->all_ANEOS[eos_unit_id_ANEOS_iron],
                   eos_mat_id_ANEOS_iron);
    parser_get_param_string(params, "EoS:planetary_ANEOS_iron_table_file",
                            eos_file);
    load_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_iron], eos_file);
    prepare_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_iron]);
    convert_units_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_iron], us);

    sprintf(param_name, "EoS:ANEOS_iron_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_ANEOS_iron,
                        mat_params_file, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_Fe85Si15", 0)) {
    set_ANEOS_Fe85Si15(&e->all_ANEOS[eos_unit_id_ANEOS_Fe85Si15],
                       eos_mat_id_ANEOS_Fe85Si15);
    parser_get_param_string(params, "EoS:planetary_ANEOS_Fe85Si15_table_file",
                            eos_file);
    load_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_Fe85Si15], eos_file);
    prepare_table_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_Fe85Si15]);
    convert_units_SESAME(&e->all_ANEOS[eos_unit_id_ANEOS_Fe85Si15], us);

    sprintf(param_name, "EoS:ANEOS_Fe85Si15_mat_params_file");
    parser_get_opt_param_string(params, param_name, mat_params_file, "NoFile");
    set_material_params(e->all_mat_params, &e->method_params, eos_mat_id_ANEOS_Fe85Si15,
                        mat_params_file, us);
  }

  // Linear EoS -- user-provided parameters
  for (int i_linear = 0; i_linear <= 9; i_linear++) {
    sprintf(param_name, "EoS:planetary_use_linear_%d", i_linear);
    if (parser_get_opt_param_int(params, param_name, 0)) {
      int mat_id = eos_type_linear * eos_type_factor + i_linear;

      sprintf(param_name, "EoS:planetary_linear_%d_param_file", i_linear);
      parser_get_param_string(params, param_name, eos_file);
      set_linear_params(&e->all_linear[i_linear],
                        (enum eos_planetary_material_id)mat_id, eos_file);

      sprintf(param_name, "EoS:linear_%d_mat_params_file", i_linear);
      parser_get_opt_param_string(params, param_name, mat_params_file,
                                  "NoFile");
      set_material_params(e->all_mat_params, &e->method_params, mat_id, mat_params_file, us);

      convert_units_linear(&e->all_linear[i_linear], us);

      sprintf(param_name, "EoS:linear_%d_mat_params_file", i_linear);
      parser_get_opt_param_string(params, param_name, mat_params_file,
                                  "NoFile");
      set_material_params(e->all_mat_params, &e->method_params, mat_id, mat_params_file, us);
    }
  }

  // Custom generic tables -- using SESAME-style tables
  for (int i_custom = 0; i_custom <= 9; i_custom++) {
    sprintf(param_name, "EoS:planetary_use_custom_%d", i_custom);
    if (parser_get_opt_param_int(params, param_name, 0)) {
      int mat_id = eos_type_custom * eos_type_factor + i_custom;
      set_custom(&e->all_custom[i_custom],
                 (enum eos_planetary_material_id)mat_id);

      sprintf(param_name, "EoS:planetary_custom_%d_table_file", i_custom);
      parser_get_param_string(params, param_name, eos_file);
      load_table_SESAME(&e->all_custom[i_custom], eos_file);
      prepare_table_SESAME(&e->all_custom[i_custom]);
      convert_units_SESAME(&e->all_custom[i_custom], us);

      sprintf(param_name, "EoS:custom_%d_mat_params_file", i_custom);
      parser_get_opt_param_string(params, param_name, mat_params_file,
                                  "NoFile");
      set_material_params(e->all_mat_params, &e->method_params, mat_id, mat_params_file, us);
    }
  }
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
