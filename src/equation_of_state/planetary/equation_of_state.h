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
#include "inline.h"
#include "physical_constants.h"
#include "restart.h"
#include "units.h"

extern struct eos_parameters eos;

/*! Material identifier flags (material_ID = type_ID * type_factor + unit_ID) */
#define eos_planetary_type_factor 100

/**
 * @brief Master type for the planetary equation of state.
 */
enum eos_planetary_type_id {

  /*! Ideal gas */
  eos_planetary_type_idg = 0,

  /*! Tillotson */
  eos_planetary_type_Til = 1,

  /*! Hubbard & MacFarlane (1980) Uranus/Neptune */
  eos_planetary_type_HM80 = 2,

  /*! SESAME */
  eos_planetary_type_SESAME = 3,

  /*! ANEOS */
  eos_planetary_type_ANEOS = 4,

  /*! Custom */
  eos_planetary_type_custom = 9,
};

/**
 * @brief Minor type for the planetary equation of state.
 */
enum eos_planetary_material_id {

  /* Ideal gas */

  /*! Default
   * (adiabatic index set at configure time by --with-adiabatic-index, default
   *  value is γ=5/3)
   */
  eos_planetary_id_idg_def = eos_planetary_type_idg * eos_planetary_type_factor,

  /* Tillotson */

  /*! Tillotson iron */
  eos_planetary_id_Til_iron =
      eos_planetary_type_Til * eos_planetary_type_factor,

  /*! Tillotson granite */
  eos_planetary_id_Til_granite =
      eos_planetary_type_Til * eos_planetary_type_factor + 1,

  /*! Tillotson water */
  eos_planetary_id_Til_water =
      eos_planetary_type_Til * eos_planetary_type_factor + 2,

  /*! Tillotson basalt */
  eos_planetary_id_Til_basalt =
      eos_planetary_type_Til * eos_planetary_type_factor + 3,

  /* Hubbard & MacFarlane (1980) Uranus/Neptune */

  /*! Hydrogen-helium atmosphere */
  eos_planetary_id_HM80_HHe =
      eos_planetary_type_HM80 * eos_planetary_type_factor,

  /*! H20-CH4-NH3 ice mix */
  eos_planetary_id_HM80_ice =
      eos_planetary_type_HM80 * eos_planetary_type_factor + 1,

  /*! SiO2-MgO-FeS-FeO rock mix */
  eos_planetary_id_HM80_rock =
      eos_planetary_type_HM80 * eos_planetary_type_factor + 2,

  /* SESAME */

  /*! SESAME iron 2140 */
  eos_planetary_id_SESAME_iron =
      eos_planetary_type_SESAME * eos_planetary_type_factor,

  /*! SESAME basalt 7530 */
  eos_planetary_id_SESAME_basalt =
      eos_planetary_type_SESAME * eos_planetary_type_factor + 1,

  /*! SESAME water 7154 */
  eos_planetary_id_SESAME_water =
      eos_planetary_type_SESAME * eos_planetary_type_factor + 2,

  /*! Senft & Stewart (2008) SESAME-like water */
  eos_planetary_id_SS08_water =
      eos_planetary_type_SESAME * eos_planetary_type_factor + 3,

  /* ANEOS */

  /*! ANEOS forsterite (Stewart et al. 2019) -- in SESAME-style tables */
  eos_planetary_id_ANEOS_forsterite =
      eos_planetary_type_ANEOS * eos_planetary_type_factor,

  /*! ANEOS iron (Stewart 2020) -- in SESAME-style tables */
  eos_planetary_id_ANEOS_iron =
      eos_planetary_type_ANEOS * eos_planetary_type_factor + 1,

  /*! ANEOS Fe85Si15 (Stewart 2020) -- in SESAME-style tables */
  eos_planetary_id_ANEOS_Fe85Si15 =
      eos_planetary_type_ANEOS * eos_planetary_type_factor + 2,
};

/* Individual EOS function headers. */
#include "hm80.h"
#include "ideal_gas.h"
#include "sesame.h"
#include "tillotson.h"

/**
 * @brief The parameters of the equation of state.
 */
struct eos_parameters {
  struct idg_params idg_def;
  struct Til_params Til_iron, Til_granite, Til_water, Til_basalt;
  struct HM80_params HM80_HHe, HM80_ice, HM80_rock;
  struct SESAME_params SESAME_iron, SESAME_basalt, SESAME_water, SS08_water;
  struct SESAME_params ANEOS_forsterite, ANEOS_iron, ANEOS_Fe85Si15;
  struct SESAME_params custom[10];
};

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_internal_energy_from_entropy(density, entropy,
                                                  &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_internal_energy_from_entropy(density, entropy,
                                                  &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_internal_energy_from_entropy(density, entropy,
                                                  &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_internal_energy_from_entropy(density, entropy,
                                                  &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_internal_energy_from_entropy(density, entropy,
                                                  &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_internal_energy_from_entropy(density, entropy,
                                                   &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_internal_energy_from_entropy(density, entropy,
                                                   &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_internal_energy_from_entropy(density, entropy,
                                                   &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_internal_energy_from_entropy(density, entropy,
                                                 &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_pressure_from_entropy(density, entropy, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_pressure_from_entropy(density, entropy, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_pressure_from_entropy(density, entropy, &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_pressure_from_entropy(density, entropy, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_pressure_from_entropy(density, entropy, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_pressure_from_entropy(density, entropy, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_pressure_from_entropy(density, entropy, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_pressure_from_entropy(density, entropy, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_pressure_from_entropy(density, entropy,
                                          &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_entropy_from_pressure(density, P, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_entropy_from_pressure(density, P, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_entropy_from_pressure(density, P, &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_entropy_from_pressure(density, P, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_entropy_from_pressure(density, P, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_entropy_from_pressure(density, P, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_entropy_from_pressure(density, P, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_entropy_from_pressure(density, P, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_entropy_from_pressure(density, P, &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_entropy_from_pressure(density, P, &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_entropy_from_pressure(density, P, &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_entropy_from_pressure(density, P, &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_entropy_from_pressure(density, P,
                                              &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_entropy_from_pressure(density, P, &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_entropy_from_pressure(density, P, &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_entropy_from_pressure(density, P, &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_soundspeed_from_entropy(density, entropy, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_soundspeed_from_entropy(density, entropy, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_soundspeed_from_entropy(density, entropy,
                                             &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_soundspeed_from_entropy(density, entropy, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_soundspeed_from_entropy(density, entropy, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_soundspeed_from_entropy(density, entropy, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_soundspeed_from_entropy(density, entropy, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_soundspeed_from_entropy(density, entropy, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_soundspeed_from_entropy(density, entropy,
                                            &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_entropy_from_internal_energy(density, u, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_entropy_from_internal_energy(density, u, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_entropy_from_internal_energy(density, u, &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_entropy_from_internal_energy(density, u, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_entropy_from_internal_energy(density, u, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_entropy_from_internal_energy(density, u, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_entropy_from_internal_energy(density, u, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_entropy_from_internal_energy(density, u, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_entropy_from_internal_energy(density, u,
                                                 &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.idg_def.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_idg_def: 1");
#endif
          return idg_pressure_from_internal_energy(density, u, &eos.idg_def);
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown material ID! mat_id = %d", mat_id);
#endif
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.Til_iron.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_Til_iron: 1");
#endif
          return Til_pressure_from_internal_energy(density, u, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.Til_granite.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_Til_granite: 1");
#endif
          return Til_pressure_from_internal_energy(density, u,
                                                   &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.Til_water.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_Til_water: 1");
#endif
          return Til_pressure_from_internal_energy(density, u, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.Til_basalt.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_Til_basalt: 1");
#endif
          return Til_pressure_from_internal_energy(density, u, &eos.Til_basalt);
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown material ID! mat_id = %d", mat_id);
#endif
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.HM80_HHe.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_HM80_HHe: 1");
#endif
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.HM80_ice.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_HM80_ice: 1");
#endif
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.HM80_rock.mat_id != mat_id)
            error("EoS not enabled. Please set EoS:planetary_use_HM80_rock: 1");
#endif
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_rock);
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown material ID! mat_id = %d", mat_id);
#endif
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.SESAME_iron.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_SESAME_iron: 1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.SESAME_basalt.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_SESAME_basalt: "
                "1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.SESAME_water.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_SESAME_water: "
                "1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.SS08_water.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_SS08_water: 1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SS08_water);
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown material ID! mat_id = %d", mat_id);
#endif
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.ANEOS_forsterite.mat_id != mat_id)
            error(
                "EoS not enabled. Please set "
                "EoS:planetary_use_ANEOS_forsterite: 1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.ANEOS_iron.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_ANEOS_iron: 1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
#ifdef SWIFT_DEBUG_CHECKS
          if (eos.ANEOS_Fe85Si15.mat_id != mat_id)
            error(
                "EoS not enabled. Please set EoS:planetary_use_ANEOS_Fe85Si15: "
                "1");
#endif
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.ANEOS_Fe85Si15);
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown material ID! mat_id = %d", mat_id);
#endif
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
#ifdef SWIFT_DEBUG_CHECKS
      if (eos.custom[i_custom].mat_id != mat_id)
        error("EoS not enabled. Please set EoS:planetary_use_custom_%d: 1",
              i_custom);
#endif
      return SESAME_pressure_from_internal_energy(density, u,
                                                  &eos.custom[i_custom]);
      break;
    }

    default:
#ifdef SWIFT_DEBUG_CHECKS
      error("Unknown material type! mat_id = %d", mat_id);
#endif
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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_internal_energy_from_pressure(density, P, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_internal_energy_from_pressure(density, P, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_internal_energy_from_pressure(density, P,
                                                   &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_internal_energy_from_pressure(density, P, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_internal_energy_from_pressure(density, P, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_internal_energy_from_pressure(density, P, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_internal_energy_from_pressure(density, P, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_internal_energy_from_pressure(density, P, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_internal_energy_from_pressure(density, P,
                                                  &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_soundspeed_from_internal_energy(density, u, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_soundspeed_from_internal_energy(density, u, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_soundspeed_from_internal_energy(density, u,
                                                     &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_soundspeed_from_internal_energy(density, u,
                                                     &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_soundspeed_from_internal_energy(density, u,
                                                     &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_soundspeed_from_internal_energy(density, u,
                                                      &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_soundspeed_from_internal_energy(density, u,
                                                      &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_soundspeed_from_internal_energy(density, u,
                                                      &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_soundspeed_from_internal_energy(density, u,
                                                    &eos.custom[i_custom]);
      break;
    }

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
      (enum eos_planetary_type_id)(mat_id / eos_planetary_type_factor);

  /* Select the material base type */
  switch (type) {

    /* Ideal gas EoS */
    case eos_planetary_type_idg:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_idg_def:
          return idg_soundspeed_from_pressure(density, P, &eos.idg_def);
          break;

        default:
          return -1.f;
      };
      break;

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_soundspeed_from_pressure(density, P, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_soundspeed_from_pressure(density, P, &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_soundspeed_from_pressure(density, P, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_soundspeed_from_pressure(density, P, &eos.Til_basalt);
          break;

        default:
          return -1.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_soundspeed_from_pressure(density, P, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_soundspeed_from_pressure(density, P, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_soundspeed_from_pressure(density, P, &eos.HM80_rock);
          break;

        default:
          return -1.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_soundspeed_from_pressure(density, P, &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_soundspeed_from_pressure(density, P,
                                                 &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_soundspeed_from_pressure(density, P, &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_soundspeed_from_pressure(density, P, &eos.SS08_water);
          break;

        default:
          return -1.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material of this type */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_pressure(density, P,
                                                 &eos.ANEOS_forsterite);
          break;

        case eos_planetary_id_ANEOS_iron:
          return SESAME_soundspeed_from_pressure(density, P, &eos.ANEOS_iron);
          break;

        case eos_planetary_id_ANEOS_Fe85Si15:
          return SESAME_soundspeed_from_pressure(density, P,
                                                 &eos.ANEOS_Fe85Si15);
          break;

        default:
          return -1.f;
      };
      break;

    /*! Generic user-provided custom tables */
    case eos_planetary_type_custom: {
      const int i_custom =
          mat_id - eos_planetary_type_custom * eos_planetary_type_factor;
      return SESAME_soundspeed_from_pressure(density, P, &eos.custom[i_custom]);
      break;
    }

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

  // Ideal gas
  if (parser_get_opt_param_int(params, "EoS:planetary_use_idg_def", 0)) {
    set_idg_def(&e->idg_def, eos_planetary_id_idg_def);
  }

  // Tillotson
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_iron", 0)) {
    set_Til_iron(&e->Til_iron, eos_planetary_id_Til_iron);
    convert_units_Til(&e->Til_iron, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_granite", 0)) {
    set_Til_granite(&e->Til_granite, eos_planetary_id_Til_granite);
    convert_units_Til(&e->Til_granite, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_water", 0)) {
    set_Til_water(&e->Til_water, eos_planetary_id_Til_water);
    convert_units_Til(&e->Til_water, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til_basalt", 0)) {
    set_Til_basalt(&e->Til_basalt, eos_planetary_id_Til_basalt);
    convert_units_Til(&e->Til_basalt, us);
  }

  // Hubbard & MacFarlane (1980)
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_HHe", 0)) {
    char HM80_HHe_table_file[PARSER_MAX_LINE_SIZE];
    set_HM80_HHe(&e->HM80_HHe, eos_planetary_id_HM80_HHe);
    parser_get_param_string(params, "EoS:planetary_HM80_HHe_table_file",
                            HM80_HHe_table_file);
    load_table_HM80(&e->HM80_HHe, HM80_HHe_table_file);
    prepare_table_HM80(&e->HM80_HHe);
    convert_units_HM80(&e->HM80_HHe, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_ice", 0)) {
    char HM80_ice_table_file[PARSER_MAX_LINE_SIZE];
    set_HM80_ice(&e->HM80_ice, eos_planetary_id_HM80_ice);
    parser_get_param_string(params, "EoS:planetary_HM80_ice_table_file",
                            HM80_ice_table_file);
    load_table_HM80(&e->HM80_ice, HM80_ice_table_file);
    prepare_table_HM80(&e->HM80_ice);
    convert_units_HM80(&e->HM80_ice, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80_rock", 0)) {
    char HM80_rock_table_file[PARSER_MAX_LINE_SIZE];
    set_HM80_rock(&e->HM80_rock, eos_planetary_id_HM80_rock);
    parser_get_param_string(params, "EoS:planetary_HM80_rock_table_file",
                            HM80_rock_table_file);
    load_table_HM80(&e->HM80_rock, HM80_rock_table_file);
    prepare_table_HM80(&e->HM80_rock);
    convert_units_HM80(&e->HM80_rock, us);
  }

  // SESAME
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_iron", 0)) {
    char SESAME_iron_table_file[PARSER_MAX_LINE_SIZE];
    set_SESAME_iron(&e->SESAME_iron, eos_planetary_id_SESAME_iron);
    parser_get_param_string(params, "EoS:planetary_SESAME_iron_table_file",
                            SESAME_iron_table_file);
    load_table_SESAME(&e->SESAME_iron, SESAME_iron_table_file);
    prepare_table_SESAME(&e->SESAME_iron);
    convert_units_SESAME(&e->SESAME_iron, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_basalt", 0)) {
    char SESAME_basalt_table_file[PARSER_MAX_LINE_SIZE];
    set_SESAME_basalt(&e->SESAME_basalt, eos_planetary_id_SESAME_basalt);
    parser_get_param_string(params, "EoS:planetary_SESAME_basalt_table_file",
                            SESAME_basalt_table_file);
    load_table_SESAME(&e->SESAME_basalt, SESAME_basalt_table_file);
    prepare_table_SESAME(&e->SESAME_basalt);
    convert_units_SESAME(&e->SESAME_basalt, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME_water", 0)) {
    char SESAME_water_table_file[PARSER_MAX_LINE_SIZE];
    set_SESAME_water(&e->SESAME_water, eos_planetary_id_SESAME_water);
    parser_get_param_string(params, "EoS:planetary_SESAME_water_table_file",
                            SESAME_water_table_file);
    load_table_SESAME(&e->SESAME_water, SESAME_water_table_file);
    prepare_table_SESAME(&e->SESAME_water);
    convert_units_SESAME(&e->SESAME_water, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SS08_water", 0)) {
    char SS08_water_table_file[PARSER_MAX_LINE_SIZE];
    set_SS08_water(&e->SESAME_water, eos_planetary_id_SS08_water);
    parser_get_param_string(params, "EoS:planetary_SS08_water_table_file",
                            SS08_water_table_file);
    load_table_SESAME(&e->SS08_water, SS08_water_table_file);
    prepare_table_SESAME(&e->SS08_water);
    convert_units_SESAME(&e->SS08_water, us);
  }

  // ANEOS -- using SESAME-style tables
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_forsterite",
                               0)) {
    char ANEOS_forsterite_table_file[PARSER_MAX_LINE_SIZE];
    set_ANEOS_forsterite(&e->ANEOS_forsterite,
                         eos_planetary_id_ANEOS_forsterite);
    parser_get_param_string(params, "EoS:planetary_ANEOS_forsterite_table_file",
                            ANEOS_forsterite_table_file);
    load_table_SESAME(&e->ANEOS_forsterite, ANEOS_forsterite_table_file);
    prepare_table_SESAME(&e->ANEOS_forsterite);
    convert_units_SESAME(&e->ANEOS_forsterite, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_iron", 0)) {
    char ANEOS_iron_table_file[PARSER_MAX_LINE_SIZE];
    set_ANEOS_iron(&e->ANEOS_iron, eos_planetary_id_ANEOS_iron);
    parser_get_param_string(params, "EoS:planetary_ANEOS_iron_table_file",
                            ANEOS_iron_table_file);
    load_table_SESAME(&e->ANEOS_iron, ANEOS_iron_table_file);
    prepare_table_SESAME(&e->ANEOS_iron);
    convert_units_SESAME(&e->ANEOS_iron, us);
  }
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS_Fe85Si15", 0)) {
    char ANEOS_Fe85Si15_table_file[PARSER_MAX_LINE_SIZE];
    set_ANEOS_Fe85Si15(&e->ANEOS_Fe85Si15, eos_planetary_id_ANEOS_Fe85Si15);
    parser_get_param_string(params, "EoS:planetary_ANEOS_Fe85Si15_table_file",
                            ANEOS_Fe85Si15_table_file);
    load_table_SESAME(&e->ANEOS_Fe85Si15, ANEOS_Fe85Si15_table_file);
    prepare_table_SESAME(&e->ANEOS_Fe85Si15);
    convert_units_SESAME(&e->ANEOS_Fe85Si15, us);
  }

  // Custom generic tables -- using SESAME-style tables
  for (int i_custom = 0; i_custom <= 9; i_custom++) {
    char param_name[PARSER_MAX_LINE_SIZE];
    sprintf(param_name, "EoS:planetary_use_custom_%d", i_custom);
    if (parser_get_opt_param_int(params, param_name, 0)) {
      char custom_table_file[PARSER_MAX_LINE_SIZE];
      int mat_id =
          eos_planetary_type_custom * eos_planetary_type_factor + i_custom;
      set_custom(&e->custom[i_custom], (enum eos_planetary_material_id)mat_id);

      sprintf(param_name, "EoS:planetary_custom_%d_table_file", i_custom);
      parser_get_param_string(params, param_name, custom_table_file);
      load_table_SESAME(&e->custom[i_custom], custom_table_file);
      prepare_table_SESAME(&e->custom[i_custom]);
      convert_units_SESAME(&e->custom[i_custom], us);
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
