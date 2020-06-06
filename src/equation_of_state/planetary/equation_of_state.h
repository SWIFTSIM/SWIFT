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

  /*! Tillotson */
  eos_planetary_type_Til = 1,

  /*! Hubbard & MacFarlane (1980) Uranus/Neptune */
  eos_planetary_type_HM80 = 2,

  /*! SESAME */
  eos_planetary_type_SESAME = 3,

  /*! ANEOS */
  eos_planetary_type_ANEOS = 4,
};

/**
 * @brief Minor type for the planetary equation of state.
 */
enum eos_planetary_material_id {

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
};

/* Individual EOS function headers. */
#include "hm80.h"
#include "sesame.h"
#include "tillotson.h"

/**
 * @brief The parameters of the equation of state.
 */
struct eos_parameters {
  struct Til_params Til_iron, Til_granite, Til_water, Til_basalt;
  struct HM80_params HM80_HHe, HM80_ice, HM80_rock;
  struct SESAME_params SESAME_iron, SESAME_basalt, SESAME_water, SS08_water, ANEOS_forsterite;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_internal_energy_from_entropy(density, entropy,
                                                     &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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

        case eos_planetary_id_SS08_water:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.SS08_water);
          break;
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_pressure_from_entropy(density, entropy,
                                              &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_entropy_from_pressure(density, P, &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_entropy(density, entropy,
                                                &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_entropy_from_internal_energy(density, u,
                                                     &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_Til_iron:
          return Til_pressure_from_internal_energy(density, u, &eos.Til_iron);
          break;

        case eos_planetary_id_Til_granite:
          return Til_pressure_from_internal_energy(density, u,
                                                   &eos.Til_granite);
          break;

        case eos_planetary_id_Til_water:
          return Til_pressure_from_internal_energy(density, u, &eos.Til_water);
          break;

        case eos_planetary_id_Til_basalt:
          return Til_pressure_from_internal_energy(density, u, &eos.Til_basalt);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_HM80_HHe:
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_HHe);
          break;

        case eos_planetary_id_HM80_ice:
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_ice);
          break;

        case eos_planetary_id_HM80_rock:
          return HM80_pressure_from_internal_energy(density, u, &eos.HM80_rock);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_SESAME_iron:
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_iron);
          break;

        case eos_planetary_id_SESAME_basalt:
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_basalt);
          break;

        case eos_planetary_id_SESAME_water:
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SESAME_water);
          break;

        case eos_planetary_id_SS08_water:
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.SS08_water);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_pressure_from_internal_energy(density, u,
                                                      &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_internal_energy_from_pressure(density, P,
                                                      &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_internal_energy(density, u,
                                                        &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

    /* Tillotson EoS */
    case eos_planetary_type_Til:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* Hubbard & MacFarlane (1980) EoS */
    case eos_planetary_type_HM80:

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* SESAME EoS */
    case eos_planetary_type_SESAME:;

      /* Select the material */
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
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    /* ANEOS -- using SESAME-style tables */
    case eos_planetary_type_ANEOS:;

      /* Select the material */
      switch (mat_id) {
        case eos_planetary_id_ANEOS_forsterite:
          return SESAME_soundspeed_from_pressure(density, P, &eos.ANEOS_forsterite);
          break;

        default:
          error("Unknown material ID! mat_id = %d", mat_id);
          return 0.f;
      };
      break;

    default:
      error("Unknown material type! mat_id = %d", mat_id);
      return 0.f;
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

  // Table file names
  char HM80_HHe_table_file[PARSER_MAX_LINE_SIZE];
  char HM80_ice_table_file[PARSER_MAX_LINE_SIZE];
  char HM80_rock_table_file[PARSER_MAX_LINE_SIZE];
  char SESAME_iron_table_file[PARSER_MAX_LINE_SIZE];
  char SESAME_basalt_table_file[PARSER_MAX_LINE_SIZE];
  char SESAME_water_table_file[PARSER_MAX_LINE_SIZE];
  char SS08_water_table_file[PARSER_MAX_LINE_SIZE];
  char ANEOS_forsterite_table_file[PARSER_MAX_LINE_SIZE];

  // Set the parameters and material IDs, load tables, etc. for each material
  // and convert to internal units
  // Tillotson
  if (parser_get_opt_param_int(params, "EoS:planetary_use_Til", 0)) {
    set_Til_iron(&e->Til_iron, eos_planetary_id_Til_iron);
    set_Til_granite(&e->Til_granite, eos_planetary_id_Til_granite);
    set_Til_water(&e->Til_water, eos_planetary_id_Til_water);
    set_Til_basalt(&e->Til_basalt, eos_planetary_id_Til_basalt);

    convert_units_Til(&e->Til_iron, us);
    convert_units_Til(&e->Til_granite, us);
    convert_units_Til(&e->Til_water, us);
    convert_units_Til(&e->Til_basalt, us);
  }

  // Hubbard & MacFarlane (1980)
  if (parser_get_opt_param_int(params, "EoS:planetary_use_HM80", 0)) {
    set_HM80_HHe(&e->HM80_HHe, eos_planetary_id_HM80_HHe);
    set_HM80_ice(&e->HM80_ice, eos_planetary_id_HM80_ice);
    set_HM80_rock(&e->HM80_rock, eos_planetary_id_HM80_rock);

    parser_get_param_string(params, "EoS:planetary_HM80_HHe_table_file",
                            HM80_HHe_table_file);
    parser_get_param_string(params, "EoS:planetary_HM80_ice_table_file",
                            HM80_ice_table_file);
    parser_get_param_string(params, "EoS:planetary_HM80_rock_table_file",
                            HM80_rock_table_file);

    load_table_HM80(&e->HM80_HHe, HM80_HHe_table_file);
    load_table_HM80(&e->HM80_ice, HM80_ice_table_file);
    load_table_HM80(&e->HM80_rock, HM80_rock_table_file);

    prepare_table_HM80(&e->HM80_HHe);
    prepare_table_HM80(&e->HM80_ice);
    prepare_table_HM80(&e->HM80_rock);

    convert_units_HM80(&e->HM80_HHe, us);
    convert_units_HM80(&e->HM80_ice, us);
    convert_units_HM80(&e->HM80_rock, us);
  }

  // SESAME
  if (parser_get_opt_param_int(params, "EoS:planetary_use_SESAME", 0)) {
    set_SESAME_iron(&e->SESAME_iron, eos_planetary_id_SESAME_iron);
    set_SESAME_basalt(&e->SESAME_basalt, eos_planetary_id_SESAME_basalt);
    set_SESAME_water(&e->SESAME_water, eos_planetary_id_SESAME_water);
    set_SS08_water(&e->SESAME_water, eos_planetary_id_SS08_water);

    parser_get_param_string(params, "EoS:planetary_SESAME_iron_table_file",
                            SESAME_iron_table_file);
    parser_get_param_string(params, "EoS:planetary_SESAME_basalt_table_file",
                            SESAME_basalt_table_file);
    parser_get_param_string(params, "EoS:planetary_SESAME_water_table_file",
                            SESAME_water_table_file);
    parser_get_param_string(params, "EoS:planetary_SS08_water_table_file",
                            SS08_water_table_file);

    load_table_SESAME(&e->SESAME_iron, SESAME_iron_table_file);
    load_table_SESAME(&e->SESAME_basalt, SESAME_basalt_table_file);
    load_table_SESAME(&e->SESAME_water, SESAME_water_table_file);
    load_table_SESAME(&e->SS08_water, SS08_water_table_file);

    prepare_table_SESAME(&e->SESAME_iron);
    prepare_table_SESAME(&e->SESAME_basalt);
    prepare_table_SESAME(&e->SESAME_water);
    prepare_table_SESAME(&e->SS08_water);

    convert_units_SESAME(&e->SESAME_iron, us);
    convert_units_SESAME(&e->SESAME_basalt, us);
    convert_units_SESAME(&e->SESAME_water, us);
    convert_units_SESAME(&e->SS08_water, us);
  }

  // ANEOS -- using SESAME-style tables
  if (parser_get_opt_param_int(params, "EoS:planetary_use_ANEOS", 0)) {
    set_ANEOS_forsterite(&e->ANEOS_forsterite, eos_planetary_id_ANEOS_forsterite);

    parser_get_param_string(params, "EoS:planetary_ANEOS_forsterite_table_file",
                            ANEOS_forsterite_table_file);

    load_table_SESAME(&e->ANEOS_forsterite, ANEOS_forsterite_table_file);

    prepare_table_SESAME(&e->ANEOS_forsterite);

    convert_units_SESAME(&e->ANEOS_forsterite, us);
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
