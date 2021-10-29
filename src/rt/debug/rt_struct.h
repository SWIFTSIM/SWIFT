/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_DEBUG_H
#define SWIFT_RT_STRUCT_DEBUG_H

/**
 * @file src/rt/debug/rt_struct.h
 * @brief Main header file for the debug radiative transfer struct.
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {

  /* data to store during entire run */

  /*! how much radiation this part received from stars during total lifetime */
  unsigned long long debug_radiation_absorbed_tot;

  /*! how many interactions this part had with stars in injection prep over
   * total lifetime */
  unsigned long long debug_iact_stars_inject_prep_tot;

  /* data to store during one time step */

  /*! how many stars this part interacted with during preparation*/
  /* Note: It's useless to write this in outputs, as it gets reset
   * at the end of every step. */
  int debug_iact_stars_inject_prep;

  /*! how many stars this part interacted with during injection*/
  /* Note: It's useless to write this in outputs, as it gets reset
   * at the end of every step. */
  int debug_iact_stars_inject;

  /*! called in a self/rt_injection task? */
  int debug_injection_check;

  /*! calls from gradient interaction loop in actual function */
  int debug_calls_iact_gradient_interaction;

  /*! calls from transport interaction loop in actual function */
  int debug_calls_iact_transport_interaction;

  /* Task completion flags */

  /*! calls from ghost1 tasks */
  int debug_injection_done;

  /*! finalised computing gradients? */
  int debug_gradients_done;

  /*! transport step done? */
  int debug_transport_done;

  /*! thermochemistry done? */
  int debug_thermochem_done;
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* data to store during entire run */

  /*! how much radiation this star emitted during total lifetime */
  unsigned long long debug_radiation_emitted_tot;

  /*! how many interactions this star had with parts during
   * injection prep over total lifetime */
  unsigned long long debug_iact_hydro_inject_prep_tot;

  /* data to store during one time step */

  /*! how many hydro particles this particle interacted with
   * during injection */
  int debug_iact_hydro_inject;

  /*! how many hydro particles this particle interacted with
   * during injection prep*/
  int debug_iact_hydro_inject_prep;

  /*! stellar photon emisison rate computed? */
  int debug_emission_rate_set;

  /*! called in a self/rt_injection task? */
  int debug_injection_check;
};

#endif /* SWIFT_RT_STRUCT_DEBUG_H */
