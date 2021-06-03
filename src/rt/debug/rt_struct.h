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
  unsigned long long
      debug_radiation_absorbed_tot; /* how much radiation this part received
                                    from stars during total lifetime */

  /* data to store during one time step */
  int debug_iact_stars_inject;    /* how many stars this partinteracted with */
  int debug_calls_iact_gradient;  /* calls from gradient interaction loop */
  int debug_calls_iact_transport; /* calls from transport interaction loop */
  int debug_injection_check;      /* called in a self/rt_injection task? */
  /* calls from gradient interaction loop in actual function */
  int debug_calls_iact_gradient_interaction;
  /* calls from transport interaction loop in actual function */
  int debug_calls_iact_transport_interaction;

  int debug_injection_done;  /* calls from ghost1 tasks */
  int debug_gradients_done;  /* finalised computing gradients? */
  int debug_transport_done;  /* transport step done? */
  int debug_thermochem_done; /* thermochemistry done? */
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* data to store during entire run */
  unsigned long long
      debug_radiation_emitted_tot; /* how much radiation this star emitted
                                      during total lifetime */

  /* data to store during one time step */
  int debug_iact_hydro_inject; /* how many hydro particles this particle
                                  interacted with */
  int debug_emission_rate_set; /* stellar photon emisison rate computed? */
  int debug_injection_check;   /* called in a self/rt_injection task? */
};

#endif /* SWIFT_RT_STRUCT_DEBUG_H */
