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
#ifndef SWIFT_RT_STRUCT_GEAR_H
#define SWIFT_RT_STRUCT_GEAR_H

/**
 * @file src/rt/GEAR/rt_struct.h
 * @brief Main header file for the GEAR M1 Closure radiative transfer struct.
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {

  /* conserved state vector */
  struct {
    float energy;
    float flux[3];
  } conserved[RT_NGROUPS];

  /* density state vector */
  struct {
    float energy;
    float flux[3];
  } density[RT_NGROUPS];

  /* Fluxes in the conservation law sense */
  struct {
    float energy;
    float flux[3];
  } flux[RT_NGROUPS];

  /* gradients of densities */
  /* for the flux[3][3] quantity:
   *    first index: x, y, z coordinate of the flux.
   *    Second index: gradient along x, y, z direction. */
  struct {
    float energy[3];
    float flux[3][3];
  } gradient[RT_NGROUPS];

  /* cell slope limiter quantities */
  /* array of length two: store min among all neighbours
   * at first index, store max among all neighbours at
   * second index */
  /* the Gizmo-style slope limiting doesn't help for RT as is,
   * so we're skipping it for now. */
  /* struct { */
  /*   float energy[2]; */
  /*   float flux[3][2]; */
  /*   [> float maxr; [> just use the hydro one <] <] */
  /* } limiter[RT_NGROUPS]; */

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* debugging data to store during entire run */
  unsigned long long
      debug_radiation_absorbed_tot; /* how much radiation this part received
                                    from stars during total lifetime */

  /* data to store during one time step */
  int debug_iact_stars_inject;    /* how many stars this part interacted with */
  int debug_calls_iact_gradient;  /* calls from gradient interaction loop */
  int debug_calls_iact_transport; /* calls from transport interaction loop */
  /* skip this for GEAR */
  /* int debug_injection_check;   [> called in a self/rt_injection task? <] */
  /* calls from gradient interaction loop in actual function */
  int debug_calls_iact_gradient_interaction;
  /* calls from transport interaction loop in actual function */
  int debug_calls_iact_transport_interaction;

  int debug_injection_done;  /* calls from ghost1 tasks */
  int debug_gradients_done;  /* finalised computing gradients? */
  int debug_transport_done;  /* transport step done? */
  int debug_thermochem_done; /* thermochemistry done? */
#endif
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* Stellar energy emission that will be injected in to gas.
   * Total energy, not density, not rate! */
  /* TODO: keep this also for RT_HYDRO_CONTROLLED_INJECTION and
   * store results with each hydro-star interaction in here */
  float emission_this_step[RT_NGROUPS];

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* data to store during entire run */
  unsigned long long
      debug_radiation_emitted_tot; /* how much radiation this star emitted
                                      during total lifetime */

  /* data to store during one time step */
  int debug_iact_hydro_inject; /* how many hydro particles this particle
                                  interacted with */
  int debug_emission_rate_set; /* stellar photon emisison rate computed? */
  /* skip this for GEAR */
  /* int debug_injection_check; [> called in a self/rt_injection task? <] */

  float debug_injected_energy[RT_NGROUPS];     /* how much energy this star
                                                  particle actually has injected
                                                  into the gas */
  float debug_injected_energy_tot[RT_NGROUPS]; /* how much energy this star
                                              particle actually has injected
                                              into the gas over the entire
                                              run*/
#endif
};

#endif /* SWIFT_RT_STRUCT_GEAR_H */
