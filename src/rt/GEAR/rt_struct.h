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

  /* Radiation state vector. */
  struct {
    float energy_density;
    float flux[3];
  } radiation[RT_NGROUPS];

  /* Fluxes in the conservation law sense */
  struct {
    float energy;
    float flux[3];
  } flux[RT_NGROUPS];

  /* Particle RT time step. */
  float flux_dt;

  /* gradients of the radiation state. */
  /* for the flux[3][3] quantity:
   *    first index: x, y, z coordinate of the flux.
   *    Second index: gradient along x, y, z direction. */
  struct {
    float energy_density[3];
    float flux[3][3];
  } gradient[RT_NGROUPS];

  /* cell slope limiter quantities */
  /* array of length two: store min among all neighbours
   * at first index, store max among all neighbours at
   * second index */
  /* the Gizmo-style slope limiting doesn't help for RT as is,
   * so we're skipping it for now. */
  /* struct { */
  /*   float energy_density[2]; */
  /*   float flux[3][2]; */
  /*   [> float maxr; [> just use the hydro one <] <] */
  /* } limiter[RT_NGROUPS]; */

  /* Data for thermochemistry */
  struct {
    float mass_fraction_HI;         /* mass fraction taken by HI */
    float mass_fraction_HII;        /* mass fraction taken by HII */
    float mass_fraction_HeI;        /* mass fraction taken by HeI */
    float mass_fraction_HeII;       /* mass fraction taken by HeII */
    float mass_fraction_HeIII;      /* mass fraction taken by HeIII */
    float number_density_electrons; /* number density of electrons */
  } tchem;

  /* Keep track of the actual mass fluxes of the gas species */
  struct {
    float HI;    /* mass fraction taken by HI */
    float HII;   /* mass fraction taken by HII */
    float HeI;   /* mass fraction taken by HeI */
    float HeII;  /* mass fraction taken by HeII */
    float HeIII; /* mass fraction taken by HeIII */
  } mass_flux;

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* debugging data to store during entire run */

  /*! how much radiation this part received from stars during total lifetime */
  unsigned long long debug_radiation_absorbed_tot;

  /* data to store during one time step */

  /*! how many stars this part interacted with during injection*/
  /* Note: It's useless to write this in outputs, as it gets reset
   * at the end of every step. */
  int debug_iact_stars_inject;

  /*! calls from gradient interaction loop in actual function */
  int debug_calls_iact_gradient_interaction;

  /*! calls from transport interaction loop in actual function */
  int debug_calls_iact_transport_interaction;

  /* Task completion flags */

  /*! part got kicked? */
  int debug_kicked;

  /*! calls from ghost1 tasks */
  int debug_injection_done;

  /*! finalised computing gradients? */
  int debug_gradients_done;

  /*! transport step done? */
  int debug_transport_done;

  /*! thermochemistry done? */
  int debug_thermochem_done;

  /* Subcycling flags */

  /*! Current subcycle wrt (last) hydro step */
  int debug_nsubcycles;

#endif
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* Stellar energy emission that will be injected in to gas.
   * Total energy, not density, not rate! */
  float emission_this_step[RT_NGROUPS];

  /*! Neighbour weigths in each octant surrounding the star */
  float octant_weights[8];

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* data to store during entire run */

  /*! how much radiation this star emitted during total lifetime */
  unsigned long long debug_radiation_emitted_tot;

  /* data to store during one time step */

  /*! how many hydro particles this particle interacted with
   * during injection */
  int debug_iact_hydro_inject;

  /*! how many hydro particles this particle interacted with
   * during injection prep*/
  int debug_iact_hydro_inject_prep;

  /*! stellar photon emisison rate computed? */
  int debug_emission_rate_set;

  /*! how much energy this star particle actually has injected into the gas */
  float debug_injected_energy[RT_NGROUPS];

  /*! how much energy this star particle actually has injected into the gas over
   * the entire run*/
  float debug_injected_energy_tot[RT_NGROUPS];

  /*! sum up total weights used during injection to compare consistency */
  float debug_psi_sum;
#endif
};

#endif /* SWIFT_RT_STRUCT_GEAR_H */
