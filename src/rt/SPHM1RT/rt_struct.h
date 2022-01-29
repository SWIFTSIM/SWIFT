/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_RT_STRUCT_SPHM1RT_H
#define SWIFT_RT_STRUCT_SPHM1RT_H

/**
 * @file src/rt/SPHM1RT/rt_struct.h
 * @brief Main header file for no radiative transfer struct.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {

  /*! time step of the gas particle */
  float dt;

  /* conserved state vector in comoving units */
  /* but comoving and physical urad and frad are the same in our convention here
   */
  /* urad: radiation energy per mass */
  /* frad: radiation flux per gas density */
  /* (they are conserved in the sense of energy/mass; assuming mass is equal) */
  struct {
    float urad;
    float frad[3];
  } conserved[RT_NGROUPS];

  /* rate of change of the conserved state vector */
  struct {
    float urad;
    float frad[3];
  } dconserved_dt[RT_NGROUPS];

  /* Store viscosity information in a separate struct. */
  struct {

    /*! Particle radiation flux divergence */
    float divf;

    /*! Particle radiation flux divergence from previous step */
    float divf_previous_step;

    /* parameter to control dissipation */
    float alpha;

  } viscosity[RT_NGROUPS];

  /* Store artificial diffusion information in a separate struct. */
  struct {

    /*! gradient of radiation energy density per gas density */
    float graduradc[3];

    /* parameter to control dissipation */
    float alpha;

  } diffusion[RT_NGROUPS];

  /* Store radiation parameter in a separate struct. */
  struct {

    /*! initial mean opacity */
    float chi[RT_NGROUPS];

    /*! reduced speed of light */
    float cred;

  } params;

  /* Store hydro information in a separate struct. */
  struct {

    /*! "Grad h" term */
    float f;

  } force;
};

/* Additional RT data in star particle struct */
struct rt_spart_data {};

#endif /* SWIFT_RT_STRUCT_SPHM1RT_H */
