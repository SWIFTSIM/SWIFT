/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_EAGLE_THERMAL_H
#define SWIFT_FEEDBACK_STRUCT_EAGLE_THERMAL_H

#include "chemistry_struct.h"
#include "rays_struct.h"

/*! The total number of rays used in stellar feedback */
#define eagle_SNII_feedback_num_of_rays FEEDBACK_NR_RAYS_SNII

/**
 * @brief Feedback fields carried by each hydro particles
 */
struct feedback_part_data {};

/**
 * @brief Extra feedback fields carried by each hydro particles
 */
struct feedback_xpart_data {};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  union {

    /**
     * @brief Values collected from the gas neighbours.
     */
    struct {

      /*! Inverse of normalisation factor used for the enrichment */
      float enrichment_weight_inv;

      /*! Total mass (unweighted) of neighbouring gas particles */
      float ngb_mass;

      /*! Integer number of neighbouring gas particles */
      int num_ngbs;

      /*! SPH-weighted density of the neighbouring gas particles (internal
       * comoving units) */
      float ngb_rho;

      /*! SPH-weighted metallicity of the neighbouring gas particles
       * (dimensionless) */
      float ngb_Z;

      /*! Total (unweighted) number gas neighbours in the stellar kernel */
      int ngb_N;

    } to_collect;

    /**
     * @brief Values to be distributed to the gas neighbours.
     *
     * WARNING: The first two elements must be the enrichment_weight and mass!!
     */
    struct {

      /*! Normalisation factor used for the enrichment */
      float enrichment_weight;

      /*! Mass released */
      float mass;

      /*! Total metal mass released */
      float total_metal_mass;

      /*! Total mass released by each element */
      float metal_mass[chemistry_element_count];

      /*! Total mass released due to SNIa */
      float mass_from_SNIa;

      /*! Total metal mass released due to SNIa */
      float metal_mass_from_SNIa;

      /*! Total iron mass released due to SNIa */
      float Fe_mass_from_SNIa;

      /*! Total mass released due to SNII */
      float mass_from_SNII;

      /*! Total metal mass released due to SNII */
      float metal_mass_from_SNII;

      /*! Total mass released due to AGB */
      float mass_from_AGB;

      /*! Total metal mass released due to AGB */
      float metal_mass_from_AGB;

      /*! Energy change due to thermal and kinetic energy of ejecta */
      float energy;

      /*! Number of SNII energy injections in thermal form */
      int SNII_num_of_thermal_energy_inj;

      /*! Change in energy from SNII feedback energy injection */
      float SNII_delta_u;

    } to_distribute;
  };

  /* Instantiate ray structs for SNII isotropic feedback  */
  struct ray_data SNII_rays[eagle_SNII_feedback_num_of_rays];
};

#endif /* SWIFT_FEEDBACK_STRUCT_EAGLE_THERMAL_H */
