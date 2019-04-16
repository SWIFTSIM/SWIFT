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
#ifndef SWIFT_FEEDBACK_STRUCT_EAGLE_H
#define SWIFT_FEEDBACK_STRUCT_EAGLE_H

struct feedback_spart_data {

  struct {

    /* Mass of ejecta */
    float mass;

    /* Total metal mass released */
    float total_metal_mass;

    /* Total mass released by element */
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

    /* Number of type Ia SNe per unit mass */
    float num_SNIa;

    /* Number of type II SNe per unit mass */
    float num_SNII;

    /* Number of SNe in timestep  */
    float num_SNe;

    /* Energy change due to thermal and kinetic energy of ejecta */
    float d_energy;

    /* Probability for heating neighbouring gas particles */
    float heating_probability;

  } to_distribute;

  /* Normalisation factor for density weight fraction for feedback (equivalent
   * to metalweight_norm in EAGLE, see eagle_enrich.c:811) */
  float density_weighted_frac_normalisation_inv;

  /* total mass (unweighted) of neighbouring gas particles */
  float ngb_mass;
};

#endif /* SWIFT_FEEDBACK_STRUCT_EAGLE_H */
