/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H

#include "interpolation.h"

/* Number of different type of companion.
   If changed, the IO needs to be updated.
 */
#define GEAR_NUMBER_TYPE_OF_COMPANION 2
#define GEAR_LABELS_SIZE 10

/**
 * @brief Model for the initial mass function.
 *
 * Describe a model such as Kroupa 2001:
 *
 * f(m) = coef[i] * pow(m, exp[i])
 */
struct initial_mass_function {

  /*! Mass limits between IMF parts (n_parts + 1 elements). */
  float *mass_limits;

  /*! Mass fraction computed at the interface between two IMF parts (n_parts + 1
   * elements). */
  float *mass_fraction;

  /*! Exponent of each IMF parts (n_parts elements). */
  float *exp;

  /*! Coefficient of each IMF parts (n_parts elements). */
  float *coef;

  /*! Number of parts (segments) in the function. */
  int n_parts;

  /*! Minimal mass contained in mass_limits, copied for more clarity. */
  float mass_min;

  /*! Maximal mass contained in mass_limits, copied for more clarity. */
  float mass_max;

  /*! Total number of stars (per mass unit) in the IMF. */
  float N_tot;

  /*! Probability to generate a star out of the continuous part of the IMF. */
  float sink_Pc;

  /*! Stellar mass of the continous part of the IMF (in solar mass). */
  float stellar_particle_mass_Msun;

  /*! Minimal mass of stars represented by discrete particles (in solar mass).
   */
  float minimal_discrete_mass_Msun;
};

/**
 * @brief Model for the stellar lifetime.
 */
struct lifetime {

  /*! Coefficients for the log10(m)^2 term */
  float quadratic[3];

  /*! Coefficients for the log10(m) term */
  float linear[3];

  /*! Coefficients for the constant term */
  float constant[3];
};

/**
 * @brief Model for SNIa.
 */
struct supernovae_ia {
  /*! Mass of each element ejected by a single supernovae */
  float yields[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! White dwarf's mass */
  float mass_white_dwarf;

  /*! Minimal mass of the progenitor */
  float mass_min_progenitor;

  /*! Maximal mass of the progenitor */
  float mass_max_progenitor;

  /*! coefficient of the initial mass function for progenitor divided by
   * progenitor_exponent */
  float progenitor_coef_exp;

  /*! exponent of the initial mass function for progenitor */
  float progenitor_exponent;

  /*! exponent of the initial mass function for binaries */
  float companion_exponent;

  struct {
    /*! Initial mass function's coeffcients */
    float coef;

    /*! Maximal mass of the companion */
    float mass_max;

    /*! Minimal mass of the companion */
    float mass_min;
  } companion[GEAR_NUMBER_TYPE_OF_COMPANION];

  /*! Energy released per supernovae */
  float energy_per_supernovae;
};

/**
 * @brief Model for SNII.
 */
struct supernovae_ii {

  /*! Yields not integrated */
  struct {
    /*! Mass fraction of metals ejected by a supernovae. */
    struct interpolation_1d yields[GEAR_CHEMISTRY_ELEMENT_COUNT];

    /*! Total mass fraction ejected. */
    struct interpolation_1d ejected_mass_processed;

    /*! Mass fraction ejected and not processed (=> with the star metallicity).
     */
    struct interpolation_1d ejected_mass_non_processed;
  } raw;

  /*! Yields integrated */
  struct {
    /*! Integrated (over the IMF) mass fraction of metals ejected by a
     * supernovae
     */
    struct interpolation_1d yields[GEAR_CHEMISTRY_ELEMENT_COUNT];

    /*! Total mass fraction ejected (integrated over the IMF) */
    struct interpolation_1d ejected_mass_processed;

    /*! Mass fraction ejected and not processed (=> with the star metallicity)
     */
    struct interpolation_1d ejected_mass_non_processed;
  } integrated;

  /*! Minimal mass for a SNII */
  float mass_min;

  /*! Maximal mass for a SNII */
  float mass_max;

  /*! exponent of the IMF */
  float exponent;

  /*! coefficient of the IMF over the exponent */
  float coef_exp;

  /*! Number of element in the interpolation array */
  int interpolation_size;

  /*! Energy released as a function of progenitor mass */
  struct interpolation_1d energy_per_progenitor_mass;
};

/**
 * @brief Model for radiation.
 */
struct radiation {

  /*! Yields not integrated, mass-only ("M" dimensionality) tables. Store
      log10(value in internal units), pychem-floored, not the raw value --
      see radiation.c's radiation_build_tables()/radiation_read_cgs_array()
      doxygen; every radiation_get_*_from_raw() getter exponentiates back. */
  struct {
    /*! Bolometric luminosity emitted. */
    struct interpolation_1d luminosities;

    /*! Ionization emission rate. */
    struct interpolation_1d dot_N_ion;

    /*! Excess-energy emission rate above the 13.6 eV HI threshold,
        dot_N_ion(m) * mean_excess_photon_energy_HI(m); see
        #radiation_get_mean_excess_photon_energy_HI_from_integral. */
    struct interpolation_1d dot_E_excess;

    /*! Same three quantities, mass x metallicity ("M,Z" dimensionality)
        variant. Populated instead of the fields above when #is_2d is set;
        see radiation_read_data(). Ships with zero test coverage until the
        PARSEC/2D validation track (plan Phase 6) exercises it. */
    struct interpolation_2d luminosities_2d;
    struct interpolation_2d dot_N_ion_2d;
    struct interpolation_2d dot_E_excess_2d;
  } raw;

  /*! Yields integrated (IMF-integrated). Mass-only ("M") tables only: a 2D
      (mass x metallicity) integrated table needs per-metallicity IMF
      integration, which no caller needs yet -- see #is_2d and
      radiation_get_luminosities_from_integral() and friends, which error
      if called on a 2D table. */
  struct {

    /*! Bolometric luminosity emitted. */
    struct interpolation_1d luminosities;

    /*! Ionization emission rate. */
    struct interpolation_1d dot_N_ion;

    /*! Excess-energy emission rate above the 13.6 eV HI threshold. */
    struct interpolation_1d dot_E_excess;
  } integrated;

  /*! Is this a mass x metallicity ("M,Z") table rather than a mass-only
      ("M") one? Set from the Data/Radiation group's own "dimensionality"
      attribute in radiation_read_data(); dispatches between the raw 1D
      and 2D fields above in the read functions, and gates the getters
      that do not yet take a metallicity argument. */
  int is_2d;

  /*! Number of element in the interpolation array (mass axis) */
  int interpolation_size;

  /*! Number of element in the interpolation array (metallicity axis,
      2D tables only) */
  int interpolation_size_metallicity;

  /*! Number of active angular (HEALPix) pixels the HII ionization budget is
      split across (1 = spherical/HEALPix disabled). Set from
      GEARFeedback:HII_angular_nside in radiation_init(). */
  int n_HII_pixels;

  /*! Is a radiation table actually loaded? False when neither
      photoionization nor radiation pressure is enabled (see
      stellar_evolution_props_init()): every other field above is then
      zero-pointered by #radiation_zero_pointers, so callers must check this
      before touching them (radiation_get_*_from_raw()/_from_integral()
      would read a zeroed interpolation table otherwise). Set to 1 at the
      end of radiation_read_data(), 0 by #radiation_zero_pointers. */
  int is_active;
};

/**
 * @brief Model for Stellar winds.
 */
struct stellar_wind {

  /*! Yields not integrated */
  struct {

    /*! energy ejected by stellar winds. */
    struct interpolation_2d ejected_energy;

    /*! mass loss from stellar winds. */
    struct interpolation_2d mass_loss;
  } raw;

  /*! Yields integrated */
  struct {
    /*! Integrated (over the IMF) energy ejected by stellar winds. */
    struct interpolation_2d ejected_energy_per_progenitor_mass;

    /*! Integrated (over the IMF) mass loss from stellar winds */
    struct interpolation_2d mass_loss_per_progenitor_mass;
  } integrated;

  /*! Minimal mass for a SW */
  float mass_min;

  /*! Maximal mass for a SW */
  float mass_max;

  /*! Minimal mass for a SW */
  float metallicity_min;

  /*! Maximal mass for a SW */
  float metallicity_max;

  /*! Number of mass element in the interpolation 2d array*/
  int interpolation_size_m;

  /*! Number of metallicity element in the interpolation 2d array*/
  int interpolation_size_z;
};

/**
 * @brief The complete stellar model.
 */
struct stellar_model {

  /*! Name of the different elements */
  char elements_name[GEAR_CHEMISTRY_ELEMENT_COUNT * GEAR_LABELS_SIZE];

  /* Solar mass abundances read from the chemistry table */
  float solar_abundances[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! The initial mass function */
  struct initial_mass_function imf;

  /*! The stellar lifetime */
  struct lifetime lifetime;

  /*! The supernovae type Ia */
  struct supernovae_ia snia;

  /*! The supernovae type II */
  struct supernovae_ii snii;

  /*! The stellar radiation */
  struct radiation rad;

  /*! The stellar wind */
  struct stellar_wind sw;

  /*! Use a discrete yields approach */
  char discrete_yields;

  /* Filename of the yields table */
  char yields_table[FILENAME_BUFFER_SIZE];

  /* Minimal gravity mass after a discrete star has completely exploded.

     This will be the mass of the gpart's friend of the star. The mass of the
     star will be 0 after it losses all its mass.

     The purpose of this is to avoid zero mass for the gravitsy
     computations. We keep the star so that we know it *existed* and we can
     extract its properties at the end of a run. If we remove the star, then
     we do not have any information about its existence.
     However, since the star is dead/inexistent, the gravity mass must be small
     so that it does not drastically alter the dynamics of the systems. */
  float discrete_star_minimal_gravity_mass;
};

#endif  // SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
