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

/*! Cap on the number of metallicity rows #radiation.longest_ms_lifetime_myr
    can hold, as a fixed-size struct member rather than a separately
    allocated pointer. Defined here, not in radiation.h (radiation.h
    includes this file, not the other way round), so it is in scope where
    #radiation itself is declared. */
#define RADIATION_MAX_METALLICITY_ROWS 64

/*! Number of provenance attributes #radiation.table_source can hold, one
    per key in #radiation_table_source_keys. */
#define RADIATION_TABLE_SOURCE_COUNT 4

/*! Size of one #radiation.table_source entry. Reported only, never parsed:
    a longer value is truncated, never an error. */
#define RADIATION_TABLE_SOURCE_SIZE 128

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
 *
 * #raw holds per-mass ("M" or "M,Z") tables; #integrated holds the
 * corresponding IMF-integrated ("per Msun of stars formed") ones. Every
 * field stores log10(value in internal units); each raw getter
 * exponentiates back. Exception: #main_sequence_lifetime_2d stores
 * log10(Myr) (see its own doxygen). Each raw quantity is a union of its 1D
 * and 2D interpolation table: #is_2d is fixed for a #radiation instance's
 * whole lifetime, so only one dimensionality is ever live. Every read site
 * gates on #is_2d before touching either member, and #radiation_zero_
 * pointers/#radiation_clean must do the same.
 */
struct radiation {

  struct {
    union {
      /*! Bolometric luminosity emitted. */
      struct interpolation_1d luminosities;

      /*! Bolometric luminosity emitted, mass x metallicity variant. */
      struct interpolation_2d luminosities_2d;
    };

    union {
      /*! Ionization emission rate. */
      struct interpolation_1d dot_N_ion;

      /*! Ionization emission rate, mass x metallicity variant. */
      struct interpolation_2d dot_N_ion_2d;
    };

    union {
      /*! Excess-energy emission rate above the 13.6 eV HI threshold,
          dot_N_ion(m) * mean_excess_photon_energy_HI(m). */
      struct interpolation_1d dot_E_excess;

      /*! Excess-energy emission rate, mass x metallicity variant. */
      struct interpolation_2d dot_E_excess_2d;
    };

    union {
      /*! Photospheric effective temperature (pychem's "Teff" dataset), a
          stellar-evolution diagnostic exposed on the snapshot's star
          particles. Raw (per-mass) only: pychem has no "Integrated_Teff"
          dataset, an IMF average of an effective temperature having no
          single-star meaning. Present only when the table carries the
          dataset (#has_teff). */
      struct interpolation_1d teff;

      /*! #teff, mass x metallicity variant. */
      struct interpolation_2d teff_2d;
    };

    union {
      /*! Non-ionizing PE band (6-11.2 eV) emission rate, read directly
          from pychem's "L_PE" dataset (required whenever #with_ISRF is
          on). Same log-log storage as #luminosities. */
      struct interpolation_1d l_pe;

      /*! #l_pe, mass x metallicity variant. */
      struct interpolation_2d l_pe_2d;
    };

    union {
      /*! Lyman-Werner band (11.2-13.6 eV) emission rate, read directly
          from pychem's "L_LW" dataset when present. See #l_pe. */
      struct interpolation_1d l_lw;

      /*! #l_lw, mass x metallicity variant. */
      struct interpolation_2d l_lw_2d;
    };

    union {
      /*! Photon-number-weighted mean Lyman-Werner photon energy of a
          single star, L_LW/Q_LW over 11.2-13.6 eV, from pychem's
          "MeanPhotonEnergyLW" dataset when present
          (#has_mean_photon_energy_lw). Held in cgs erg, NOT internal
          energy units, matching the same mixed-unit convention
          #dot_E_excess/#dot_N_ion's ratio already produces for
          mean_excess_photon_energy_HI. Same log-log storage as
          #luminosities.

          Below pychem's own LW mass floor, where L_LW is forced to zero
          and Q_LW with it, the table holds the 11.2-13.6 eV band
          midpoint (12.4 eV in erg) rather than a measured mean: a
          finite, in-band placeholder, not a stellar-atmosphere result. */
      struct interpolation_1d mean_photon_energy_lw;

      /*! #mean_photon_energy_lw, mass x metallicity variant. */
      struct interpolation_2d mean_photon_energy_lw_2d;
    };

    /*! Main-sequence duration (TAMS age minus ZAMS age), mass x
        metallicity, 2D-only. Stored as log10(Myr), NOT log10(internal
        units) like the fields above: it is only ever compared against a
        star's age in Myr, and GEAR's own Poirier lifetime model
        (lifetime_get_log_lifetime_from_mass()) already works natively in
        log10(Myr). Caps Q_H/DotEExcess to 0 once a star exceeds this
        window; #luminosities_2d is not capped by it. */
    struct interpolation_2d main_sequence_lifetime_2d;

    /*! (Z, age) -> mass inverse of #main_sequence_lifetime_2d: the mass
        whose PARSEC main-sequence lifetime equals a population's age at a
        given metallicity. Population-feedback cap only (see
        stellar_evolution_compute_preSN_feedback_spart()). Age axis
        identity-resampled at the table's native resolution; Z axis keeps
        the table's own native metallicity nodes, like every other 2D
        field. Stored log10(Msun); 2D-only. */
    struct interpolation_2d main_sequence_lifetime_inverse_2d;
  } raw;

  /*! IMF-integrated yields, read directly from pychem's precomputed,
      number-weighted "Integrated_*" datasets; never integrated on the
      SWIFT side. Same union-per-field layout as #raw. */
  struct {

    union {
      /*! Bolometric luminosity emitted. */
      struct interpolation_1d luminosities;

      /*! Bolometric luminosity emitted, mass x metallicity variant. */
      struct interpolation_2d luminosities_2d;
    };

    union {
      /*! Ionization emission rate. */
      struct interpolation_1d dot_N_ion;

      /*! Ionization emission rate, mass x metallicity variant. */
      struct interpolation_2d dot_N_ion_2d;
    };

    union {
      /*! Excess-energy emission rate above the 13.6 eV HI threshold. */
      struct interpolation_1d dot_E_excess;

      /*! Excess-energy emission rate, mass x metallicity variant. */
      struct interpolation_2d dot_E_excess_2d;
    };

    union {
      /*! IMF-integrated PE emission rate per Msun of stars formed, from
          pychem's "Integrated_L_PE" dataset (required whenever
          #with_ISRF is on). Linear (un-logged), unlike #raw's log10
          storage. */
      struct interpolation_1d l_pe;

      /*! #l_pe, mass x metallicity variant. */
      struct interpolation_2d l_pe_2d;
    };

    union {
      /*! IMF-integrated Lyman-Werner emission rate per Msun of stars
          formed, from pychem's "Integrated_L_LW" dataset. See #l_pe. */
      struct interpolation_1d l_lw;

      /*! #l_lw, mass x metallicity variant. */
      struct interpolation_2d l_lw_2d;
    };

    union {
      /*! Photon-number-weighted mean Lyman-Werner photon energy of the
          whole population formed between the IMF's own mass_min and each
          tabulated mass, from pychem's "Integrated_MeanPhotonEnergyLW"
          dataset. In cgs erg, like #raw.mean_photon_energy_lw.

          Unlike every other "Integrated_*" dataset this is INTENSIVE, not
          per Msun of stars formed: it is the ratio
          Integrated_L_LW/Integrated_Q_LW of two cumulative integrals, so
          it must not be rescaled by a star particle's birth mass. It is
          also not a cumulative quantity a consumer may difference across
          a mass window: a difference of two ratios is not the ratio over
          the window, and pychem does not export Integrated_Q_LW, so no
          window mean is recoverable. Only the [mass_min, m] value the
          dataset itself tabulates is meaningful. Stored logged, unlike
          the linear "Integrated_*" fields above, because it is read as a
          value at a single mass rather than as a difference. */
      struct interpolation_1d mean_photon_energy_lw;

      /*! #mean_photon_energy_lw, mass x metallicity variant. */
      struct interpolation_2d mean_photon_energy_lw_2d;
    };
  } integrated;

  /*! Is this a mass x metallicity ("M,Z") table rather than mass-only
      ("M")? Selects the live member of each union above. */
  char is_2d;

  /*! Number of element in the interpolation array (mass axis) */
  int interpolation_size;

  /*! Number of active angular (HEALPix) pixels the HII ionization budget is
      split across (1 = spherical/HEALPix disabled). Set from
      GEARFeedback:HII_angular_nside in radiation_init(). */
  int n_HII_pixels;

  /*! Longest tabulated MS lifetime (Myr) per native metallicity row,
      sharing its index space with #main_sequence_lifetime_inverse_2d's own
      x axis, which holds the same native log10(Z) nodes. Reduced at read
      time from MainSequenceLifetimeInverseExcluded (2D tables only).
      FLT_MAX (not INFINITY: -ffast-math disallows it) for a row with no
      excluded cells. */
  float longest_ms_lifetime_myr[RADIATION_MAX_METALLICITY_ROWS];

  /*! Longest tabulated age across every metallicity row (Myr), from the
      "Age" dataset's "age_max_myr" attribute (2D tables only, 0
      otherwise). A population older than this has no grid point at all in
      #main_sequence_lifetime_inverse_2d, unlike #longest_ms_lifetime_myr
      (which covers ages beyond one row's own tabulated lifetime, still
      in-grid). */
  float age_max_myr;

  /*! Is a radiation table actually loaded? False when neither
      photoionization nor radiation pressure is enabled: every other field
      is then zero-pointered by #radiation_zero_pointers, so callers must
      check this before touching them. Set to 1 at the end of
      radiation_read_data(), 0 by #radiation_zero_pointers. */
  int is_active;

  /*! Is the local Lyman-Werner/PE feedback
      (GEARFeedback:with_interstellar_radiation_field) on? Set in
      radiation_init(), before radiation_read_data() is
      called: the latter requires L_PE/L_LW/Integrated_L_PE/
      Integrated_L_LW in the table whenever this is set (error() otherwise;
      no Teff-based fallback). Persists across restart as a plain scalar
      (#radiation_dump/#radiation_restore), since params is NULL on
      restart. */
  char with_ISRF;

  /*! Does the loaded table carry a "Teff" dataset? File-derived like
      #is_2d: re-probed on every read, including restart, never
      round-tripped. A table generated before pychem exported Teff has
      none, and must still load: every read of #raw.teff is gated on this
      flag. */
  char has_teff;

  /*! The radiation table's own provenance attributes, one entry per key in
      #radiation_table_source_keys and in that order. GEARFeedback:yields_table
      names a file, and the name says nothing about which photon budget the run
      used: a Pop II and a Pop III table can each be staged under any name, and
      the choice moves Q_H. An entry is empty when the table does not carry
      that attribute; all of them are optional. File-derived, re-read on every
      read including restart. */
  char table_source[RADIATION_TABLE_SOURCE_COUNT][RADIATION_TABLE_SOURCE_SIZE];

  /*! Lowest mass the table tabulates (Msun), as read from the file,
      before the resampling onto #interpolation_size points. */
  float table_mass_min;

  /*! Highest mass the table tabulates (Msun). See #table_mass_min. */
  float table_mass_max;

  /*! Number of native mass points the table carries, i.e. the file's own
      "nm" attribute, not the resampled #interpolation_size. */
  int table_n_mass;

  /*! Number of native metallicity rows (2D tables only, 0 otherwise). */
  int table_n_metallicity;

  /*! Lowest tabulated metallicity, in mass fraction (2D tables only, 0
      otherwise). */
  float table_metallicity_min;

  /*! Highest tabulated metallicity, in mass fraction. See
      #table_metallicity_min. */
  float table_metallicity_max;

  /*! Does this table carry pychem's "MeanPhotonEnergyLW" and
      "Integrated_MeanPhotonEnergyLW" datasets? Both are optional: a table
      generated before pychem exported them leaves
      #raw.mean_photon_energy_lw / #integrated.mean_photon_energy_lw
      unbuilt, and the mean photon energy is then simply not reported. It
      is a diagnostic either way: the H2 photodissociation rate reads
      #RADIATION_SIGMA_H2_OVER_E_LW_CGS. */
  char has_mean_photon_energy_lw;
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

  /*! Minimal metallicity for a SW */
  float metallicity_min;

  /*! Maximal metallicity for a SW */
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

  /*! Solar mass abundances read from the chemistry table */
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

  /*! Filename of the yields table */
  char yields_table[FILENAME_BUFFER_SIZE];

  /*! Floor on a discrete star's gravity mass once it has fully exploded:
      the gravity solver cannot handle an exactly-zero mass. The star is
      kept (not removed) so its properties stay available at the end of a
      run; this mass is kept small enough not to perturb the dynamics. */
  float discrete_star_minimal_gravity_mass;
};

#endif  // SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
