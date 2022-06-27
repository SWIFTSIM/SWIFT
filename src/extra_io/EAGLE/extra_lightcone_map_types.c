/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "black_holes.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "lightcone/lightcone_map.h"
#include "part.h"
#include "stars.h"

/* This object's header */
#include "lightcone/lightcone_map_types.h"

/* Required for the xrays */
#include "extra_io.h"
#include "io_properties.h"

/**
 * @brief Determine time since AGN injection at scale factor a_cross
 *
 * Returns -1 if no there has been no AGN injection. Result is in
 * internal time units.
 *
 * @param xp the #xpart for which we're evaluating the time
 * @param c the #cosmology struct
 * @param a_cross expansion factor at which the particle crosses the lightcone
 */
INLINE static double get_time_since_AGN_injection(const struct xpart *xp,
                                                  const struct cosmology *c,
                                                  double a_cross) {

  /* Check for the case where there has been no AGN injection yet */
  const double last_AGN_injection_scale_factor =
      xp->tracers_data.last_AGN_injection_scale_factor;
  if (last_AGN_injection_scale_factor < 0.0) return -1.0;

  /* Check for heating after lightcone crossing - possible if heated on current
   * time step? */
  if (last_AGN_injection_scale_factor > a_cross) return 0.0;

  /* Find time since the last injection in internal units */
  const double last_AGN_injection_time =
      cosmology_get_time_since_big_bang(c, last_AGN_injection_scale_factor);
  const double time_at_crossing = cosmology_get_time_since_big_bang(c, a_cross);
  const double delta_time = time_at_crossing - last_AGN_injection_time;

  return delta_time;
}

INLINE static int exclude_particle(
    const struct lightcone_props *lightcone_props, const struct engine *e,
    const struct part *p, const struct xpart *xp, double a_cross) {

  /* Get AGN heating temperature */
  const double AGN_delta_T = e->black_holes_properties->AGN_delta_T_desired;

  /* Check if we need to exclude this particle due to recent AGN heating */
  if (lightcone_props->xray_maps_recent_AGN_injection_exclusion_time > 0) {
    const double t = get_time_since_AGN_injection(xp, e->cosmology, a_cross);
    if (t >= 0 &&
        t < lightcone_props->xray_maps_recent_AGN_injection_exclusion_time) {

      /* Check if it is within the exclusion temperature range */
      const double temp_min =
          AGN_delta_T * lightcone_props->xray_maps_recent_AGN_min_temp_factor;
      const double temp_max =
          AGN_delta_T * lightcone_props->xray_maps_recent_AGN_max_temp_factor;
      const double part_temp = cooling_get_temperature(
          e->physical_constants, e->hydro_properties, e->internal_units,
          e->cosmology, e->cooling_func, p, xp);
      if (part_temp > temp_min && part_temp < temp_max) return 1;
    }
  }

  /* Not excluding this particle */
  return 0;
}

/**
 * @brief Make a healpix map of projected erosita-low intrinsic photon flux in
 * each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_erosita_low_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_erosita_low_intrinsic_photons);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        (1 + z_cross));  // photon luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/* erosita-low energy flux */

/**
 * @brief Make a healpix map of projected erosita-low intrinsic energy flux in
 * each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_erosita_low_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_erosita_low_intrinsic_energies);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        pow((1 + z_cross), 2));  // energy luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/* erosita_high photon flux */

/**
 * @brief Make a healpix map of projected erosita-high intrinsic photon flux in
 * each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_erosita_high_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_erosita_high_intrinsic_photons);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        (1 + z_cross));  // photon luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/* erosita-high energy flux */

/**
 * @brief Make a healpix map of projected erosita-high intrinsic energy flux in
 * each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_erosita_high_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_erosita_high_intrinsic_energies);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        pow((1 + z_cross), 2));  // energy luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/* ROSAT photon flux */

/**
 * @brief Make a healpix map of projected ROSAT intrinsic photon flux in each
 * pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_rosat_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_ROSAT_intrinsic_photons);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        (1 + z_cross));  // photon luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/* ROSAT energy flux */

/**
 * @brief Make a healpix map of projected ROSAT intrinsic energy flux in each
 * pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_xray_rosat_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      const double z_cross = (1 / a_cross) - 1;
      const double cdist_cross =
          sqrt(pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2));

      const double luminosity = extra_io_get_xray_fluxes(
          p, xp, e, xray_band_types_ROSAT_intrinsic_energies);

      const double flux =
          luminosity / (4 * M_PI * pow(cdist_cross, 2) *
                        pow((1 + z_cross), 2));  // energy luminosity distance

      return flux;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Make a healpix map of the compton y parameter
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_compton_y_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cooling_function_data *cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      double y_compton = cooling_get_ycompton(phys_const, hydro_props, us,
                                              cosmo, cool_func, p, xp);
      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      /* This angular diameter distance is only correct for flat cosmologies */
#ifdef SWIFT_DEBUG_CHECKS
      if (fabs(e->cosmology->Omega_k) > 0.001)
        error("only implemented for flat cosmology");
#endif

      double pixel_size_2 = lightcone_props->pixel_area_steradians;
      double y_for_map =
          y_compton / (pixel_size_2 * angular_diameter_distance_2);

      return y_for_map;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
      break;
  }
}

/**
 * @brief Make a healpix map of the doppler b parameter
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_doppler_b_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cooling_function_data *cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      double n_e = cooling_get_electron_density(phys_const, hydro_props, us,
                                                cosmo, cool_func, p, xp);

      double rho = hydro_get_physical_density(p, cosmo);

      double m = hydro_get_mass(p);

      const double c = phys_const->const_speed_light_c;

      const double sigma_thompson = phys_const->const_thomson_cross_section;

      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;
      double angular_diameter_distance = sqrt(angular_diameter_distance_2);

      /* This angular diameter distance is only correct for flat cosmologies */
#ifdef SWIFT_DEBUG_CHECKS
      if (fabs(e->cosmology->Omega_k) > 0.001)
        error("only implemented for flat cosmology");
#endif

      /* beware: we need comoving radial peculiar velocities, so no a_cross
         factor is required (reminder: p->v = a^2 \dot{x})! */
      double radial_velocity =
          (p->v[0] * x_cross[0] + p->v[1] * x_cross[1] + p->v[2] * x_cross[2]) /
          angular_diameter_distance;

      double pixel_size_2 = lightcone_props->pixel_area_steradians;

      double b_for_map = n_e * m * sigma_thompson * radial_velocity /
                         (pixel_size_2 * angular_diameter_distance_2 * rho * c);

      return b_for_map;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
      break;
  }
}

/**
 * @brief Make a healpix map of the dispersion meassure
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_dispersion_meassure_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cooling_function_data *cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

      /* Exclude recently AGN heated particles */
      if (exclude_particle(lightcone_props, e, p, xp, a_cross)) return 0.0;

      double n_e = cooling_get_electron_density(phys_const, hydro_props, us,
                                                cosmo, cool_func, p, xp);

      double rho = hydro_get_physical_density(p, cosmo);

      double m = hydro_get_mass(p);

      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      /* This angular diameter distance is only correct for flat cosmologies */
#ifdef SWIFT_DEBUG_CHECKS
      if (fabs(e->cosmology->Omega_k) > 0.001)
        error("only implemented for flat cosmology");
#endif

      double pixel_size_2 = lightcone_props->pixel_area_steradians;
      /* an additional a_cross = 1/(1+z) is part of the integrand */
      double dm_for_map = n_e * m * a_cross /
                          (pixel_size_2 * angular_diameter_distance_2 * rho);

      return dm_for_map;
    } break;
    default:
      /* Not gas, nothing to do */
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
      break;
  }
}
