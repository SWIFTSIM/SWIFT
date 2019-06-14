/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_STARFORMATION_LOGGER_H
#define SWIFT_GEAR_STARFORMATION_LOGGER_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cell.h"
#include "hydro.h"
#include "part.h"
#include "star_formation_logger_struct.h"
#include "units.h"


/**
 * @brief Update the stellar quantities in the current cell after creating
 * the new star particle spart sp.
 *
 * @param time_step, the current time step of the simulation
 * @param sp new created star particle
 * @param sf the star_formation_history struct of the current cell
 */
INLINE static void star_formation_logger_log_new_spart(
    const struct spart *sp, struct star_formation_history *sf) {

  /* Add mass of created sparticle to the total stellar mass in this cell */
  sf->stellar_mass += sp->mass;

  /* Increase the number of stars */
  sf->number_of_stars += 1;

  /* No need to deal with the integrated quantities, only the engine's one is updated */
}

/**
 * @brief Initialize the star formation history struct in the case the cell is
 * inactive
 *
 * @param sf the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_logger_log_inactive_cell(
    struct star_formation_history *sf) {

  /* Initialize the stellar mass to zero */
  sf->stellar_mass = 0.f;

  /* initialize number of stars to zero*/
  sf->number_of_stars = 0;

}
/**
 * @brief Initialize the star formation history structure in the #engine
 *
 * @param sfh The pointer to the star formation history structure
 */
INLINE static void star_formation_logger_init(
    struct star_formation_history *sfh) {
  /* Initialize the collecting SFH structure to zero */
  sfh->new_stellar_mass = 0.f;
  sfh->number_new_stars = 0;
}
/**
 * @brief Initialize the star formation history struct in the #engine when
 * starting the simulation.
 *
 * @param sfh the star_formation_history struct we want to initialize
 */
INLINE static void star_formation_logger_first_init(
    struct star_formation_history *sfh) {
  /* Initialize all values to zero */
  sfh->new_stellar_mass = 0.f;
  sfh->number_new_stars = 0;
  sfh->total_number_stars = 0;
  sfh->total_stellar_mass = 0.f;
}

/**
 * @brief add a star formation history struct to an other star formation history
 * struct
 *
 * @param sf_add the star formation struct which we want to add to the star
 * formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void star_formation_logger_add(
    struct star_formation_history *sf_update,
    const struct star_formation_history *sf_add) {

  /* Update the SFH structure */
  sf_update->number_new_stars += sf_add->number_new_stars;
  sf_update->new_stellar_mass += sf_add->new_stellar_mass;
}

/**
 * @brief add a star formation history struct to an other star formation history
 * struct
 *
 * @param sf_add the star formation struct which we want to add to the star
 * formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void star_formation_logger_add_to_engine(
    struct star_formation_history *sf_update,
    const struct star_formation_history *sf_add) {

  /* Remove self contribution */
  sf_update->total_number_stars -= sf_update->number_new_stars;
  sf_update->total_stellar_mass -= sf_update->new_stellar_mass;

  /* Update the SFH structure */
  sf_update->number_new_stars += sf_add->number_new_stars;
  sf_update->new_stellar_mass += sf_add->new_stellar_mass;
  sf_update->total_number_stars += sf_add->number_new_stars;
  sf_update->total_stellar_mass += sf_add->new_stellar_mass;

}

/**
za * @brief Write the final SFH to a file
 *
 * @param fp The file to write to.
 * @param time the simulation time (time since Big Bang) in internal units.
 * @param a the scale factor.
 * @param z the redshift.
 * @param sf the #star_formation_history struct.
 * @param step The time-step of the simulation.
 */
INLINE static void star_formation_logger_write_to_log_file(
    FILE *fp, const double time, const double a, const double z,
    struct star_formation_history sf, const int step) {

  fprintf(fp, "%6d %16e %12.7f %14e %14ld %14e %14ld %14e\n", step, time, a, z,
          sf.total_number_stars, sf.total_stellar_mass, sf.number_new_stars, sf.new_stellar_mass);
}
/**
 * @brief Initialize the SFH logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void star_formation_logger_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

  /* Write some general text to the logger file */
  fprintf(fp, "# Star Formation History Logger file\n");
  fprintf(fp, "######################################################\n");
  fprintf(fp, "# The quantities are all given in internal physical units!\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0) Simulation step\n");
  fprintf(fp,
          "# (1) Time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp, "#     Unit = %e seconds\n", us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e yr or %e Myr\n", 1.f / phys_const->const_year,
          1.f / phys_const->const_year / 1e6);
  fprintf(fp, "# (2) Scale factor (no unit)\n");
  fprintf(fp, "# (3) Redshift     (no unit)\n");
  fprintf(fp,
          "# (4) Total number of stars formed in the simulation (no unit)\n");
  fprintf(fp, "# (5) Total stellar mass formed in the simulation.\n");
  fprintf(fp, "#     Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(fp, "#     Unit = %e solar mass\n",
          1.f / phys_const->const_solar_mass);
  fprintf(fp, "# (6) Number of stars formed in the current time step (no unit).\n");
  fprintf(fp, "# (7) Mass of stars formed in the current time step.\n");
  fprintf(fp, "#     Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(fp, "#     Unit = %e solar mass\n",
          1.f / phys_const->const_solar_mass);
  fprintf(fp,
          "# (0)         (1)            (2)          (3)                   (4) "
          "      (5)        (6)        (7)\n");
}

/**
 * @brief Add the SFR tracer to the total active SFR of this cell
 *
 * @param p the #part
 * @param xp the #xpart
 *
 * Nothing to do here
 *
 * @param sf the SFH logger struct
 * @param dt_star The length of the time-step in physical internal units.
 */
INLINE static void star_formation_logger_log_active_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf, const double dt_star) {}

/**
 * @brief Add the SFR tracer to the total inactive SFR of this cell as long as
 * the SFR tracer is larger than 0
 *
 * Nothing to do here
 *
 * @param p the #part
 * @param xp the #xpart
 * @param sf the SFH logger struct
 */
INLINE static void star_formation_logger_log_inactive_part(
    const struct part *p, const struct xpart *xp,
    struct star_formation_history *sf) {}

#endif /* SWIFT_GEAR_STARFORMATION_LOGGER_H */
