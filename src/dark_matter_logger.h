/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_DARK_MATTER_LOGGER_H
#define SWIFT_DARK_MATTER_LOGGER_H

/* Local includes */
#include "cell.h"
#include "dark_matter_logger_struct.h"
#include "units.h"


/**
 * @brief log a SIDM event
 *
 * @param dmi the dmpart of the SIDM event
 * @param cosmo the cosmology struct
 * @param num_events number of events per time step
 */
INLINE static void dark_matter_log_num_events(struct sidm_history *sidm_history, const int num_events) {
        sidm_history->num_kicks += num_events;
}

/**
 * @brief log a SIDM event
 *
 * @param dmi the dmpart of the SIDM event
 * @param cosmo the cosmology struct
 * @param num_events number of events per time step
 */
INLINE static void dark_matter_log_num_active_parts(struct sidm_history *sidm_history, const int num_parts) {
    sidm_history->n_parts_active += num_parts;
}

/**
 * @brief log a SIDM event
 *
 * @param dmi the dmpart of the SIDM event
 * @param cosmo the cosmology struct
 * @param num_events number of events per time step
 */
INLINE static void dark_matter_log_total_kinetic_energy(
           struct sidm_history *sidm_history,
           const double energy_before, const double energy_after,
           const double p_before[3], const double p_after[3]) {
    sidm_history->K_before += energy_before;
    sidm_history->K_after += energy_after;
    sidm_history->p_before[0] += p_before[0];
    sidm_history->p_before[1] += p_before[1];
    sidm_history->p_before[2] += p_before[2];
    sidm_history->p_after[0] += p_after[0];
    sidm_history->p_after[1] += p_after[1];
    sidm_history->p_after[2] += p_after[2];

}

/**
 * @brief add a star formation history struct to an other star formation history
 * struct
 *
 * @param sf_add the star formation struct which we want to add to the star
 * formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void dark_matter_logger_add(
             struct sidm_history *sh_update,
             const struct sidm_history *sh_add) {
    
    /* Update the SIDM history structure */
    sh_update->num_kicks += sh_add->num_kicks;
    sh_update->K_before += sh_add->K_before;
    sh_update->K_after += sh_add->K_after;
    sh_update->n_parts_active += sh_add->n_parts_active;
    sh_update->p_before[0] += sh_add->p_before[0];
    sh_update->p_before[1] += sh_add->p_before[1];
    sh_update->p_before[2] += sh_add->p_before[2];
    sh_update->p_after[0] += sh_add->p_after[0];
    sh_update->p_after[1] += sh_add->p_after[1];
    sh_update->p_after[2] += sh_add->p_after[2];
}

/**
 * @brief add a star formation history struct to the engine star formation
 * history accumulator struct
 *
 * @param sf_add the star formation accumulator struct which we want to add to
 * the star formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void dark_matter_logger_add_to_accumulator(
        struct sidm_history_accumulator *sh_update,
        const struct sidm_history *sh_add) {
    
    /* Update the SIDM history structure */
    sh_update->num_kicks = sh_add->num_kicks;
    sh_update->K_before = sh_add->K_before;
    sh_update->K_after = sh_add->K_after;
    sh_update->p_before[0] = sh_add->p_before[0];
    sh_update->p_before[1] = sh_add->p_before[1];
    sh_update->p_before[2] = sh_add->p_before[2];
    sh_update->p_after[0] = sh_add->p_after[0];
    sh_update->p_after[1] = sh_add->p_after[1];
    sh_update->p_after[2] = sh_add->p_after[2];
    sh_update->n_parts_active = sh_add->n_parts_active;

}


/**
 * @brief Initialize the SIDM history structure in the #engine
 *
 * @param sh The pointer to the sidm history structure
 */
INLINE static void dark_matter_logger_accumulator_init(
       struct sidm_history_accumulator *sh) {
    /* Initialize the collecting SIDM structure to zero */
    sh->num_kicks = 0.f;
    sh->K_before = 0.f;
    sh->K_after = 0.f;
    sh->n_parts_active = 0.f;
    sh->p_before[0] = 0.f;
    sh->p_before[1] = 0.f;
    sh->p_before[2] = 0.f;
    sh->p_after[0] = 0.f;
    sh->p_after[1] = 0.f;
    sh->p_after[2] = 0.f;
}

/**
 * @brief Initialize the SIDM history structure in the #engine
 *
 * @param sfh The pointer to the SIDM history structure
 */
INLINE static void dark_matter_logger_init(struct sidm_history *sh) {
    /* Initialize the collecting SIDM structure to zero */
    sh->num_kicks = 0.f;
    sh->K_before = 0.f;
    sh->K_after = 0.f;
    sh->n_parts_active = 0.f;
    sh->p_before[0] = 0.f;
    sh->p_before[1] = 0.f;
    sh->p_before[2] = 0.f;
    sh->p_after[0] = 0.f;
    sh->p_after[1] = 0.f;
    sh->p_after[2] = 0.f;

}

/**
 * @brief Write the final status to a file
 *
 * @param fp The file to write to.
 * @param time the simulation time (time since Big Bang) in internal units.
 * @param a the scale factor.
 * @param z the redshift.
 * @param sh the #sidm_history struct.
 * @param step The time-step of the simulation.
 */
INLINE static void dark_matter_write_to_log_file(
           FILE *fp, const double time, const double a, const double z,
           struct sidm_history_accumulator sh, const int step) {
    
    fprintf(fp, "%6d %16e %12.7f %12.7f %7d %7d %17.7f %17.7f %17.7f %17.7f %17.7f %17.7f %17.7f %17.7f\n", step, time, a, z,
            sh.num_kicks, sh.n_parts_active, sh.K_before, sh.K_after,
            sh.p_before[0], sh.p_before[1], sh.p_before[2],
            sh.p_after[0], sh.p_after[1], sh.p_after[2]);
}

/**
 * @brief Initialize the SIDM logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void dark_matter_logger_init_log_file(
               FILE *fp, const struct unit_system *restrict us,
               const struct phys_const *phys_const) {
    
    /* Write some general text to the logger file */
    fprintf(fp, "# Self-interacting DM History Logger file\n");
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
    fprintf(fp, "# (4) Total number of SIDM kicks in the current time-step.\n");
    fprintf(fp, "# (5) Total number of active particles.\n");
    fprintf(fp, "# (6) Total kinetic energy before kicks (int. units).\n");
    fprintf(fp, "# (7) Total kinetic energy after kicks (int. units).\n");
    fprintf(fp, "# (8-9-10) Linear momentum (x, y, z) before kicks (int. units).\n");
    fprintf(fp, "# (11-12-13) Linear momentum (x, y, z) after kicks (int. units).\n");
    fprintf(fp, "#\n");
    fprintf(fp,
            "# (0)         (1)            (2)          (3)          (4)         "
            "(5)         (6)         (7)         (8)         (9)         (10)"
            "         (11)         (12)         (13)\n");
    fprintf(fp,
            "#            Time             a            z         N. kicks     N. active p."
            "K (before)   K (after)  p[x] (before)  p[y] (before)  p[z] (before)"
            "p[x] (after)  p[y] (after)  p[z] (after)\n");
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication for the SIDM events logger
 *
 * @param e the engine we are running
 */
INLINE static void dark_matter_logger_MPI_Reduce(const struct engine *e) {}
#endif


#endif
