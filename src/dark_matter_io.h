/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_DARK_MATTER_IO_H
#define SWIFT_DARK_MATTER_IO_H

#include "io_properties.h"
#include "kernel_gravity.h"

/**
 * @brief Returns the current co-moving softening of a particle
 *
 * @param gp The particle of interest
 * @param grav_props The global gravity properties.
 */
__attribute__((always_inline)) INLINE static float get_softening(const struct dmpart* dmp, const struct gravity_props* restrict grav_props) {
    return dmp->gpart->epsilon;
}


INLINE static void convert_dmpart_pos(const struct engine* e,
                                     const struct dmpart* gp, double* ret) {
    
    const struct space* s = e->s;
    if (s->periodic) {
        ret[0] = box_wrap(gp->x[0] - s->pos_dithering[0], 0.0, s->dim[0]);
        ret[1] = box_wrap(gp->x[1] - s->pos_dithering[1], 0.0, s->dim[1]);
        ret[2] = box_wrap(gp->x[2] - s->pos_dithering[2], 0.0, s->dim[2]);
    } else {
        ret[0] = gp->x[0];
        ret[1] = gp->x[1];
        ret[2] = gp->x[2];
    }
}

INLINE static void convert_dmpart_vel(const struct engine* e,
                                     const struct dmpart* gp, float* ret) {
    
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const struct cosmology* cosmo = e->cosmology;
    const integertime_t ti_current = e->ti_current;
    const double time_base = e->time_base;
    
    const integertime_t ti_beg = get_integer_time_begin(ti_current, gp->time_bin);
    const integertime_t ti_end = get_integer_time_end(ti_current, gp->time_bin);
    
    /* Get time-step since the last kick */
    float dt_kick_grav;
    if (with_cosmology) {
        dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
        dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    } else {
        dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    }
    
    /* Extrapolate the velocities to the current time */
    ret[0] = gp->v_full[0] + gp->gpart->a_grav[0] * dt_kick_grav;
    ret[1] = gp->v_full[1] + gp->gpart->a_grav[1] * dt_kick_grav;
    ret[2] = gp->v_full[2] + gp->gpart->a_grav[2] * dt_kick_grav;
    
    /* Conversion from internal units to peculiar velocities */
    ret[0] *= cosmo->a_inv;
    ret[1] *= cosmo->a_inv;
    ret[2] *= cosmo->a_inv;
}

INLINE static void convert_dmpart_velocity_dispersion(const struct engine* e,
                                                      const struct dmpart* dmp, float* ret) {

    const struct cosmology* cosmo = e->cosmology;

    /* Conversion from internal to physical units */
    ret[0] = dmp->velocity_dispersion * cosmo->a_inv;
}


INLINE static void convert_dmpart_soft(const struct engine* e,
                                      const struct dmpart* dmp, float* ret) {
    
    ret[0] = kernel_gravity_softening_plummer_equivalent_inv *
    get_softening(dmp, e->gravity_properties);
}

/**
 * @brief Specifies which g-particle fields to read from a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void darkmatter_read_as_dmparticles(struct dmpart* dmparts,
                                             struct io_props* list,
                                             int* num_fields) {
    
    /* Say how much we want to read */
    *num_fields = 6;
    
    /* List what we want to read */
    list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                  UNIT_CONV_LENGTH, dmparts, x);
    list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                  UNIT_CONV_SPEED, dmparts, v_full);
    list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                  dmparts, mass);
    list[3] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                  UNIT_CONV_NO_UNITS, dmparts, id_or_neg_offset);
    list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                  UNIT_CONV_LENGTH, dmparts, h);
    list[5] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                  UNIT_CONV_DENSITY, dmparts, rho);

}

/**
 * @brief Specifies which g-particle fields to write to a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void darkmatter_write_dmparts(const struct dmpart* dmparts,
                                              struct io_props* list,
                                              int* num_fields) {
    
    /* Say how much we want to write */
    *num_fields = 5;
    
    /* List what we want to write */
    list[0] = io_make_output_field_convert_dmpart(
                                                 "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, dmparts,
                                                 convert_dmpart_pos, "Co-moving position of the dark matter particles");
    
    list[1] = io_make_output_field_convert_dmpart(
                                                 "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, dmparts, convert_dmpart_vel,
                                                 "Peculiar velocities of the dark matter. This is a * dx/dt where x is the "
                                                 "co-moving position of the particles.");
    
    list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                   dmparts, mass, "Masses of the particles");
    
    list[3] = io_make_output_field(
                                   "ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   id_or_neg_offset, "Unique ID of the particles");
    
    list[4] = io_make_output_field_convert_dmpart(
                                                 "Softenings", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, dmparts, convert_dmpart_soft,
                                                 "Co-moving Plummer-equivalent softening lengths of the particles.");

}

/**
 * @brief Specifies which dm particle fields to write to a dataset
 *
 * @param gparts The dm particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static int sidm_write_dmparts(const struct dmpart* dmparts,
                                              struct io_props* list) {
    
    
    /* List what we want to write */
    list[0] = io_make_output_field("SIDM_events", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.number_of_sidm_events, "Total number of DM-DM collisions the particle has had");
    
    list[1] = io_make_output_field("SIDM_search_radius", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, dmparts, h,
                                   "Co-moving smoothing lengths (FWHM of the kernel) of the DM particles");
    
    list[2] = io_make_output_field("Densities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f, dmparts, rho,
                                   "Co-moving mass densities of the particles");
    
    list[3] = io_make_output_field("SIDM_probability", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   sidm_probability, "Average probability of particle scattering with another (internal units)");
    
    list[4] = io_make_output_field("Number_of_neighbours", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   num_neighbours, "Number of neighbours");
    
    list[5] = io_make_output_field("Time_step_size", FLOAT, 1, UNIT_CONV_TIME, 1.f, dmparts,
                                   time_step_size, "Average probability of particle scattering with another (internal units)");
    
    list[6] = io_make_output_field("Max_SIDM_events_per_timestep", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.max_sidm_events_per_timestep, "Maximum number of DM-DM collisions the particle has had in a single timestep");

    list[7] = io_make_output_field("Cross_section", DOUBLE, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.sigma,"DM particle cross section [internal units]");

    list[8] = io_make_output_field_convert_dmpart("Velocity_dispersion", FLOAT, 1, UNIT_CONV_SPEED, 0.f, dmparts, convert_dmpart_velocity_dispersion,"Velocity dispersion at the position of the DM particle [internal units]");

  return 9;
    
}


#endif /* SWIFT_DARK_MATTER_IO_H */
