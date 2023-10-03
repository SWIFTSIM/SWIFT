/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_MHD_H
#define SWIFT_NONE_MHD_H

/**
 * @brief Returns the magnetic energy contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0) {

  return 0.f;
}

/**
 * @brief Returns the magnetic helicity contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_helicity(
    const struct part *p, const struct xpart *xp) {

  return 0.f;
}

/**
 * @brief Returns the magnetic cross-helicity contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_cross_helicity(
    const struct part *p, const struct xpart *xp) {

  return 0.f;
}

/**
 * @brief Returns the magnetic field divergence of the particle.
 *
 * This is (div B) / (B / h) and is hence dimensionless.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  return 0.f;
}

/**
 * @brief Computes the MHD time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float mhd_compute_timestep(
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo,
    const float mu_0) {

  return FLT_MAX;
}

/**
 * @brief Compute the signal velocity between two gas particles
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float mhd_signal_velocity(
    const float dx[3], const struct part *pi, const struct part *pj,
    const float mu_ij, const float beta) {

  error("Calling an MHD signal velocity when compiling without MHD!");
  return -1.f;
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density loop over neighbours. Typically, all fields of the
 * density sub-structure of a particle get zeroed in here.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_init_part(
    struct part *p) {}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * Additional quantities such as velocity gradients will also get the final
 * terms added to them here.
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_density(
    struct part *p, const struct cosmology *cosmo) {}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_gradient(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props) {}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after mhd_prepare_gradient.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_reset_gradient(
    struct part *p) {}

/**
 * @brief Finishes the gradient calculation.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void mhd_end_gradient(
    struct part *p) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * In the desperate case where a particle has no neighbours (likely because
 * of the h_max ceiling), set the particle fields to something sensible to avoid
 * NaNs in the next calculations.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_part_has_no_neighbours(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_force(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float dt_alpha, const float mu_0) {}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_reset_acceleration(
    struct part *p) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model
 */
__attribute__((always_inline)) INLINE static void mhd_reset_predicted_values(
    struct part *p, const struct xpart *xp, const struct cosmology *cosmo) {}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Note the different time-step sizes used for the different quantities as they
 * include cosmological factors.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void mhd_predict_extra(
    struct part *p, const struct xpart *xp, const float dt_drift,
    const float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the force and accelerations by the appropiate constants
 * and add the self-contribution term. In most cases, there is little
 * to do here.
 *
 * Cosmological terms are also added/multiplied here.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_force(
    struct part *p, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float mu_0) {}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void mhd_kick_extra(
    struct part *p, struct xpart *xp, const float dt_therm, const float dt_grav,
    const float dt_hydro, const float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {}

/**
 * @brief Converts MHD quantities of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy in the case
 * of hydro for instance.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 */
__attribute__((always_inline)) INLINE static void mhd_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props) {}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_first_init_part(
    struct part *p, struct xpart *xp, const struct mhd_global_data *mhd_data,
    const double Lsize) {

  mhd_reset_acceleration(p);
  mhd_init_part(p);
}

#endif /* SWIFT_NONE_MHD_H */
