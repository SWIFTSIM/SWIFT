#ifndef SWIFT_POTENTIALS_H
#define SWIFT_POTENTIALS_H

#include "physical_constants_cgs.h"
#include "physical_constants.h"
#include "units.h"

/* External Potential Constants */
struct external_potential {
  struct {
    double x, y, z;
    double mass;
  } point_mass;
};

/* Properties of Point Mass */
#ifdef EXTERNAL_POTENTIAL_POINTMASS

#define External_Potential_X (50000 * PARSEC_IN_CGS / const_unit_length_in_cgs)
#define External_Potential_Y (50000 * PARSEC_IN_CGS / const_unit_length_in_cgs)
#define External_Potential_Z (50000 * PARSEC_IN_CGS / const_unit_length_in_cgs)
#define External_Potential_Mass \
  (1e10 * SOLAR_MASS_IN_CGS / const_unit_mass_in_cgs)


/**
 * @brief Computes the time-step due to the acceleration from a point mass
 *
 * @param phys_cont The physical constants in internal units.
 * @param gp Pointer to the g-particle data.
 */
__attribute__((always_inline))
    INLINE static float external_gravity_pointmass_timestep(
        const struct phys_const* const phys_const,
        const struct gpart* const g) {

  const double G_newton = phys_const->newton_gravity;
  const float dx = g->x[0] - External_Potential_X;
  const float dy = g->x[1] - External_Potential_Y;
  const float dz = g->x[2] - External_Potential_Z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);
  const float drdv = (g->x[0] - External_Potential_X) * (g->v_full[0]) +
                     (g->x[1] - External_Potential_Y) * (g->v_full[1]) +
                     (g->x[2] - External_Potential_Z) * (g->v_full[2]);
  const float dota_x = G_newton * External_Potential_Mass * rinv * rinv * rinv *
                       (-g->v_full[0] + 3.f * rinv * rinv * drdv * dx);
  const float dota_y = G_newton * External_Potential_Mass * rinv * rinv * rinv *
                       (-g->v_full[1] + 3.f * rinv * rinv * drdv * dy);
  const float dota_z = G_newton * External_Potential_Mass * rinv * rinv * rinv *
                       (-g->v_full[2] + 3.f * rinv * rinv * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return 0.03f * sqrtf(a_2 / dota_2);
}


/**
 * @brief Computes the gravitational acceleration of a particle due to a point
 * mass
 *
 * @param phys_cont The physical constants in internal units.
 * @param gp Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_pointmass(
    const struct phys_const* const phys_const, struct gpart* g) {

  const double G_newton = phys_const->newton_gravity;
  const float dx = g->x[0] - External_Potential_X;
  const float dy = g->x[1] - External_Potential_Y;
  const float dz = g->x[2] - External_Potential_Z;
  const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz);

  g->a_grav[0] += -G_newton * External_Potential_Mass * dx * rinv * rinv * rinv;
  g->a_grav[1] += -G_newton * External_Potential_Mass * dy * rinv * rinv * rinv;
  g->a_grav[2] += -G_newton * External_Potential_Mass * dz * rinv * rinv * rinv;
}
#endif /* EXTERNAL_POTENTIAL_POINTMASS */

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
void initPotentialProperties(struct UnitSystem* us,
                             struct external_potential* potential);

#endif /* SWIFT_POTENTIALS_H */
