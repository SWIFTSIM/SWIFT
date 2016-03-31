#ifndef SWIFT_GRAVITY_CONST_H
#define SWIFT_GRAVITY_CONST_H

/* External Potential Constants */

/* Properties of Point Mass */
#ifdef EXTERNAL_POTENTIAL_POINTMASS
#define External_Potential_X  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Y  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Z  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Mass (1e10 * SOLAR_MASS_IN_CGS / const_unit_mass_in_cgs)
#endif


__attribute__((always_inline)) INLINE static float external_gravity_pointmass_timestep(
    struct gpart* g) {

  /* Currently no limit is imposed */
  const float dx   = g->x[0]-External_Potential_X;
  const float dy   = g->x[1]-External_Potential_Y;
  const float dz   = g->x[2]-External_Potential_Z;
  const float rinv = 1.f / sqrtf(dx*dx + dy*dy + dz*dz);
  const float drdv = (g->x[0]-External_Potential_X) * (g->v_full[0]) + (g->x[1]-External_Potential_Y) * (g->v_full[1]) + (g->x[2]-External_Potential_Z) * (g->v_full[2]); 
  const float dota_x = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[0] + 3.f * rinv * rinv * drdv * dx);
  const float dota_y = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[1] + 3.f * rinv * rinv * drdv * dy);
  const float dota_z = const_G * External_Potential_Mass * rinv * rinv * rinv * (-g->v_full[2] + 3.f * rinv * rinv * drdv * dz);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2    = g->a_grav_external[0] * g->a_grav_external[0] + g->a_grav_external[1] * g->a_grav_external[1] + g->a_grav_external[2] * g->a_grav_external[2];
  
  return 0.03f * sqrtf(a_2/dota_2);
}

__attribute__((always_inline)) INLINE static void external_gravity_pointmass(struct gpart *g)
{
  const float dx   = g->x[0]-External_Potential_X;
  const float dy   = g->x[1]-External_Potential_Y;
  const float dz   = g->x[2]-External_Potential_Z;
  const float rinv = 1.f / sqrtf(dx*dx + dy*dy + dz*dz);
  

  g->a_grav_external[0] += - const_G *  External_Potential_Mass * dx * rinv * rinv * rinv;
  g->a_grav_external[1] += - const_G *  External_Potential_Mass * dy * rinv * rinv * rinv;
  g->a_grav_external[2] += - const_G *  External_Potential_Mass * dz * rinv * rinv * rinv;
}



#endif /* SWIFT_GRAVITY_CONST_H */

