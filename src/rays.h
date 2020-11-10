/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Evgenii Chaikin (chaikin@strw.leidenuniv.nl)
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
#ifndef SWIFT_RAYS_H
#define SWIFT_RAYS_H

/**
 * @brief For a given feedback scheme, sets all fields in the ray struct
 * array at their default values
 *
 * @param rays An array of ray structs
 * @param max_number_of_rays Maximum number of rays in the feedback scheme
 */
__attribute__((always_inline)) INLINE static void ray_init(
    struct ray_data *rays, const int max_number_of_rays) {

  /* Set all fields in the ray struct at their default values */
  for (int i = 0; i < max_number_of_rays; i++) {
    rays[i].id_min_length = -1;
    rays[i].min_length = FLT_MAX;
    rays[i].mass = 0.f;
  }
}

/**
 * @brief Reset the ray(s) that point to the particle with a given id
 *
 * @param rays An array of ray structs
 * @param max_number_of_rays Maximum number of rays in the feedback scheme
 * @param gas_part_id The id of the gas particle the ray is pointing to
 */
__attribute__((always_inline)) INLINE static void ray_reset_part_id(
    struct ray_data *rays, const int max_number_of_rays,
    const long long gas_part_id) {

  /* Loop over all rays the provided stellar/BH particle has */
  for (int i = 0; i < max_number_of_rays; i++) {

    /* Find the ray carrying the gas-partile id equal to gas_part_id */
    if (rays[i].id_min_length == gas_part_id) {

      /* Reset all fields in this ray */
      rays[i].id_min_length = -1;
      rays[i].min_length = FLT_MAX;
      rays[i].mass = 0.f;
    }
  }
}

/**
 * @brief For a given feedback scheme, sets all fields in the ray_extra
 * struct array at their default values
 *
 * @param rays_ext An array of ray_extra structs
 * @param max_number_of_rays Maximum number of rays in the feedback scheme
 */
__attribute__((always_inline)) INLINE static void ray_extra_init(
    struct ray_data_extra *rays_ext, const int max_number_of_rays) {

  /* Set all fields in the ray struct at their default values */
  for (int i = 0; i < max_number_of_rays; i++) {
    for (int j = 0; j < 3; j++) {
      rays_ext[i].x[j] = 0.f;
      rays_ext[i].v[j] = 0.f;
    }
  }
}

/**
 * @brief Returns the arc length on a sphere of radius r_sphere
 * between two points with angular coordinates (theta_1, phi_1) and (theta_2,
 * phi_2)
 *
 * @param theta_1 Polar angle of point 1; \theta \in [-\pi/2, \pi/2]
 * @param phi_1 Azimuthal angle of point 1; \phi \in [-\pi, \pi)
 * @param theta_2 Polar angle of point 2; \theta \in [-\pi/2, \pi/2]
 * @param phi_2 Azimuthal angle of point 2; \phi \in [-\pi, \pi)
 * @param r_sphere Radius of the sphere on which the arc length between the two
 * points is computed
 */
__attribute__((always_inline)) INLINE static float ray_arclength(
    const double theta_1, const double phi_1, const double theta_2,
    const double phi_2, const float r_sphere) {

  const double delta_theta = theta_2 - theta_1;
  const double delta_phi = phi_2 - phi_1;

  const double sine_theta = sin(delta_theta / 2.0);
  const double sine_phi = sin(delta_phi / 2.0);

  const float arc_length =
      asin(sqrt(sine_theta * sine_theta +
                cos(theta_1) * cos(theta_2) * sine_phi * sine_phi));

  return 2.f * r_sphere * arc_length;
}

/**
 * @brief For a given ray, minimise the arc length between the stellar/BH
 * particle and the gas neighbours in the isotropic feedback and store relevant
 * properties of the gas particle that has the minimum arc length with the ray
 *
 * @param dx Comoving vector separating both particles (si - pj)
 * @param r Comoving distance between the two particles
 * @param ray Ray data
 * @param mirror_particle_switch Mirror particle switch
 * @param gas_part_id ID of the gas particle
 * @param rand_theta_gen Random number to generate \theta_ray
 * @param rand_phi_gen Random number to generate \phi_ray
 * @param m Gas particle mass
 * @param ray_ext Extra ray data
 * @param v Gas particle velocity
 */
__attribute__((always_inline)) INLINE static void ray_minimise_arclength(
    const float *dx, const float r, struct ray_data *ray,
    const int mirror_particle_switch, const long long gas_part_id,
    const double rand_theta_gen, const double rand_phi_gen, const float m,
    struct ray_data_extra *ray_ext, const float *v) {

  /* Angular coordinates of the particle with respect to the star/BH */
  const double theta_j = acos(-dx[2] / r);
  const double phi_j = atan2(-dx[1], -dx[0]);

  /* Transform the random number from [0,1[ to [-1, 1[ */
  const double cos_theta_ray = 2. * rand_theta_gen - 1.;

  /* Get the \theta angle */
  double theta_ray = acos(cos_theta_ray);

  /* Transform the random number from [0,1[ to [-pi, pi[ */
  double phi_ray = 2.0 * M_PI * rand_phi_gen - M_PI;

  /* Flip the angles if it is a mirror particle */
  if (mirror_particle_switch == 1) {
    theta_ray = M_PI - theta_ray;
    phi_ray = phi_ray - copysign(M_PI, phi_ray);
  }

  /* Calculate the arc length on a unit sphere between the gas particle
   * and this ray. Note that shift by -\pi/2 is required for the \theta's */
  const float new_length = ray_arclength(
      theta_ray - M_PI_2, phi_ray, theta_j - M_PI_2, phi_j, /*r_sphere=*/1.f);

  /* If the new arc length is smaller than the older value, store
   * the new one alongside the particle id and its other relevant
   * properties that will be used later */
  if (new_length < ray->min_length) {

    ray->min_length = new_length;
    ray->id_min_length = gas_part_id;
    ray->mass = m;

    /* In stellar feedback, we do kicks so we need to store additional data */
    if (mirror_particle_switch != -1) {

      /* Position and velocities are neIf prob <1. and we are running with eded
       * in SNII kinetic feedback to exactly conserve momentum and energy.
       * That's because in a pair of two particles, the first one needs to know
       * the properties of the other one, and vice versa */

      ray_ext->x[0] = -dx[0];
      ray_ext->x[1] = -dx[1];
      ray_ext->x[2] = -dx[2];

      ray_ext->v[0] = v[0];
      ray_ext->v[1] = v[1];
      ray_ext->v[2] = v[2];
    }
  }
}

/**
 * @brief Compute the kick velocity in stellar kinetic feedback
 *
 * @param v_kick The kick velocty vector to be computed
 * @param v_kick_abs Modulus of the kick velocity vector
 * @param ray_ext_true Extra ray data of the true gas particle
 * @param ray_ext_mirr Extra ray data of the mirror gas particle
 * @param mirror_particle_switch Mirror particle switch
 * @param energy_pair SN energy per pair
 * @param cosmo The cosmological model
 * @param current_mass Current mass of the gas particle
 * @param v_star Velocity of the stellar particle
 * @param rand_theta_gen Random number to generate \theta_ray
 * @param rand_phi_gen Random number to generate \phi_ray
 * @param mass_true Unaffected mass of the true particle
 * @param mass_mirror Unaffected mass of the mirror particle
 */
__attribute__((always_inline)) INLINE static void
ray_kinetic_feedback_compute_kick_velocity(
    float *v_kick, float *v_kick_abs, const struct ray_data_extra *ray_ext_true,
    const struct ray_data_extra *ray_ext_mirr, const int mirror_particle_switch,
    const double energy_pair, const struct cosmology *cosmo,
    const double current_mass, const float *v_star, const double rand_theta_gen,
    const double rand_phi_gen, const double mass_true,
    const double mass_mirror) {

  /* Transform the random number from [0,1[ to [-1, 1[ */
  const double cos_theta_ray = 2. * rand_theta_gen - 1.;

  /* Transform the random number from [0,1[ to [-pi, pi[ */
  const double phi_ray = 2.0 * M_PI * rand_phi_gen - M_PI;

  /* To conserve angular momentum -- along with linear momentum and energy
   * -- we need to kick the two particles in the pair along the line that
   * connects these two particles. In the code below, we start with
   * finding the equation defining this line. The line equation(-s) reads:
   * (x-x_0)/a = (y-y_0)/b = (z-z_0)/c . We thus need to find the
   * coefficients a, b and c. Since we already know two points on this
   * line, which are the particles' positions, this is a straightforward
   * exercise to do. */

  /* Coefficients defining the line connecting the two particles */
  double a, b, c;

  /* Modulus of the vector {a, b, c} that defines the line
   * connecting the original particle with the mirror particle */
  double normalisation;

  /* The kick angles (with respect to the stellar particle) */
  double phi_kick, theta_kick;

  /* If we have no degeneracy along the x coordinate, we can divide
   * by (x_2 - x_1) to compute the coefficients b and c */
  if (ray_ext_true->x[0] != ray_ext_mirr->x[0]) {

    /* Due to the freedom in normalisation, the coefficient a can be
     * set to unity */
    a = 1.0;

    /* Compute b = b/1.0 = b/a = (y_2 - y_1) / (x_2 - x_1) */
    b = (ray_ext_true->x[1] - ray_ext_mirr->x[1]) /
        (ray_ext_true->x[0] - ray_ext_mirr->x[0]);

    /* Compute c = c/1.0 = c/a = (z_2 - z_1) / (x_2 - x_1) */
    c = (ray_ext_true->x[2] - ray_ext_mirr->x[2]) /
        (ray_ext_true->x[0] - ray_ext_mirr->x[0]);

    /* Compute the modulus of {a, b, c} and the kick angles */
    normalisation = sqrt(1.0 + b * b + c * c);
    phi_kick = atan2(b, 1.0);
    theta_kick = acos(c / normalisation);
  }

  /* If x2 and x1 are the same, we cannot divide by (x_2 - x_1)
   * If we have no degeneracy along the y coordinate, we can divide
   * by (y_2 - y_1) to compute the coefficient c */
  else if (ray_ext_true->x[1] != ray_ext_mirr->x[1]) {

    /* Since x_2 = x_1, the line is perpendicular to the x axis so
     * that a = 0 and \phi_kick = \pi/2 */
    a = 0.0;

    /* Due to the freedom in normalisation, the coefficient b can be
     * set to unity */
    b = 1.0;

    /* Compute c = c/1.0 = c/b = (z_2 - z_1) / (y_2 - y_1) */
    c = (ray_ext_true->x[2] - ray_ext_mirr->x[2]) /
        (ray_ext_true->x[1] - ray_ext_mirr->x[1]);

    /* Compute the modulus of {a ,b ,c} and the kick angles */
    normalisation = sqrt(1.0 + c * c);
    phi_kick = M_PI_2;
    theta_kick = acos(c / normalisation);
  }

  /* y_2 - y_1 = 0 as well as x_2 - x_1 = 0 */
  else {

    /* The line is parallel to the z axis, the coefficients a and b
     * are zero. The angles \phi_kick = 0.0, and \theta_kick = 0.0 */
    a = 0.0;
    b = 0.0;

    /* Due to the freedom in normalisation, the coefficient c can be
     * set to unity */
    c = 1.0;

    /* Compute the modulus of {a ,b ,c} and the kick angles */
    normalisation = 1.0;
    phi_kick = 0.0;
    theta_kick = 0.0;
  }

  /* By now, by computing the equation of the line connecting the two
   * particles in the pair, we have found the direction of the kick.
   * However, for a given particle in the pair, we know this direction
   * only up to a sign, i.e. we as of yet do not know whether we should
   * kick the particle "up" or "down" the line. In order to get the sign
   * right we need to do some additional calculations. Namely, we need
   * to compare the direction of the line we found with that of the
   * original and mirror rays. In total, we have four cases:
   *
   * 1 and 2: IF the (theta_kick, phi_kick) direction is closest to the
   * direction of the original (mirror) ray AND we have caught the
   * original (mirror) particle THEN kick along the (theta_kick,
   * phi_kick) direction. 3 and 4: IF the (theta_kick, phi_kick)
   * direction is closest to the direction of the original (mirror) ray
   * BUT we have caught the mirror (original) particle THEN kick along
   * the opposite (\pi - theta_kick, phi_kick - pi*sign(phi_kick))
   * direction */

  /* Thus far we've only got cos\theta. But in order to compute arc
  lengths used in the algorithm described above we need \theta */
  const double theta_ray = acos(cos_theta_ray);

  /* Compute the arc lengh between the original ray and the line
   * defining the direction of the kicks. We shift \theta coordinates
   * by -\pi/2 because we have \theta \in [0, \pi[ while the function
   * wants \theta \in [-\pi/2, \pi/2[*/
  const double arclength_true =
      ray_arclength(theta_kick - M_PI_2, phi_kick, theta_ray - M_PI_2, phi_ray,
                    /*r_sphere=*/1.f);

  /* \theta and \phi angles of the mirror ray, pointing in the
   * direction opposite to the original one. We need those to compute
   * the arc length below */
  const double theta_ray_mirror = M_PI - theta_ray;
  const double phi_ray_mirror = phi_ray - copysign(M_PI, phi_ray);

  /* Compute the arc lengh between the mirror ray and the line defining
   * the direction of the kicks */
  const double arclength_mirror =
      ray_arclength(theta_kick - M_PI_2, phi_kick, theta_ray_mirror - M_PI_2,
                    phi_ray_mirror, /*r_sphere=*/1.f);

  /* Find the minimum arc length between the two arc lengths computed
   * above */
  const double arclength_min = min(arclength_true, arclength_mirror);

  /* Find out whether the minimal arc length is that with the original
  ray or the mirror ray */
  const int mirror_ray_switch = (arclength_min == arclength_mirror);

  /* Compute the sign of the kick depending on which case we are in: 1,
   * 2, 3 or 4 */
  const double kick_sign =
      (mirror_particle_switch == mirror_ray_switch) ? 1.0 : -1.0;

  /* Compute the normal vector of the kick */
  const double n_kick[3] = {kick_sign * a / normalisation,   // x,
                            kick_sign * b / normalisation,   // y
                            kick_sign * c / normalisation};  // z

  /* Compute mass weights that are used below */
  /* Note that mass_true and current_mass may differ if, for example,
   * at this time-step the considered here gas particle has already
   * received some ejecta mass from another star. To compute the weights,
   * we need to use the "unaffected" mass, mass_true, because it does not
   * change when viewed from either of the particles in the pair. At the
   * same time, in the definition of mass_weight, we do divide by
   * current_mass and not by mass_true since this is the weight for the
   * velocity v, which -- for the considered gas particle -- is related to
   * the momentum p as p = current_mass * v and not via mass_true. The
   * differences between current_mass and mass_true, if ever exist, are
   * obviously very minor. Nonetheless, these are not negligible because
   * we require the exact conservation of linear momentum. */
  const double m_alpha =
      sqrt(mass_true * mass_mirror) / (mass_true + mass_mirror);
  const double mass_weight = sqrt(mass_true * mass_mirror) / current_mass;

  /* Relative velocity between the gas particle and the stellar particle
   */
  double v_gas_star[3] = {ray_ext_true->v[0] - v_star[0],
                          ray_ext_true->v[1] - v_star[1],
                          ray_ext_true->v[2] - v_star[2]};

  /* Relative velocity between the mirror gas particle and the stellar
   * particle */
  double v_gas_mirror_star[3] = {ray_ext_mirr->v[0] - v_star[0],
                                 ray_ext_mirr->v[1] - v_star[1],
                                 ray_ext_mirr->v[2] - v_star[2]};

  /* Divide the velocities by the cosmic scale factor
  to get peculiar velocities in proper coordinates */
  for (int j = 0; j < 3; j++) {
    v_gas_star[j] /= cosmo->a;
    v_gas_mirror_star[j] /= cosmo->a;
  }

  /* Compute scalar product between v_gas_star and n */
  const double v_cos_theta = v_gas_star[0] * n_kick[0] +
                             v_gas_star[1] * n_kick[1] +
                             v_gas_star[2] * n_kick[2];

  /* Compute scalar product between v_gas_mirror_star and n */
  const double v_mirror_cos_theta = v_gas_mirror_star[0] * n_kick[0] +
                                    v_gas_mirror_star[1] * n_kick[1] +
                                    v_gas_mirror_star[2] * n_kick[2];

  /* Compute the characteristic kick velocity corresponding to the kinetic
   * energy per pair. */
  const double SNII_delta_v =
      sqrt(2.0 * energy_pair / (mass_true + mass_mirror));

  /* Compute the correction to the energy and momentum due to relative
   * star-gas motion. If it is the mirror particle multiply by the minus
   * sign. If there is no correction then alpha = 0 and beta = 1 */
  const double correction_sign = (mirror_particle_switch) ? -1.0 : 1.0;

  const double alpha = correction_sign * m_alpha *
                       (v_cos_theta - v_mirror_cos_theta) / SNII_delta_v;
  const double beta = sqrt(alpha * alpha + 1.0) - alpha;

  /* Compute the physical kick velocity in internal units */
  (*v_kick_abs) = SNII_delta_v * mass_weight * beta;

  /* The kick velocity vector */
  v_kick[0] = (*v_kick_abs) * n_kick[0];
  v_kick[1] = (*v_kick_abs) * n_kick[1];
  v_kick[2] = (*v_kick_abs) * n_kick[2];
}

/**
 * @brief Minimises the distance between the stellar/BH particle and
 * its neighbours. In the resulting ray array, the gas neighbours will
 * be sorted based on their separation from the star/BH. The left-most
 * element will store the information of the particle with the smallest
 * distance to the star/BH, and the right-most with the largest.
 *
 * @param r Comoving distance between the two particles
 * @param ray Ray data
 * @param N_ray_arr Size of ray array
 * @param gas_part_id ID of the gas particle
 * @param m Gas particle mass
 */
__attribute__((always_inline)) INLINE static void ray_minimise_distance(
    const float r, struct ray_data *ray, const int N_ray_arr,
    const long long gas_part_id, const float m) {

  /* Id of the element of the ray array that we want to update */
  int elem_update_id = -1;

  /* Loop from left to right. Find the left-most element whose
   * current min length is larger than r. Note that N_ray_arr is
   * equal to min(how many times this function has been called
   * at this time-step, the maximum number of rays per particle) */
  for (int i = 0; i < N_ray_arr; i++) {
    if (r < ray[i].min_length) {
      elem_update_id = i;
      break;
    }
  }

  /* If found something to update */
  if (elem_update_id != -1) {

    /* Loop from right to left. Before updating elem_update_id, move
     * everything to the right by one element. */
    for (int i = N_ray_arr - 2; i >= elem_update_id; i--) {
      ray[i + 1].min_length = ray[i].min_length;
      ray[i + 1].id_min_length = ray[i].id_min_length;
      ray[i + 1].mass = ray[i].mass;
    }

    /* Update elem_update_id */
    ray[elem_update_id].min_length = r;
    ray[elem_update_id].id_min_length = gas_part_id;
    ray[elem_update_id].mass = m;
  }
}

#endif /* SWIFT_RAYS_H */
