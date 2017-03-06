/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MULTIPOLE_H
#define SWIFT_MULTIPOLE_H

/* Some standard headers. */
#include <math.h>
#include <string.h>

/* Includes. */
//#include "active.h"
#include "align.h"
#include "const.h"
#include "error.h"
#include "gravity_derivatives.h"
#include "inline.h"
#include "kernel_gravity.h"
#include "minmax.h"
#include "part.h"

#define multipole_align 128

struct acc_tensor {

  double F_000;
};

struct pot_tensor {

  double F_000;
};

struct multipole {

  float mass;

  /*! Bulk velocity */
  float vel[3];
};

/**
 * @brief The multipole expansion of a mass distribution.
 */
struct gravity_tensors {

  /*! Linking pointer for "memory management". */
  struct gravity_tensors *next;

  /*! Centre of mass of the matter dsitribution */
  double CoM[3];

  /*! The actual content */
  struct {

    /*! Multipole mass */
    struct multipole m_pole;

    /*! Field tensor for acceleration along x */
    struct acc_tensor a_x;

    /*! Field tensor for acceleration along y */
    struct acc_tensor a_y;

    /*! Field tensor for acceleration along z */
    struct acc_tensor a_z;

    /*! Field tensor for the potential */
    struct pot_tensor pot;

#ifdef SWIFT_DEBUG_CHECKS

    /* Total mass this gpart interacted with */
    double mass_interacted;

#endif
  };
} SWIFT_STRUCT_ALIGN;

/**
 * @brief Reset the data of a #multipole.
 *
 * @param m The #multipole.
 */
INLINE static void gravity_reset(struct gravity_tensors *m) {

  /* Just bzero the struct. */
  bzero(m, sizeof(struct gravity_tensors));
}

INLINE static void gravity_field_tensors_init(struct gravity_tensors *m) {

  bzero(&m->a_x, sizeof(struct acc_tensor));
  bzero(&m->a_y, sizeof(struct acc_tensor));
  bzero(&m->a_z, sizeof(struct acc_tensor));
  bzero(&m->pot, sizeof(struct pot_tensor));

#ifdef SWIFT_DEBUG_CHECKS
  m->mass_interacted = 0.;
#endif
}

/**
 * @brief Adds field tensrs to other ones (i.e. does la += lb).
 *
 * @param la The gravity tensors to add to.
 * @param lb The gravity tensors to add.
 */
INLINE static void gravity_field_tensors_add(struct gravity_tensors *la,
                                             const struct gravity_tensors *lb) {
  la->a_x.F_000 += lb->a_x.F_000;
  la->a_y.F_000 += lb->a_y.F_000;
  la->a_z.F_000 += lb->a_z.F_000;

#ifdef SWIFT_DEBUG_CHECKS
  if (lb->mass_interacted == 0.f) error("Adding tensors that did not interact");
  la->mass_interacted += lb->mass_interacted;
#endif
}

/**
 * @brief Prints the content of a #multipole to stdout.
 *
 * Note: Uses directly printf(), not a call to message().
 *
 * @param m The #multipole to print.
 */
INLINE static void gravity_multipole_print(const struct multipole *m) {

  // printf("CoM= [%12.5e %12.5e %12.5e\n", m->CoM[0], m->CoM[1], m->CoM[2]);
  printf("Mass= %12.5e\n", m->mass);
  printf("Vel= [%12.5e %12.5e %12.5e\n", m->vel[0], m->vel[1], m->vel[2]);
}

/**
 * @brief Adds a #multipole to another one (i.e. does ma += mb).
 *
 * @param ma The multipole to add to.
 * @param mb The multipole to add.
 */
INLINE static void gravity_multipole_add(struct multipole *ma,
                                         const struct multipole *mb) {

  const float mass = ma->mass + mb->mass;
  const float imass = 1.f / mass;

  /* Add the bulk velocities */
  ma->vel[0] = (ma->vel[0] * ma->mass + mb->vel[0] * mb->mass) * imass;
  ma->vel[1] = (ma->vel[1] * ma->mass + mb->vel[1] * mb->mass) * imass;
  ma->vel[2] = (ma->vel[2] * ma->mass + mb->vel[2] * mb->mass) * imass;

  /* Add the masses */
  ma->mass = mass;
}

/**
 * @brief Verifies whether two #multipole's are equal or not.
 *
 * @param ma The first #multipole.
 * @param mb The second #multipole.
 * @param tolerance The maximal allowed difference for the fields.
 * @return 1 if the multipoles are equal 0 otherwise.
 */
INLINE static int gravity_multipole_equal(const struct multipole *ma,
                                          const struct multipole *mb,
                                          double tolerance) {

  /* Check CoM */
  /* if (fabs(ma->CoM[0] - mb->CoM[0]) / fabs(ma->CoM[0] + mb->CoM[0]) >
   * tolerance) */
  /*   return 0; */
  /* if (fabs(ma->CoM[1] - mb->CoM[1]) / fabs(ma->CoM[1] + mb->CoM[1]) >
   * tolerance) */
  /*   return 0; */
  /* if (fabs(ma->CoM[2] - mb->CoM[2]) / fabs(ma->CoM[2] + mb->CoM[2]) >
   * tolerance) */
  /*   return 0; */

  /* Check bulk velocity (if non-zero)*/
  if (fabsf(ma->vel[0] + mb->vel[0]) > 0.f &&
      fabsf(ma->vel[0] - mb->vel[0]) / fabsf(ma->vel[0] + mb->vel[0]) >
          tolerance)
    return 0;
  if (fabsf(ma->vel[1] + mb->vel[1]) > 0.f &&
      fabsf(ma->vel[1] - mb->vel[1]) / fabsf(ma->vel[1] + mb->vel[1]) >
          tolerance)
    return 0;
  if (fabsf(ma->vel[2] + mb->vel[2]) > 0.f &&
      fabsf(ma->vel[2] - mb->vel[2]) / fabsf(ma->vel[2] + mb->vel[2]) >
          tolerance)
    return 0;

  /* Check mass */
  if (fabsf(ma->mass - mb->mass) / fabsf(ma->mass + mb->mass) > tolerance)
    return 0;

  /* All is good */
  return 1;
}

/**
 * @brief Drifts a #multipole forward in time.
 *
 * @param m The #multipole.
 * @param dt The drift time-step.
 */
INLINE static void gravity_multipole_drift(struct gravity_tensors *m,
                                           double dt) {

  /* Move the whole thing according to bulk motion */
  m->CoM[0] += m->m_pole.vel[0];
  m->CoM[1] += m->m_pole.vel[1];
  m->CoM[2] += m->m_pole.vel[2];
}

/**
 * @brief Applies the forces due to particles j onto particles i directly.
 *
 * @param gparts_i The #gpart to update.
 * @param gcount_i The number of particles to update.
 * @param gparts_j The #gpart that source the gravity field.
 * @param gcount_j The number of sources.
 */
INLINE static void gravity_P2P(struct gpart *gparts_i, int gcount_i,
                               const struct gpart *gparts_j, int gcount_j) {}

/**
 * @brief Constructs the #multipole of a bunch of particles around their
 * centre of mass.
 *
 * Corresponds to equation (28c).
 *
 * @param m The #multipole (content will  be overwritten).
 * @param gparts The #gpart.
 * @param gcount The number of particles.
 */
INLINE static void gravity_P2M(struct gravity_tensors *m,
                               const struct gpart *gparts, int gcount) {

#if const_gravity_multipole_order >= 2
#error "Implementation of P2M kernel missing for this order."
#endif

  /* Temporary variables */
  double mass = 0.0;
  double com[3] = {0.0, 0.0, 0.0};
  float vel[3] = {0.f, 0.f, 0.f};

  /* Collect the particle data. */
  for (int k = 0; k < gcount; k++) {
    const float m = gparts[k].mass;

    mass += m;
    com[0] += gparts[k].x[0] * m;
    com[1] += gparts[k].x[1] * m;
    com[2] += gparts[k].x[2] * m;
    vel[0] += gparts[k].v_full[0] * m;
    vel[1] += gparts[k].v_full[1] * m;
    vel[2] += gparts[k].v_full[2] * m;
  }

  const double imass = 1.0 / mass;

  /* Store the data on the multipole. */
  m->m_pole.mass = mass;
  m->CoM[0] = com[0] * imass;
  m->CoM[1] = com[1] * imass;
  m->CoM[2] = com[2] * imass;
  m->m_pole.vel[0] = vel[0] * imass;
  m->m_pole.vel[1] = vel[1] * imass;
  m->m_pole.vel[2] = vel[2] * imass;
}

/**
 * @brief Creates a copy of #multipole shifted to a new location.
 *
 * Corresponds to equation (28d).
 *
 * @param m_a The #multipole copy (content will  be overwritten).
 * @param m_b The #multipole to shift.
 * @param pos_a The position to which m_b will be shifted.
 * @param pos_b The current postion of the multipole to shift.
 * @param periodic Is the calculation periodic ?
 */
INLINE static void gravity_M2M(struct multipole *m_a,
                               const struct multipole *m_b,
                               const double pos_a[3], const double pos_b[3],
                               int periodic) {

  m_a->mass = m_b->mass;

  m_a->vel[0] = m_b->vel[0];
  m_a->vel[1] = m_b->vel[1];
  m_a->vel[2] = m_b->vel[2];
}
/**
 * @brief Compute the field tensors due to a multipole.
 *
 * Corresponds to equation (28b).
 *
 * @param l_a The field tensor to compute.
 * @param m_b The multipole creating the field.
 * @param pos_a The position of the field tensor.
 * @param pos_b The position of the multipole.
 * @param periodic Is the calculation periodic ?
 */
INLINE static void gravity_M2L(struct gravity_tensors *l_a,
                               const struct multipole *m_b,
                               const double pos_a[3], const double pos_b[3],
                               int periodic) {

  double dx, dy, dz;
  if (periodic) {
    dx = box_wrap(pos_a[0] - pos_b[0], 0., 1.);
    dy = box_wrap(pos_a[1] - pos_b[1], 0., 1.);
    dz = box_wrap(pos_a[2] - pos_b[2], 0., 1.);
  } else {
    dx = pos_a[0] - pos_b[0];
    dy = pos_a[1] - pos_b[1];
    dz = pos_a[2] - pos_b[2];
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  const double r_inv = 1. / sqrt(r2);

  /* 1st order multipole term */
  l_a->a_x.F_000 += D_100(dx, dy, dz, r_inv) * m_b->mass;
  l_a->a_y.F_000 += D_010(dx, dy, dz, r_inv) * m_b->mass;
  l_a->a_z.F_000 += D_001(dx, dy, dz, r_inv) * m_b->mass;

#ifdef SWIFT_DEBUG_CHECKS
  l_a->mass_interacted += m_b->mass;
#endif
}

/**
 * @brief Creates a copy of #acc_tensor shifted to a new location.
 *
 * Corresponds to equation (28e).
 *
 * @param l_a The #acc_tensor copy (content will  be overwritten).
 * @param l_b The #acc_tensor to shift.
 * @param pos_a The position to which m_b will be shifted.
 * @param pos_b The current postion of the multipole to shift.
 * @param periodic Is the calculation periodic ?
 */
INLINE static void gravity_L2L(struct gravity_tensors *l_a,
                               const struct gravity_tensors *l_b,
                               const double pos_a[3], const double pos_b[3],
                               int periodic) {

  l_a->a_x.F_000 = l_b->a_x.F_000;
  l_a->a_y.F_000 = l_b->a_y.F_000;
  l_a->a_z.F_000 = l_b->a_z.F_000;

#ifdef SWIFT_DEBUG_CHECKS
  if (l_b->mass_interacted == 0.f)
    error("Shifting tensors that did not interact");
  l_a->mass_interacted = l_b->mass_interacted;
#endif
}

/**
 * @brief Applies the  #acc_tensor to a  #gpart.
 *
 * Corresponds to equation (28a).
 */
INLINE static void gravity_L2P(const struct gravity_tensors *l,
                               struct gpart *gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (l->mass_interacted == 0.f) error("Interacting with empty field tensor");
  gp->mass_interacted += l->mass_interacted;
#endif

  // message("a=[%e %e %e]", l->a_x.F_000, l->a_y.F_000, l->a_z.F_000);

  /* 0th order interaction */
  gp->a_grav[0] += l->a_x.F_000;
  gp->a_grav[1] += l->a_y.F_000;
  gp->a_grav[2] += l->a_z.F_000;
}

#endif /* SWIFT_MULTIPOLE_H */
