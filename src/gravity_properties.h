/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GRAVITY_PROPERTIES
#define SWIFT_GRAVITY_PROPERTIES

/* Config parameters. */
#include <config.h>

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Forward declarations */
struct cosmology;
struct phys_const;
struct swift_params;

/**
 * @brief Contains all the constants and parameters of the self-gravity scheme
 */
struct gravity_props {

  /* -------------- Softening for the regular particles ---------------- */

  /*! Co-moving softening length for high-res. DM particles at the current
   * redshift.  */
  float epsilon_DM_cur;

  /*! Co-moving softening length for high-res. baryon particles at the current
   * redshift.  */
  float epsilon_baryon_cur;

  /* -------------- Softening for the background DM -------------------- */

  /*! Conversion factor from cbrt of particle mass to softening assuming
   * a constant fraction of the mean inter-particle separation at that mass. */
  float epsilon_background_fac;

  /* -------------- Softening for the neutrino DM ---------------------- */

  /*! Co-moving softening length for neutrino DM particles at the current
   * redshift.  */
  float epsilon_nu_cur;

  /* -------------- Properties of the FFM gravity ---------------------- */

  /*! What MAC are we currently using? */
  int use_advanced_MAC;

  /*! Are we using the adaptive opening angle? (as read from param file) */
  int use_adaptive_tolerance;

  /*! Are we using the Gadget adaptive opening angle? (as read from param file)
   */
  int use_gadget_tolerance;

  /*! Accuracy parameter of the advanced MAC */
  float adaptive_tolerance;

  /*! Tree opening angle (Multipole acceptance criterion) */
  double theta_crit;

  /*! Are we allowing tree gravity below softening? */
  int use_tree_below_softening;

  /*! Are we applying long-range truncation to the forces in the MAC? */
  int consider_truncation_in_MAC;

  /* ------------- Properties of the softened gravity ------------------ */

  /*! Co-moving softening length for for high-res. DM particles */
  float epsilon_DM_comoving;

  /*! Maximal softening length in physical coordinates for the high-res.
   * DM particles */
  float epsilon_DM_max_physical;

  /*! Co-moving softening length for for high-res. baryon particles */
  float epsilon_baryon_comoving;

  /*! Maximal softening length in physical coordinates for the high-res.
   * baryon particles */
  float epsilon_baryon_max_physical;

  /*! Co-moving softening length for for neutrino DM particles */
  float epsilon_nu_comoving;

  /*! Maximal softening length in physical coordinates for the neutrino
   * DM particles */
  float epsilon_nu_max_physical;

  /*! Fraction of the mean inter particle separation corresponding to the
   * co-moving softening length of the low-res. particles (DM + baryons) */
  float mean_inter_particle_fraction_high_res;

  /*! Maximal comoving softening in the case of adaptive softening for gas */
  float max_adaptive_softening;

  /*! Minimal comoving softening in the case of adaptive softening for gas */
  float min_adaptive_softening;

  /* ------------- Properties of the time integration  ----------------- */

  /*! Frequency of tree-rebuild in units of #gpart updates. */
  float rebuild_frequency;

  /*! Fraction of active #gparts needed to trigger a tree-rebuild */
  float rebuild_active_fraction;

  /*! Time integration dimensionless multiplier */
  float eta;

  /* ------------- Properties of the mesh-based gravity ---------------- */

  /*! Periodic long-range mesh side-length */
  int mesh_size;

  /*! Whether mesh is distributed between MPI ranks when we use MPI  */
  int distributed_mesh;

  /*! Whether or not to use local patches rather than
   * direct atomic writes to the mesh when running without MPI */
  int mesh_uses_local_patches;

  /*! Mesh smoothing scale in units of top-level cell size */
  float a_smooth;

  /*! Distance below which the truncated mesh force is Newtonian in units of
   * a_smooth */
  float r_cut_min_ratio;

  /*! Distance above which the truncated mesh force is negligible in units of
   * a_smooth */
  float r_cut_max_ratio;

  /*! Long-range gravity mesh scale. */
  float r_s;

  /*! Inverse of the long-range gravity mesh scale. */
  float r_s_inv;

  /* ------------- Physical constants ---------------------------------- */

  /*! Gravitational constant (in internal units, copied from the physical
   * constants) */
  float G_Newton;
};

void gravity_props_print(const struct gravity_props *p);
void gravity_props_init(struct gravity_props *p, struct swift_params *params,
                        const struct phys_const *phys_const,
                        const struct cosmology *cosmo, const int with_cosmology,
                        const int with_external_potential,
                        const int has_baryons, const int has_DM,
                        const int has_neutrinos, const int is_zoom_simulation,
                        const int periodic, const double dim[3],
                        const int cdim[3]);
void gravity_props_update(struct gravity_props *p,
                          const struct cosmology *cosmo);
void gravity_props_update_MAC_choice(struct gravity_props *p);
#if defined(HAVE_HDF5)
void gravity_props_print_snapshot(hid_t h_grpsph,
                                  const struct gravity_props *p);
#endif

/* Dump/restore. */
void gravity_props_struct_dump(const struct gravity_props *p, FILE *stream);
void gravity_props_struct_restore(struct gravity_props *p, FILE *stream);

#endif /* SWIFT_GRAVITY_PROPERTIES */
