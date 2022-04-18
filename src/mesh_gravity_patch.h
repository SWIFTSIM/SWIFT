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

#ifndef SWIFT_MESH_GRAVITY_PATCH_H
#define SWIFT_MESH_GRAVITY_PATCH_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "align.h"
#include "error.h"
#include "inline.h"

/* Forward declarations */
struct cell;

/**
 * @brief Data structure for a patch of mesh covering a cell
 */
struct pm_mesh_patch {

  /*! Size of the full mesh */
  int N;

  /*! Inverse of the mesh cell size */
  double fac;

  /*! Minimum of coordinate range into which particles should be wrapped */
  double wrap_min[3];

  /*! Maximum of coordinate range into which particles should be wrapped */
  double wrap_max[3];

  /*! Size of the mesh in each dimension */
  int mesh_size[3];

  /*! Minimum integer coordinate of the mesh in each dimension */
  int mesh_min[3];

  /*! Maximum integer coordinate of the mesh in each dimension */
  int mesh_max[3];

  /*! Pointer to the mesh data */
  double *mesh;
};

void pm_mesh_patch_init(struct pm_mesh_patch *patch, const struct cell *cell,
                        const int N, const double fac, const double dim[3],
                        const int boundary_size);

void pm_mesh_patch_zero(struct pm_mesh_patch *patch);

void pm_mesh_patch_clean(struct pm_mesh_patch *patch);

/**
 * @brief Return the array index in the patch corresponding to
 * coordinates in the full mesh
 *
 * @param patch Pointer to the patch
 * @param i Integer x coordinate in the full mesh
 * @param j Integer y coordinate in the full mesh
 * @param k Integer z coordinate in the full mesh
 *
 */
__attribute__((always_inline)) INLINE static int pm_mesh_patch_index(
    const struct pm_mesh_patch *patch, const int i, const int j, const int k) {

#ifdef SWIFT_DEBUG_CHECKS
  if (i < 0 || i >= patch->mesh_size[0])
    error("Coordinate in local mesh out of range!");
  if (j < 0 || j >= patch->mesh_size[1])
    error("Coordinate in local mesh out of range!");
  if (k < 0 || k >= patch->mesh_size[2])
    error("Coordinate in local mesh out of range!");
#endif

  return (i * patch->mesh_size[1] * patch->mesh_size[2]) +
         (j * patch->mesh_size[2]) + k;
}

/**
 * @brief Cloud in cell evaluation of the mesh patch
 *
 * @param patch Pointer to the patch
 * @param i Integer x coordinate in the mesh patch
 * @param j Integer y coordinate in the mesh patch
 * @param k Integer z coordinate in the mesh patch
 * @param tx CIC parameter in the x direction
 * @param ty CIC parameter in the y direction
 * @param tz CIC parameter in the z direction
 * @param dx CIC parameter in the x direction
 * @param dy CIC parameter in the y direction
 * @param dz CIC parameter in the z direction
 *
 */
__attribute__((always_inline)) INLINE static double pm_mesh_patch_CIC_get(
    const struct pm_mesh_patch *patch, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz) {

  /* Remind the compiler that the arrays are nicely aligned */
  swift_declare_aligned_ptr(const double, mesh, patch->mesh,
                            SWIFT_CACHE_ALIGNMENT);

  double temp;

  /* Classic 3D CIC */
  temp = mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 0)] * tx * ty * tz;
  temp += mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 1)] * tx * ty * dz;
  temp += mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 0)] * tx * dy * tz;
  temp += mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 1)] * tx * dy * dz;
  temp += mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 0)] * dx * ty * tz;
  temp += mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 1)] * dx * ty * dz;
  temp += mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 0)] * dx * dy * tz;
  temp += mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 1)] * dx * dy * dz;
  return temp;
}

/**
 * @brief Cloud in cell assignment to the mesh patch
 *
 * @param patch Pointer to the patch
 * @param i Integer x coordinate in the mesh patch
 * @param j Integer y coordinate in the mesh patch
 * @param k Integer z coordinate in the mesh patch
 * @param tx CIC parameter in the x direction
 * @param ty CIC parameter in the y direction
 * @param tz CIC parameter in the z direction
 * @param dx CIC parameter in the x direction
 * @param dy CIC parameter in the y direction
 * @param dz CIC parameter in the z direction
 * @param value The value to set
 */
__attribute__((always_inline)) INLINE static void pm_mesh_patch_CIC_set(
    const struct pm_mesh_patch *patch, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz, const double value) {

  /* Remind the compiler that the arrays are nicely aligned */
  swift_declare_aligned_ptr(double, mesh, patch->mesh, SWIFT_CACHE_ALIGNMENT);

  /* Classic 3D CIC */
  mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 0)] += value * tx * ty * tz;
  mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 1)] += value * tx * ty * dz;
  mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 0)] += value * tx * dy * tz;
  mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 1)] += value * tx * dy * dz;
  mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 0)] += value * dx * ty * tz;
  mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 1)] += value * dx * ty * dz;
  mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 0)] += value * dx * dy * tz;
  mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 1)] += value * dx * dy * dz;
}

#endif
