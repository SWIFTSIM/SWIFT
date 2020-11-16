/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 James S. Wills (james.s.willis@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_SORT_PART_H
#define SWIFT_SORT_PART_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "inline.h"

/**
 * @brief Entry in a list of sorted indices.
 */
struct sort_entry {

  /*! Distance on the axis */
  float d;

  /*! Particle index */
  int i;
};

/* Orientation of the cell pairs */
static const double runner_shift[13][3] = {
    {5.773502691896258e-01, 5.773502691896258e-01, 5.773502691896258e-01},
    {7.071067811865475e-01, 7.071067811865475e-01, 0.0},
    {5.773502691896258e-01, 5.773502691896258e-01, -5.773502691896258e-01},
    {7.071067811865475e-01, 0.0, 7.071067811865475e-01},
    {1.0, 0.0, 0.0},
    {7.071067811865475e-01, 0.0, -7.071067811865475e-01},
    {5.773502691896258e-01, -5.773502691896258e-01, 5.773502691896258e-01},
    {7.071067811865475e-01, -7.071067811865475e-01, 0.0},
    {5.773502691896258e-01, -5.773502691896258e-01, -5.773502691896258e-01},
    {0.0, 7.071067811865475e-01, 7.071067811865475e-01},
    {0.0, 1.0, 0.0},
    {0.0, 7.071067811865475e-01, -7.071067811865475e-01},
    {0.0, 0.0, 1.0},
};

/* Does the axis need flipping ? */
static const char runner_flip[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Map shift vector to sortlist. */
static const int sortlistID[27] = {
    /* ( -1 , -1 , -1 ) */ 0,
    /* ( -1 , -1 ,  0 ) */ 1,
    /* ( -1 , -1 ,  1 ) */ 2,
    /* ( -1 ,  0 , -1 ) */ 3,
    /* ( -1 ,  0 ,  0 ) */ 4,
    /* ( -1 ,  0 ,  1 ) */ 5,
    /* ( -1 ,  1 , -1 ) */ 6,
    /* ( -1 ,  1 ,  0 ) */ 7,
    /* ( -1 ,  1 ,  1 ) */ 8,
    /* (  0 , -1 , -1 ) */ 9,
    /* (  0 , -1 ,  0 ) */ 10,
    /* (  0 , -1 ,  1 ) */ 11,
    /* (  0 ,  0 , -1 ) */ 12,
    /* (  0 ,  0 ,  0 ) */ 0,
    /* (  0 ,  0 ,  1 ) */ 12,
    /* (  0 ,  1 , -1 ) */ 11,
    /* (  0 ,  1 ,  0 ) */ 10,
    /* (  0 ,  1 ,  1 ) */ 9,
    /* (  1 , -1 , -1 ) */ 8,
    /* (  1 , -1 ,  0 ) */ 7,
    /* (  1 , -1 ,  1 ) */ 6,
    /* (  1 ,  0 , -1 ) */ 5,
    /* (  1 ,  0 ,  0 ) */ 4,
    /* (  1 ,  0 ,  1 ) */ 3,
    /* (  1 ,  1 , -1 ) */ 2,
    /* (  1 ,  1 ,  0 ) */ 1,
    /* (  1 ,  1 ,  1 ) */ 0};

/* Ratio of particles interacting assuming a uniform distribution */
static const float sid_scale[13] = {0.1897f, 0.4025f, 0.1897f, 0.4025f, 0.5788f,
                                    0.4025f, 0.1897f, 0.4025f, 0.1897f, 0.4025f,
                                    0.5788f, 0.4025f, 0.5788f};

/* Sid flags for every sub-pair of a self task. */
static const int sub_sid_flag[7][8] = {
    {-1, 12, 10, 9, 4, 3, 1, 0},     {-1, -1, 11, 10, 5, 4, 2, 1},
    {-1, -1, -1, 12, 7, 6, 4, 3},    {-1, -1, -1, -1, 8, 7, 5, 4},
    {-1, -1, -1, -1, -1, 12, 10, 9}, {-1, -1, -1, -1, -1, -1, 11, 10},
    {-1, -1, -1, -1, -1, -1, -1, 12}};

/**
 * @brief Determines whether a pair of cells are corner to corner.
 *
 * @param sid sort ID
 *
 * @return 1 if corner to corner, 0 otherwise.
 */
__attribute__((always_inline, const)) INLINE static int sort_is_corner(
    const int sid) {
  return (sid == 0 || sid == 2 || sid == 6 || sid == 8);
}

/**
 * @brief Determines whether a pair of cells are edge to edge.
 *
 * @param sid sort ID
 *
 * @return 1 if edge to edge, 0 otherwise.
 */
__attribute__((always_inline, const)) INLINE static int sort_is_edge(
    const int sid) {
  return (sid == 1 || sid == 3 || sid == 5 || sid == 7 || sid == 9 ||
          sid == 11);
}

/**
 * @brief Determines whether a pair of cells are face to face.
 *
 * @param sid sort ID
 *
 * @return 1 if face to face, 0 otherwise.
 */
__attribute__((always_inline, const)) INLINE static int sort_is_face(
    const int sid) {
  return (sid == 4 || sid == 10 || sid == 12);
}

/**
 * @brief Returns the position of the cell interface on the axis linking
 * two neighbouring cells ci and cj.
 *
 * @param sid The direction of interaction
 * @param cell_loc The location of cj.
 * @param cell_width The width of the cells.
 */
INLINE static double sort_get_cell_min_dist(const int sid,
                                            const double cell_loc[3],
                                            const double cell_width[3]) {

  double pos[3];

  switch (sid) {

    case 0: /* ( -1 , -1 , -1 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2];
      break;
    case 1: /* ( -1 , -1 ,  0 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1];
      pos[2] = 0.;
      break;
    case 2: /* ( -1 , -1 ,  1 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 3: /* ( -1 ,  0 , -1 ) */
      pos[0] = cell_loc[0];
      pos[1] = 0.;
      pos[2] = cell_loc[2];
      break;
    case 4: /* ( -1 ,  0 ,  0 ) */
      pos[0] = cell_loc[0];
      pos[1] = 0.;
      pos[2] = 0.;
      break;
    case 5: /* ( -1 ,  0 ,  1 ) */
      pos[0] = cell_loc[0];
      pos[1] = 0.;
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 6: /* ( -1 ,  1 , -1 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2];
      break;
    case 7: /* ( -1 ,  1 ,  0 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = 0.;
      break;
    case 8: /* ( -1 ,  1 ,  1 ) */
      pos[0] = cell_loc[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 9: /* ( 0 , -1 , -1 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2];
      break;
    case 10: /* ( 0 , -1 ,  0 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1];
      pos[2] = 0.;
      break;
    case 11: /* ( 0 , -1 ,  1 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 12: /* ( 0 ,  0 , -1 ) */
      pos[0] = 0.;
      pos[1] = 0.;
      pos[2] = cell_loc[2];
      break;
    case 13: /* ( 0 ,  0 ,  0 ) */
      pos[0] = 0.;
      pos[1] = 0.;
      pos[2] = 0.;
      break;
    case 14: /* ( 0 ,  0 ,  1 ) */
      pos[0] = 0.;
      pos[1] = 0.;
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 15: /* ( 0 ,  1 , -1 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2];
      break;
    case 16: /* ( 0 ,  1 ,  0 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = 0.;
      break;
    case 17: /* ( 0 ,  1 ,  1 ) */
      pos[0] = 0.;
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 18: /* ( 1 , -1 , -1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2];
      break;
    case 19: /* ( 1 , -1 ,  0 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1];
      pos[2] = 0.;
      break;
    case 20: /* ( 1 , -1 ,  1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 21: /* ( 1 ,  0 , -1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = 0.;
      pos[2] = cell_loc[2];
      break;
    case 22: /* ( 1 ,  0 ,  0 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = 0.;
      pos[2] = 0.;
      break;
    case 23: /* ( 1 ,  0 ,  1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = 0.;
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    case 24: /* ( 1 ,  1 , -1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2];
      break;
    case 25: /* ( 1 ,  1 ,  0 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = 0.;
      break;
    case 26: /* ( 1 ,  1 ,  1 ) */
      pos[0] = cell_loc[0] + cell_width[0];
      pos[1] = cell_loc[1] + cell_width[1];
      pos[2] = cell_loc[2] + cell_width[2];
      break;
    default:
      error("Invalid sid");
      pos[0] = 0.;
      pos[1] = 0.;
      pos[2] = 0.;
  }

  return pos[0] * runner_shift[sid][0] + pos[1] * runner_shift[sid][1] +
         pos[2] * runner_shift[sid][2];
}

#endif /* SWIFT_SORT_PART_H */
