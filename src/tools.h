/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

#ifndef SWIFT_TOOL_H
#define SWIFT_TOOL_H

#include "cell.h"
#include "gravity_properties.h"
#include "part.h"
#include "physical_constants.h"
#include "runner.h"

void factor(int value, int *f1, int *f2);
void density_dump(int N);
void pairs_single_grav(double *dim, long long int pid,
                       struct gpart *restrict gparts, const struct part *parts,
                       int N, int periodic);
void pairs_single_density(double *dim, long long int pid,
                          struct part *restrict parts, int N, int periodic);

void pairs_all_density(struct runner *r, struct cell *ci, struct cell *cj);
void self_all_density(struct runner *r, struct cell *ci);
void pairs_all_gradient(struct runner *r, struct cell *ci, struct cell *cj);
void self_all_gradient(struct runner *r, struct cell *ci);
void pairs_all_force(struct runner *r, struct cell *ci, struct cell *cj);
void self_all_force(struct runner *r, struct cell *ci);
void pairs_all_stars_density(struct runner *r, struct cell *ci,
                             struct cell *cj);
void self_all_stars_density(struct runner *r, struct cell *ci);

void pairs_n2(double *dim, struct part *restrict parts, int N, int periodic);

double random_uniform(double a, double b);
void shuffle_particles(struct part *parts, const int count);
void shuffle_sparticles(struct spart *sparts, const int scount);
void gravity_n2(struct gpart *gparts, const int gcount,
                const struct phys_const *constants,
                const struct gravity_props *gravity_properties, float rlr);
int compare_values(double a, double b, double threshold, double *absDiff,
                   double *absSum, double *relDiff);
int compare_particles(struct part *a, struct part *b, double threshold);

long get_maxrss(void);

char *trim_leading(char *s);
char *trim_trailing(char *s);
char *trim_both(char *s);

#endif /* SWIFT_TOOL_H */
