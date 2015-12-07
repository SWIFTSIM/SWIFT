/*******************************************************************************
 * This file is part of SWIFT.
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
#ifndef SWIFT_PARTITION_H
#define SWIFT_PARTITION_H

#include "space.h"
#include "cell.h"

int part_pick_random(struct space *s, int nregions, float *samplelist);
void part_split_random(struct space *s, int nregions, float *samplelist);

void part_pick_vector(struct space *s, int nregions, int *samplecells);
void part_split_vector(struct space *s, int nregions, int *samplecells);

void part_pick_metis(struct space *s, int nregions, int *weight, int *celllist);
void part_split_metis(struct space *s, int nregions, int *celllist);

#endif /* SWIFT_POISSON_DISC_H */
