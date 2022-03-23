//
// Created by yuyttenh on 23/03/22.
//

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_H

void runner_dopair_grid_construction(struct runner *restrict r,
                                     struct cell *restrict ci,
                                     struct cell *restrict cj) {
  /* TODO */
}

void runner_doself_grid_construction(struct runner *restrict r,
                                     struct cell *restrict ci) {
  /* TODO */
}

void runner_dopair_subset_grid_construction(struct runner *restrict r,
                                            struct cell *restrict ci,
                                            struct part *restrict parts_i,
                                            int *restrict ind, int count,
                                            struct cell *restrict cj) {
  /* TODO */
}

void runner_doself_subset_grid_construction(struct runner *restrict r,
                                            struct cell *restrict ci,
                                            struct part *restrict parts_i,
                                            int *restrict ind, int count) {
  /* TODO */
}

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_H
