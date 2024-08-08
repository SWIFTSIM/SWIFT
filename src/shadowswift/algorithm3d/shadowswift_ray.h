//
// Created by yuyttenh on 10/05/22.
//

#ifndef SWIFTSIM_GITLAB_SHADOWSWIFT_RAY_H
#define SWIFTSIM_GITLAB_SHADOWSWIFT_RAY_H

struct shadowswift_ray {
  double origin[3];
  double direction[3];
  unsigned long origin_ul[3];
  unsigned long end_ul[3];
};

inline static void shadowswift_ray_init(struct shadowswift_ray* r,
                                        const double* origin, const double* end,
                                        const unsigned long* origin_ul,
                                        const unsigned long* end_ul) {
  for (int i = 0; i < 3; i++) {
    r->origin[i] = origin[i];
    r->direction[i] = end[i] - origin[i];
    r->origin_ul[i] = origin_ul[i];
    r->end_ul[i] = end_ul[i];
  }
  double norm = sqrt(r->direction[0] * r->direction[0] +
                     r->direction[1] * r->direction[1] +
                     r->direction[2] * r->direction[2]);
  r->direction[0] /= norm;
  r->direction[1] /= norm;
  r->direction[2] /= norm;
}

#endif  // SWIFTSIM_GITLAB_SHADOWSWIFT_RAY_H
