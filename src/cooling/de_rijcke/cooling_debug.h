//
// Created by yuyttenh on 10/01/23.
//

#ifndef SWIFTSIM_COOLING_DEBUG_DE_RIJCKE_H
#define SWIFTSIM_COOLING_DEBUG_DE_RIJCKE_H

__attribute__((always_inline)) INLINE static void cooling_debug_particle(
    const struct part* p, const struct xpart* xp) {

  if (xp != NULL) {
    warning("[PID%lld] cooling_xpart_data:", p->id);
    warning("[PID%lld] radiated_energy = %.3e", p->id,
            xp->cooling_data.radiated_energy);
  }
}

#endif  // SWIFTSIM_COOLING_DEBUG_DE_RIJCKE_H
