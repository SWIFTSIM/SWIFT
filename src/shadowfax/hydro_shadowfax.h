#ifndef SWIFT_HYDRO_SHADOWFAX_H
#define SWIFT_HYDRO_SHADOWFAX_H

__attribute__((always_inline)) INLINE static void hydro_shadowfax_flux_exchange(
    struct part *restrict pi, struct part *restrict pj, double *midpoint, double surface_area) {
    
    ++pi->voronoi.nface;
    ++pj->voronoi.nface;
}

#endif /* SWIFT_HYDRO_SHADOWFAX_H */
