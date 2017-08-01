/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "runner_doiact_vec.h"

/* Local headers. */
#include "active.h"

#ifdef WITH_VECTORIZATION
/**
 * @brief Compute the vector remainder interactions from the secondary cache.
 *
 * @param int_cache (return) secondary #cache of interactions between two
 * particles.
 * @param icount Interaction count.
 * @param rhoSum (return) #vector holding the cumulative sum of the density
 * update on pi.
 * @param rho_dhSum (return) #vector holding the cumulative sum of the density
 * gradient update on pi.
 * @param wcountSum (return) #vector holding the cumulative sum of the wcount
 * update on pi.
 * @param wcount_dhSum (return) #vector holding the cumulative sum of the wcount
 * gradient update on pi.
 * @param div_vSum (return) #vector holding the cumulative sum of the divergence
 * update on pi.
 * @param curlvxSum (return) #vector holding the cumulative sum of the curl of
 * vx update on pi.
 * @param curlvySum (return) #vector holding the cumulative sum of the curl of
 * vy update on pi.
 * @param curlvzSum (return) #vector holding the cumulative sum of the curl of
 * vz update on pi.
 * @param v_hi_inv #vector of 1/h for pi.
 * @param v_vix #vector of x velocity of pi.
 * @param v_viy #vector of y velocity of pi.
 * @param v_viz #vector of z velocity of pi.
 * @param icount_align (return) Interaction count after the remainder
 * interactions have been performed, should be a multiple of the vector length.
 */
__attribute__((always_inline)) INLINE static void calcRemInteractions(
    struct c2_cache *const int_cache, const int icount, vector *rhoSum,
    vector *rho_dhSum, vector *wcountSum, vector *wcount_dhSum,
    vector *div_vSum, vector *curlvxSum, vector *curlvySum, vector *curlvzSum,
    vector v_hi_inv, vector v_vix, vector v_viy, vector v_viz,
    int *icount_align) {

#ifdef HAVE_AVX512_F
  KNL_MASK_16 knl_mask, knl_mask2;
#endif
  vector int_mask, int_mask2;

  /* Work out the number of remainder interactions and pad secondary cache. */
  *icount_align = icount;
  int rem = icount % (NUM_VEC_PROC * VEC_SIZE);
  if (rem != 0) {
    int pad = (NUM_VEC_PROC * VEC_SIZE) - rem;
    *icount_align += pad;

/* Initialise masks to true. */
#ifdef HAVE_AVX512_F
    knl_mask = 0xFFFF;
    knl_mask2 = 0xFFFF;
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
#else
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
#endif
    /* Pad secondary cache so that there are no contributions in the interaction
     * function. */
    for (int i = icount; i < *icount_align; i++) {
      int_cache->mq[i] = 0.f;
      int_cache->r2q[i] = 1.f;
      int_cache->dxq[i] = 0.f;
      int_cache->dyq[i] = 0.f;
      int_cache->dzq[i] = 0.f;
      int_cache->vxq[i] = 0.f;
      int_cache->vyq[i] = 0.f;
      int_cache->vzq[i] = 0.f;
    }

    /* Zero parts of mask that represent the padded values.*/
    if (pad < VEC_SIZE) {
#ifdef HAVE_AVX512_F
      knl_mask2 = knl_mask2 >> pad;
#else
      for (int i = VEC_SIZE - pad; i < VEC_SIZE; i++) int_mask2.i[i] = 0;
#endif
    } else {
#ifdef HAVE_AVX512_F
      knl_mask = knl_mask >> (VEC_SIZE - rem);
      knl_mask2 = 0;
#else
      for (int i = rem; i < VEC_SIZE; i++) int_mask.i[i] = 0;
      int_mask2.v = vec_setzero();
#endif
    }

    /* Perform remainder interaction and remove remainder from aligned
     * interaction count. */
    *icount_align = icount - rem;
    runner_iact_nonsym_2_vec_density(
        &int_cache->r2q[*icount_align], &int_cache->dxq[*icount_align],
        &int_cache->dyq[*icount_align], &int_cache->dzq[*icount_align],
        v_hi_inv, v_vix, v_viy, v_viz, &int_cache->vxq[*icount_align],
        &int_cache->vyq[*icount_align], &int_cache->vzq[*icount_align],
        &int_cache->mq[*icount_align], rhoSum, rho_dhSum, wcountSum,
        wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum, int_mask,
        int_mask2,
#ifdef HAVE_AVX512_F
        knl_mask, knl_mask2);
#else
        0, 0);
#endif
  }
}

/**
 * @brief Left-packs the values needed by an interaction into the secondary
 * cache (Supports AVX, AVX2 and AVX512 instruction sets).
 *
 * @param mask Contains which particles need to interact.
 * @param pjd Index of the particle to store into.
 * @param v_r2 #vector of the separation between two particles squared.
 * @param v_dx #vector of the x separation between two particles.
 * @param v_dy #vector of the y separation between two particles.
 * @param v_dz #vector of the z separation between two particles.
 * @param v_mj #vector of the mass of particle pj.
 * @param v_vjx #vector of x velocity of pj.
 * @param v_vjy #vector of y velocity of pj.
 * @param v_vjz #vector of z velocity of pj.
 * @param cell_cache #cache of all particles in the cell.
 * @param int_cache (return) secondary #cache of interactions between two
 * particles.
 * @param icount Interaction count.
 * @param rhoSum #vector holding the cumulative sum of the density update on pi.
 * @param rho_dhSum #vector holding the cumulative sum of the density gradient
 * update on pi.
 * @param wcountSum #vector holding the cumulative sum of the wcount update on
 * pi.
 * @param wcount_dhSum #vector holding the cumulative sum of the wcount gradient
 * update on pi.
 * @param div_vSum #vector holding the cumulative sum of the divergence update
 * on pi.
 * @param curlvxSum #vector holding the cumulative sum of the curl of vx update
 * on pi.
 * @param curlvySum #vector holding the cumulative sum of the curl of vy update
 * on pi.
 * @param curlvzSum #vector holding the cumulative sum of the curl of vz update
 * on pi.
 * @param v_hi_inv #vector of 1/h for pi.
 * @param v_vix #vector of x velocity of pi.
 * @param v_viy #vector of y velocity of pi.
 * @param v_viz #vector of z velocity of pi.
 */
__attribute__((always_inline)) INLINE static void storeInteractions(
    const int mask, const int pjd, vector *v_r2, vector *v_dx, vector *v_dy,
    vector *v_dz, vector *v_mj, vector *v_vjx, vector *v_vjy, vector *v_vjz,
    const struct cache *const cell_cache, struct c2_cache *const int_cache,
    int *icount, vector *rhoSum, vector *rho_dhSum, vector *wcountSum,
    vector *wcount_dhSum, vector *div_vSum, vector *curlvxSum,
    vector *curlvySum, vector *curlvzSum, vector v_hi_inv, vector v_vix,
    vector v_viy, vector v_viz) {

/* Left-pack values needed into the secondary cache using the interaction mask.
 */
#if defined(HAVE_AVX2) || defined(HAVE_AVX512_F)
  int pack = 0;

#ifdef HAVE_AVX512_F
  pack += __builtin_popcount(mask);
  VEC_LEFT_PACK(v_r2->v, mask, &int_cache->r2q[*icount]);
  VEC_LEFT_PACK(v_dx->v, mask, &int_cache->dxq[*icount]);
  VEC_LEFT_PACK(v_dy->v, mask, &int_cache->dyq[*icount]);
  VEC_LEFT_PACK(v_dz->v, mask, &int_cache->dzq[*icount]);
  VEC_LEFT_PACK(v_mj->v, mask, &int_cache->mq[*icount]);
  VEC_LEFT_PACK(v_vjx->v, mask, &int_cache->vxq[*icount]);
  VEC_LEFT_PACK(v_vjy->v, mask, &int_cache->vyq[*icount]);
  VEC_LEFT_PACK(v_vjz->v, mask, &int_cache->vzq[*icount]);
#else
  vector v_mask;
  VEC_FORM_PACKED_MASK(mask, v_mask.m, pack);

  VEC_LEFT_PACK(v_r2->v, v_mask.m, &int_cache->r2q[*icount]);
  VEC_LEFT_PACK(v_dx->v, v_mask.m, &int_cache->dxq[*icount]);
  VEC_LEFT_PACK(v_dy->v, v_mask.m, &int_cache->dyq[*icount]);
  VEC_LEFT_PACK(v_dz->v, v_mask.m, &int_cache->dzq[*icount]);
  VEC_LEFT_PACK(v_mj->v, v_mask.m, &int_cache->mq[*icount]);
  VEC_LEFT_PACK(v_vjx->v, v_mask.m, &int_cache->vxq[*icount]);
  VEC_LEFT_PACK(v_vjy->v, v_mask.m, &int_cache->vyq[*icount]);
  VEC_LEFT_PACK(v_vjz->v, v_mask.m, &int_cache->vzq[*icount]);

#endif /* HAVE_AVX512_F */

  (*icount) += pack;
#else
  /* Quicker to do it serially in AVX rather than use intrinsics. */
  for (int bit_index = 0; bit_index < VEC_SIZE; bit_index++) {
    if (mask & (1 << bit_index)) {
      /* Add this interaction to the queue. */
      int_cache->r2q[*icount] = v_r2->f[bit_index];
      int_cache->dxq[*icount] = v_dx->f[bit_index];
      int_cache->dyq[*icount] = v_dy->f[bit_index];
      int_cache->dzq[*icount] = v_dz->f[bit_index];
      int_cache->mq[*icount] = cell_cache->m[pjd + bit_index];
      int_cache->vxq[*icount] = cell_cache->vx[pjd + bit_index];
      int_cache->vyq[*icount] = cell_cache->vy[pjd + bit_index];
      int_cache->vzq[*icount] = cell_cache->vz[pjd + bit_index];

      (*icount)++;
    }
  }

#endif /* defined(HAVE_AVX2) || defined(HAVE_AVX512_F) */

  /* Flush the c2 cache if it has reached capacity. */
  if (*icount >= (C2_CACHE_SIZE - (NUM_VEC_PROC * VEC_SIZE))) {

    int icount_align = *icount;

    /* Peform remainder interactions. */
    calcRemInteractions(int_cache, *icount, rhoSum, rho_dhSum, wcountSum,
                        wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum,
                        v_hi_inv, v_vix, v_viy, v_viz, &icount_align);

    vector int_mask, int_mask2;
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);

    /* Perform interactions. */
    for (int pjd = 0; pjd < icount_align; pjd += (NUM_VEC_PROC * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(
          &int_cache->r2q[pjd], &int_cache->dxq[pjd], &int_cache->dyq[pjd],
          &int_cache->dzq[pjd], v_hi_inv, v_vix, v_viy, v_viz,
          &int_cache->vxq[pjd], &int_cache->vyq[pjd], &int_cache->vzq[pjd],
          &int_cache->mq[pjd], rhoSum, rho_dhSum, wcountSum, wcount_dhSum,
          div_vSum, curlvxSum, curlvySum, curlvzSum, int_mask, int_mask2, 0, 0);
    }

    /* Reset interaction count. */
    *icount = 0;
  }
}

/**
 * @brief Populates the arrays max_di and max_dj with the maximum distances of
 * particles into their neighbouring cells. Also finds the first pi that
 * interacts with any particle in cj and the last pj that interacts with any
 * particle in ci.
 *
 * @param ci #cell pointer to ci
 * @param cj #cell pointer to cj
 * @param sort_i #entry array for particle distance in ci
 * @param sort_j #entry array for particle distance in cj
 * @param dx_max maximum particle movement allowed in cell
 * @param rshift cutoff shift
 * @param hi_max Maximal smoothing length in cell ci
 * @param hj_max Maximal smoothing length in cell cj
 * @param di_max Maximal position on the axis that can interact in cell ci
 * @param dj_min Minimal position on the axis that can interact in cell ci
 * @param max_di array to hold the maximum distances of pi particles into cell
 * cj
 * @param max_dj array to hold the maximum distances of pj particles into cell
 * cj
 * @param init_pi first pi to interact with a pj particle
 * @param init_pj last pj to interact with a pi particle
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void populate_max_d_no_cache(
    const struct cell *ci, const struct cell *cj,
    const struct entry *restrict sort_i, const struct entry *restrict sort_j,
    const float dx_max, const float rshift, const double hi_max,
    const double hj_max, const double di_max, const double dj_min,
    float *max_di, float *max_dj, int *init_pi, int *init_pj,
    const struct engine *e) {

  const struct part *restrict parts_i = ci->parts;
  const struct part *restrict parts_j = cj->parts;

  int first_pi = 0, last_pj = cj->count - 1;

  /* Find the first active particle in ci to interact with any particle in cj.
   */
  /* Populate max_di with distances. */
  int active_id = ci->count - 1;
  for (int k = ci->count - 1; k >= 0; k--) {
    const struct part *pi = &parts_i[sort_i[k].i];
    const float d = sort_i[k].d + dx_max;

    // max_di[k] = d + h * kernel_gamma - rshift;
    max_di[k] = d + hi_max;

    /* If the particle is out of range set the index to
     * the last active particle within range. */
    if (d + hi_max < dj_min) {
      first_pi = active_id;
      break;
    } else {
      if (part_is_active(pi, e)) active_id = k;
    }
  }

  /* Find the maximum distance of pi particles into cj.*/
  for (int k = first_pi + 1; k < ci->count; k++) {
    max_di[k] = fmaxf(max_di[k - 1], max_di[k]);
  }

  /* Find the last particle in cj to interact with any particle in ci. */
  /* Populate max_dj with distances. */
  active_id = 0;
  for (int k = 0; k < cj->count; k++) {
    const struct part *pj = &parts_j[sort_j[k].i];
    const float d = sort_j[k].d - dx_max;

    /*TODO: don't think rshift should be taken off here, waiting on Pedro. */
    // max_dj[k] = d - h * kernel_gamma - rshift;
    max_dj[k] = d - hj_max;

    /* If the particle is out of range set the index to
     * the last active particle within range. */
    if (d - hj_max > di_max) {
      last_pj = active_id;
      break;
    } else {
      if (part_is_active(pj, e)) active_id = k;
    }
  }

  /* Find the maximum distance of pj particles into ci.*/
  for (int k = 1; k <= last_pj; k++) {
    max_dj[k] = fmaxf(max_dj[k - 1], max_dj[k]);
  }

  *init_pi = first_pi;
  *init_pj = last_pj;
}
#endif /* WITH_VECTORIZATION */

/**
 * @brief Compute the cell self-interaction (non-symmetric) using vector
 * intrinsics with one particle pi at a time.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE void runner_doself1_density_vec(
    struct runner *r, struct cell *restrict c) {

#ifdef WITH_VECTORIZATION
  const struct engine *e = r->e;
  int doi_mask;
  struct part *restrict pi;
  int count_align;
  int num_vec_proc = NUM_VEC_PROC;

  struct part *restrict parts = c->parts;
  const int count = c->count;

  vector v_hi, v_vix, v_viy, v_viz, v_hig2, v_r2;

  TIMER_TIC

  if (!cell_is_active(c, e)) return;

  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  /* Get the particle cache from the runner and re-allocate
   * the cache if it is not big enough for the cell. */
  struct cache *restrict cell_cache = &r->ci_cache;

  if (cell_cache->count < count) {
    cache_init(cell_cache, count);
  }

  /* Read the particles from the cell and store them locally in the cache. */
  cache_read_particles(c, cell_cache);

  /* Create secondary cache to store particle interactions. */
  struct c2_cache int_cache;
  int icount = 0, icount_align = 0;

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    pi = &parts[pid];

    /* Is the ith particle active? */
    if (!part_is_active(pi, e)) continue;

    vector pix, piy, piz;

    const float hi = cell_cache->h[pid];

    /* Fill particle pi vectors. */
    pix.v = vec_set1(cell_cache->x[pid]);
    piy.v = vec_set1(cell_cache->y[pid]);
    piz.v = vec_set1(cell_cache->z[pid]);
    v_hi.v = vec_set1(hi);
    v_vix.v = vec_set1(cell_cache->vx[pid]);
    v_viy.v = vec_set1(cell_cache->vy[pid]);
    v_viz.v = vec_set1(cell_cache->vz[pid]);

    const float hig2 = hi * hi * kernel_gamma2;
    v_hig2.v = vec_set1(hig2);

    /* Reset cumulative sums of update vectors. */
    vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum,
        curlvySum, curlvzSum;

    /* Get the inverse of hi. */
    vector v_hi_inv;

    v_hi_inv = vec_reciprocal(v_hi);

    rhoSum.v = vec_setzero();
    rho_dhSum.v = vec_setzero();
    wcountSum.v = vec_setzero();
    wcount_dhSum.v = vec_setzero();
    div_vSum.v = vec_setzero();
    curlvxSum.v = vec_setzero();
    curlvySum.v = vec_setzero();
    curlvzSum.v = vec_setzero();

    /* Pad cache if there is a serial remainder. */
    count_align = count;
    int rem = count % (num_vec_proc * VEC_SIZE);
    if (rem != 0) {
      int pad = (num_vec_proc * VEC_SIZE) - rem;

      count_align += pad;

      /* Set positions to the same as particle pi so when the r2 > 0 mask is
       * applied these extra contributions are masked out.*/
      for (int i = count; i < count_align; i++) {
        cell_cache->x[i] = pix.f[0];
        cell_cache->y[i] = piy.f[0];
        cell_cache->z[i] = piz.f[0];
      }
    }

    vector pjx, pjy, pjz;
    vector pjvx, pjvy, pjvz, mj;
    vector pjx2, pjy2, pjz2;
    vector pjvx2, pjvy2, pjvz2, mj2;

    /* Find all of particle pi's interacions and store needed values in the
     * secondary cache.*/
    for (int pjd = 0; pjd < count_align; pjd += (num_vec_proc * VEC_SIZE)) {

      /* Load 2 sets of vectors from the particle cache. */
      pjx.v = vec_load(&cell_cache->x[pjd]);
      pjy.v = vec_load(&cell_cache->y[pjd]);
      pjz.v = vec_load(&cell_cache->z[pjd]);
      pjvx.v = vec_load(&cell_cache->vx[pjd]);
      pjvy.v = vec_load(&cell_cache->vy[pjd]);
      pjvz.v = vec_load(&cell_cache->vz[pjd]);
      mj.v = vec_load(&cell_cache->m[pjd]);

      pjx2.v = vec_load(&cell_cache->x[pjd + VEC_SIZE]);
      pjy2.v = vec_load(&cell_cache->y[pjd + VEC_SIZE]);
      pjz2.v = vec_load(&cell_cache->z[pjd + VEC_SIZE]);
      pjvx2.v = vec_load(&cell_cache->vx[pjd + VEC_SIZE]);
      pjvy2.v = vec_load(&cell_cache->vy[pjd + VEC_SIZE]);
      pjvz2.v = vec_load(&cell_cache->vz[pjd + VEC_SIZE]);
      mj2.v = vec_load(&cell_cache->m[pjd + VEC_SIZE]);

      /* Compute the pairwise distance. */
      vector v_dx_tmp, v_dy_tmp, v_dz_tmp;
      vector v_dx_tmp2, v_dy_tmp2, v_dz_tmp2, v_r2_2;

      v_dx_tmp.v = vec_sub(pix.v, pjx.v);
      v_dx_tmp2.v = vec_sub(pix.v, pjx2.v);
      v_dy_tmp.v = vec_sub(piy.v, pjy.v);
      v_dy_tmp2.v = vec_sub(piy.v, pjy2.v);
      v_dz_tmp.v = vec_sub(piz.v, pjz.v);
      v_dz_tmp2.v = vec_sub(piz.v, pjz2.v);

      v_r2.v = vec_mul(v_dx_tmp.v, v_dx_tmp.v);
      v_r2_2.v = vec_mul(v_dx_tmp2.v, v_dx_tmp2.v);
      v_r2.v = vec_fma(v_dy_tmp.v, v_dy_tmp.v, v_r2.v);
      v_r2_2.v = vec_fma(v_dy_tmp2.v, v_dy_tmp2.v, v_r2_2.v);
      v_r2.v = vec_fma(v_dz_tmp.v, v_dz_tmp.v, v_r2.v);
      v_r2_2.v = vec_fma(v_dz_tmp2.v, v_dz_tmp2.v, v_r2_2.v);

/* Form a mask from r2 < hig2 and r2 > 0.*/
#ifdef HAVE_AVX512_F
      // KNL_MASK_16 doi_mask, doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi_mask_check, doi_mask2, doi_mask2_check;

      doi_mask_check = vec_cmp_gt(v_r2.v, vec_setzero());
      doi_mask = vec_cmp_lt(v_r2.v, v_hig2.v);

      doi_mask2_check = vec_cmp_gt(v_r2_2.v, vec_setzero());
      doi_mask2 = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      doi_mask = doi_mask & doi_mask_check;
      doi_mask2 = doi_mask2 & doi_mask2_check;

#else
      vector v_doi_mask, v_doi_mask_check, v_doi_mask2, v_doi_mask2_check;
      int doi_mask2;

      /* Form r2 > 0 mask and r2 < hig2 mask. */
      v_doi_mask_check.v = vec_cmp_gt(v_r2.v, vec_setzero());
      v_doi_mask.v = vec_cmp_lt(v_r2.v, v_hig2.v);

      /* Form r2 > 0 mask and r2 < hig2 mask. */
      v_doi_mask2_check.v = vec_cmp_gt(v_r2_2.v, vec_setzero());
      v_doi_mask2.v = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      /* Combine two masks and form integer mask. */
      doi_mask = vec_cmp_result(vec_and(v_doi_mask.v, v_doi_mask_check.v));
      doi_mask2 = vec_cmp_result(vec_and(v_doi_mask2.v, v_doi_mask2_check.v));
#endif /* HAVE_AVX512_F */

      /* If there are any interactions left pack interaction values into c2
       * cache. */
      if (doi_mask) {
        storeInteractions(doi_mask, pjd, &v_r2, &v_dx_tmp, &v_dy_tmp, &v_dz_tmp,
                          &mj, &pjvx, &pjvy, &pjvz, cell_cache, &int_cache,
                          &icount, &rhoSum, &rho_dhSum, &wcountSum,
                          &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum,
                          &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      }
      if (doi_mask2) {
        storeInteractions(
            doi_mask2, pjd + VEC_SIZE, &v_r2_2, &v_dx_tmp2, &v_dy_tmp2,
            &v_dz_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2, cell_cache, &int_cache,
            &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum,
            &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      }
    }

    /* Perform padded vector remainder interactions if any are present. */
    calcRemInteractions(&int_cache, icount, &rhoSum, &rho_dhSum, &wcountSum,
                        &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum,
                        &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz,
                        &icount_align);

    /* Initialise masks to true in case remainder interactions have been
     * performed. */
    vector int_mask, int_mask2;
#ifdef HAVE_AVX512_F
    KNL_MASK_16 knl_mask = 0xFFFF;
    KNL_MASK_16 knl_mask2 = 0xFFFF;
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
#else
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
#endif

    /* Perform interaction with 2 vectors. */
    for (int pjd = 0; pjd < icount_align; pjd += (num_vec_proc * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(
          &int_cache.r2q[pjd], &int_cache.dxq[pjd], &int_cache.dyq[pjd],
          &int_cache.dzq[pjd], v_hi_inv, v_vix, v_viy, v_viz,
          &int_cache.vxq[pjd], &int_cache.vyq[pjd], &int_cache.vzq[pjd],
          &int_cache.mq[pjd], &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum,
          &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, int_mask, int_mask2,
#ifdef HAVE_AVX512_F
          knl_mask, knl_mask2);
#else
          0, 0);
#endif
    }

    /* Perform horizontal adds on vector sums and store result in particle pi.
     */
    VEC_HADD(rhoSum, pi->rho);
    VEC_HADD(rho_dhSum, pi->density.rho_dh);
    VEC_HADD(wcountSum, pi->density.wcount);
    VEC_HADD(wcount_dhSum, pi->density.wcount_dh);
    VEC_HADD(div_vSum, pi->density.div_v);
    VEC_HADD(curlvxSum, pi->density.rot_v[0]);
    VEC_HADD(curlvySum, pi->density.rot_v[1]);
    VEC_HADD(curlvzSum, pi->density.rot_v[2]);

    /* Reset interaction count. */
    icount = 0;
  } /* loop over all particles. */

  TIMER_TOC(timer_doself_density);
#endif /* WITH_VECTORIZATION */
}

/**
 * @brief Compute the density interactions between a cell pair (non-symmetric)
 * using vector intrinsics.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair
 * @param shift The shift vector to apply to the particles in ci.
 */
void runner_dopair1_density_vec(struct runner *r, struct cell *ci,
                                struct cell *cj, const int sid,
                                const double *shift) {

#ifdef WITH_VECTORIZATION
  const struct engine *restrict e = r->e;

  vector v_hi, v_vix, v_viy, v_viz, v_hig2;

  TIMER_TIC;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->sort[sid];
  const struct entry *restrict sort_j = cj->sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->count; pid++) {
    const struct part *p = &ci->parts[sort_i[pid].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->dx_max_sort >
        1.0e-4 * max(fabsf(d), ci->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->dx_max_sort=%e "
          "ci->dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->dx_max_sort,
          ci->dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->count; pjd++) {
    const struct part *p = &cj->parts[sort_j[pjd].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_j[pjd].d) - cj->dx_max_sort >
        1.0e-4 * max(fabsf(d), cj->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->dx_max_sort=%e "
          "cj->dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->dx_max_sort,
          cj->dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get some other useful values. */
  const int count_i = ci->count;
  const int count_j = cj->count;
  const double hi_max = ci->h_max * kernel_gamma - rshift;
  const double hj_max = cj->h_max * kernel_gamma;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->dx_max_sort + cj->dx_max_sort);

  /* Check if any particles are active and return if there are not. */
  int numActive = 0;
  for (int pid = count_i - 1;
       pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {
    struct part *restrict pi = &parts_i[sort_i[pid].i];
    if (part_is_active(pi, e)) {
      numActive++;
      break;
    }
  }

  if (!numActive) {
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {
      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      if (part_is_active(pj, e)) {
        numActive++;
        break;
      }
    }
  }

  if (numActive == 0) return;

  /* Get both particle caches from the runner and re-allocate
   * them if they are not big enough for the cells. */
  struct cache *restrict ci_cache = &r->ci_cache;
  struct cache *restrict cj_cache = &r->cj_cache;

  if (ci_cache->count < count_i) {
    cache_init(ci_cache, count_i);
  }
  if (cj_cache->count < count_j) {
    cache_init(cj_cache, count_j);
  }

  int first_pi, last_pj;
  float *max_di __attribute__((aligned(sizeof(float) * VEC_SIZE)));
  float *max_dj __attribute__((aligned(sizeof(float) * VEC_SIZE)));

  max_di = r->ci_cache.max_d;
  max_dj = r->cj_cache.max_d;

  /* Find particles maximum distance into cj, max_di[] and ci, max_dj[]. */
  /* Also find the first pi that interacts with any particle in cj and the last
   * pj that interacts with any particle in ci. */
  populate_max_d_no_cache(ci, cj, sort_i, sort_j, dx_max, rshift, hi_max,
                          hj_max, di_max, dj_min, max_di, max_dj, &first_pi,
                          &last_pj, e);

  /* Find the maximum index into cj that is required by a particle in ci. */
  /* Find the maximum index into ci that is required by a particle in cj. */
  float di, dj;
  int max_ind_j = count_j - 1;
  int max_ind_i = 0;

  dj = sort_j[max_ind_j].d;
  while (max_ind_j > 0 && max_di[count_i - 1] < dj) {
    max_ind_j--;

    dj = sort_j[max_ind_j].d;
  }

  di = sort_i[max_ind_i].d;
  while (max_ind_i < count_i - 1 && max_dj[0] > di) {
    max_ind_i++;

    di = sort_i[max_ind_i].d;
  }

  /* Limits of the outer loops. */
  int first_pi_loop = first_pi;
  int last_pj_loop = last_pj;

  /* Take the max/min of both values calculated to work out how many particles
   * to read into the cache. */
  last_pj = max(last_pj, max_ind_j);
  first_pi = min(first_pi, max_ind_i);

  /* Read the needed particles into the two caches. */
  int first_pi_align = first_pi;
  int last_pj_align = last_pj;
  cache_read_two_partial_cells_sorted(ci, cj, ci_cache, cj_cache, sort_i,
                                      sort_j, shift, &first_pi_align,
                                      &last_pj_align, 1);

  /* Get the number of particles read into the ci cache. */
  int ci_cache_count = count_i - first_pi_align;

  if (cell_is_active(ci, e)) {

    /* Loop over the parts in ci until nothing is within range in cj. */
    for (int pid = count_i - 1; pid >= first_pi_loop && max_ind_j >= 0; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      if (!part_is_active(pi, e)) continue;

      /* Set the cache index. */
      int ci_cache_idx = pid - first_pi_align;

      /* Skip this particle if no particle in cj is within range of it. */
      const float hi = ci_cache->h[ci_cache_idx];
      const double di_test =
          sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di_test < dj_min) continue;

      /* Determine the exit iteration of the interaction loop. */
      dj = sort_j[max_ind_j].d;
      while (max_ind_j > 0 && max_di[pid] < dj) {
        max_ind_j--;

        dj = sort_j[max_ind_j].d;
      }
      int exit_iteration = max_ind_j + 1;

      const float hig2 = hi * hi * kernel_gamma2;

      vector pix, piy, piz;

      /* Fill particle pi vectors. */
      pix.v = vec_set1(ci_cache->x[ci_cache_idx]);
      piy.v = vec_set1(ci_cache->y[ci_cache_idx]);
      piz.v = vec_set1(ci_cache->z[ci_cache_idx]);
      v_hi.v = vec_set1(hi);
      v_vix.v = vec_set1(ci_cache->vx[ci_cache_idx]);
      v_viy.v = vec_set1(ci_cache->vy[ci_cache_idx]);
      v_viz.v = vec_set1(ci_cache->vz[ci_cache_idx]);

      v_hig2.v = vec_set1(hig2);

      /* Reset cumulative sums of update vectors. */
      vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum,
          curlvySum, curlvzSum;

      /* Get the inverse of hi. */
      vector v_hi_inv;

      v_hi_inv = vec_reciprocal(v_hi);

      rhoSum.v = vec_setzero();
      rho_dhSum.v = vec_setzero();
      wcountSum.v = vec_setzero();
      wcount_dhSum.v = vec_setzero();
      div_vSum.v = vec_setzero();
      curlvxSum.v = vec_setzero();
      curlvySum.v = vec_setzero();
      curlvzSum.v = vec_setzero();

      /* Pad the exit iteration if there is a serial remainder. */
      int exit_iteration_align = exit_iteration;
      int rem = exit_iteration % VEC_SIZE;
      if (rem != 0) {
        int pad = VEC_SIZE - rem;

        if (exit_iteration_align + pad <= last_pj_align + 1)
          exit_iteration_align += pad;
      }

      vector pjx, pjy, pjz;

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < exit_iteration_align; pjd += VEC_SIZE) {

        /* Get the cache index to the jth particle. */
        int cj_cache_idx = pjd;

        vector v_dx, v_dy, v_dz, v_r2;

#ifdef SWIFT_DEBUG_CHECKS
        if (cj_cache_idx % VEC_SIZE != 0 || cj_cache_idx < 0) {
          error("Unaligned read!!! cj_cache_idx=%d", cj_cache_idx);
        }
#endif

        /* Load 2 sets of vectors from the particle cache. */
        pjx.v = vec_load(&cj_cache->x[cj_cache_idx]);
        pjy.v = vec_load(&cj_cache->y[cj_cache_idx]);
        pjz.v = vec_load(&cj_cache->z[cj_cache_idx]);

        /* Compute the pairwise distance. */
        v_dx.v = vec_sub(pix.v, pjx.v);
        v_dy.v = vec_sub(piy.v, pjy.v);
        v_dz.v = vec_sub(piz.v, pjz.v);

        v_r2.v = vec_mul(v_dx.v, v_dx.v);
        v_r2.v = vec_fma(v_dy.v, v_dy.v, v_r2.v);
        v_r2.v = vec_fma(v_dz.v, v_dz.v, v_r2.v);

        vector v_doi_mask;
        int doi_mask;

        /* Form r2 < hig2 mask. */
        v_doi_mask.v = vec_cmp_lt(v_r2.v, v_hig2.v);

        /* Form integer mask. */
        doi_mask = vec_cmp_result(v_doi_mask.v);

        /* If there are any interactions perform them. */
        if (doi_mask)
          runner_iact_nonsym_1_vec_density(
              &v_r2, &v_dx, &v_dy, &v_dz, v_hi_inv, v_vix, v_viy, v_viz,
              &cj_cache->vx[cj_cache_idx], &cj_cache->vy[cj_cache_idx],
              &cj_cache->vz[cj_cache_idx], &cj_cache->m[cj_cache_idx], &rhoSum,
              &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum,
              &curlvySum, &curlvzSum, v_doi_mask,
#ifdef HAVE_AVX512_F
              knl_mask);
#else
              0);
#endif

      } /* loop over the parts in cj. */

      /* Perform horizontal adds on vector sums and store result in particle pi.
      */
      VEC_HADD(rhoSum, pi->rho);
      VEC_HADD(rho_dhSum, pi->density.rho_dh);
      VEC_HADD(wcountSum, pi->density.wcount);
      VEC_HADD(wcount_dhSum, pi->density.wcount_dh);
      VEC_HADD(div_vSum, pi->density.div_v);
      VEC_HADD(curlvxSum, pi->density.rot_v[0]);
      VEC_HADD(curlvySum, pi->density.rot_v[1]);
      VEC_HADD(curlvzSum, pi->density.rot_v[2]);

    } /* loop over the parts in ci. */
  }

  if (cell_is_active(cj, e)) {

    /* Loop over the parts in cj until nothing is within range in ci. */
    for (int pjd = 0; pjd <= last_pj_loop && max_ind_i < count_i; pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      if (!part_is_active(pj, e)) continue;

      /* Set the cache index. */
      int cj_cache_idx = pjd;

      /*TODO: rshift term. */
      /* Skip this particle if no particle in ci is within range of it. */
      const float hj = cj_cache->h[cj_cache_idx];
      const double dj_test =
          sort_j[pjd].d - hj * kernel_gamma - dx_max - rshift;
      if (dj_test > di_max) continue;

      /* Determine the exit iteration of the interaction loop. */
      di = sort_i[max_ind_i].d;
      while (max_ind_i < count_i - 1 && max_dj[pjd] > di) {
        max_ind_i++;

        di = sort_i[max_ind_i].d;
      }
      int exit_iteration = max_ind_i;

      const float hjg2 = hj * hj * kernel_gamma2;

      vector pjx, pjy, pjz;
      vector v_hj, v_vjx, v_vjy, v_vjz, v_hjg2;

      /* Fill particle pi vectors. */
      pjx.v = vec_set1(cj_cache->x[cj_cache_idx]);
      pjy.v = vec_set1(cj_cache->y[cj_cache_idx]);
      pjz.v = vec_set1(cj_cache->z[cj_cache_idx]);
      v_hj.v = vec_set1(hj);
      v_vjx.v = vec_set1(cj_cache->vx[cj_cache_idx]);
      v_vjy.v = vec_set1(cj_cache->vy[cj_cache_idx]);
      v_vjz.v = vec_set1(cj_cache->vz[cj_cache_idx]);

      v_hjg2.v = vec_set1(hjg2);

      /* Reset cumulative sums of update vectors. */
      vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum,
          curlvySum, curlvzSum;

      /* Get the inverse of hj. */
      vector v_hj_inv;

      v_hj_inv = vec_reciprocal(v_hj);

      rhoSum.v = vec_setzero();
      rho_dhSum.v = vec_setzero();
      wcountSum.v = vec_setzero();
      wcount_dhSum.v = vec_setzero();
      div_vSum.v = vec_setzero();
      curlvxSum.v = vec_setzero();
      curlvySum.v = vec_setzero();
      curlvzSum.v = vec_setzero();

      vector pix, piy, piz;

      /* Convert exit iteration to cache indices. */
      int exit_iteration_align = exit_iteration - first_pi_align;

      /* Pad the exit iteration align so cache reads are aligned. */
      int rem = exit_iteration_align % VEC_SIZE;
      if (exit_iteration_align < VEC_SIZE) {
        exit_iteration_align = 0;
      } else
        exit_iteration_align -= rem;

      /* Loop over the parts in ci. */
      for (int ci_cache_idx = exit_iteration_align;
           ci_cache_idx < ci_cache_count; ci_cache_idx += VEC_SIZE) {

#ifdef SWIFT_DEBUG_CHECKS
        if (ci_cache_idx % VEC_SIZE != 0 || ci_cache_idx < 0) {
          error("Unaligned read!!! ci_cache_idx=%d", ci_cache_idx);
        }
#endif

        vector v_dx, v_dy, v_dz, v_r2;

        /* Load 2 sets of vectors from the particle cache. */
        pix.v = vec_load(&ci_cache->x[ci_cache_idx]);
        piy.v = vec_load(&ci_cache->y[ci_cache_idx]);
        piz.v = vec_load(&ci_cache->z[ci_cache_idx]);

        /* Compute the pairwise distance. */
        v_dx.v = vec_sub(pjx.v, pix.v);
        v_dy.v = vec_sub(pjy.v, piy.v);
        v_dz.v = vec_sub(pjz.v, piz.v);

        v_r2.v = vec_mul(v_dx.v, v_dx.v);
        v_r2.v = vec_fma(v_dy.v, v_dy.v, v_r2.v);
        v_r2.v = vec_fma(v_dz.v, v_dz.v, v_r2.v);

        vector v_doj_mask;
        int doj_mask;

        /* Form r2 < hig2 mask. */
        v_doj_mask.v = vec_cmp_lt(v_r2.v, v_hjg2.v);

        /* Form integer mask. */
        doj_mask = vec_cmp_result(v_doj_mask.v);

        /* If there are any interactions perform them. */
        if (doj_mask)
          runner_iact_nonsym_1_vec_density(
              &v_r2, &v_dx, &v_dy, &v_dz, v_hj_inv, v_vjx, v_vjy, v_vjz,
              &ci_cache->vx[ci_cache_idx], &ci_cache->vy[ci_cache_idx],
              &ci_cache->vz[ci_cache_idx], &ci_cache->m[ci_cache_idx], &rhoSum,
              &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum,
              &curlvySum, &curlvzSum, v_doj_mask,
#ifdef HAVE_AVX512_F
              knl_mask);
#else
              0);
#endif

      } /* loop over the parts in ci. */

      /* Perform horizontal adds on vector sums and store result in particle pj.
       */
      VEC_HADD(rhoSum, pj->rho);
      VEC_HADD(rho_dhSum, pj->density.rho_dh);
      VEC_HADD(wcountSum, pj->density.wcount);
      VEC_HADD(wcount_dhSum, pj->density.wcount_dh);
      VEC_HADD(div_vSum, pj->density.div_v);
      VEC_HADD(curlvxSum, pj->density.rot_v[0]);
      VEC_HADD(curlvySum, pj->density.rot_v[1]);
      VEC_HADD(curlvzSum, pj->density.rot_v[2]);

    } /* loop over the parts in cj. */

    TIMER_TOC(timer_dopair_density);
  }

#endif /* WITH_VECTORIZATION */
}
