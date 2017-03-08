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

#include "active.h"

/* This object's header. */
#include "runner_doiact_vec.h"

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

  if (!cell_is_drifted(c, e)) cell_drift_particles(c, e);

  /* Get the particle cache from the runner and re-allocate
   * the cache if it is not big enough for the cell. */
  struct cache *restrict cell_cache = &r->par_cache;

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
      v_dy_tmp.v = vec_sub(piy.v, pjy.v);
      v_dz_tmp.v = vec_sub(piz.v, pjz.v);
      v_dx_tmp2.v = vec_sub(pix.v, pjx2.v);
      v_dy_tmp2.v = vec_sub(piy.v, pjy2.v);
      v_dz_tmp2.v = vec_sub(piz.v, pjz2.v);

      v_r2.v = vec_mul(v_dx_tmp.v, v_dx_tmp.v);
      v_r2.v = vec_fma(v_dy_tmp.v, v_dy_tmp.v, v_r2.v);
      v_r2.v = vec_fma(v_dz_tmp.v, v_dz_tmp.v, v_r2.v);
      v_r2_2.v = vec_mul(v_dx_tmp2.v, v_dx_tmp2.v);
      v_r2_2.v = vec_fma(v_dy_tmp2.v, v_dy_tmp2.v, v_r2_2.v);
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
 * @brief Compute the cell self-interaction (non-symmetric) using vector
 * intrinsics with two particle pis at a time.
 *
 * CURRENTLY BROKEN DO NOT USE.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE void runner_doself1_density_vec_2(
    struct runner *r, struct cell *restrict c) {

#ifdef WITH_VECTORIZATION
  const struct engine *e = r->e;
  int doi_mask;
  int doi2_mask;
  struct part *restrict pi;
  struct part *restrict pi2;
  int count_align;

  vector v_hi, v_vix, v_viy, v_viz, v_hig2, v_r2;
  vector v_hi2, v_vix2, v_viy2, v_viz2, v_hig2_2, v2_r2;

  TIMER_TIC

  if (!cell_is_active(c, e)) return;

  if (!cell_is_drifted(c, e)) cell_drift_particles(c, e);

  /* TODO: Need to find two active particles, not just one. */

  struct part *restrict parts = c->parts;
  const int count = c->count;

  /* Get the particle cache from the runner and re-allocate
   * the cache if it is not big enough for the cell. */
  struct cache *restrict cell_cache = &r->par_cache;

  if (cell_cache->count < count) {
    cache_init(cell_cache, count);
  }

  /* Read the particles from the cell and store them locally in the cache. */
  cache_read_particles(c, &r->par_cache);

  /* Create two secondary caches. */
  int icount = 0, icount_align = 0;
  struct c2_cache int_cache;

  int icount2 = 0, icount_align2 = 0;
  struct c2_cache int_cache2;

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid += 2) {

    /* Get a pointer to the ith particle and next i particle. */
    pi = &parts[pid];
    pi2 = &parts[pid + 1];

    /* Is the ith particle active? */
    if (!part_is_active(pi, e)) continue;

    vector pix, piy, piz;
    vector pix2, piy2, piz2;

    const float hi = cell_cache->h[pid];
    const float hi2 = cell_cache->h[pid + 1];

    /* Fill pi position vector. */
    pix.v = vec_set1(cell_cache->x[pid]);
    piy.v = vec_set1(cell_cache->y[pid]);
    piz.v = vec_set1(cell_cache->z[pid]);
    v_hi.v = vec_set1(hi);
    v_vix.v = vec_set1(cell_cache->vx[pid]);
    v_viy.v = vec_set1(cell_cache->vy[pid]);
    v_viz.v = vec_set1(cell_cache->vz[pid]);

    pix2.v = vec_set1(cell_cache->x[pid + 1]);
    piy2.v = vec_set1(cell_cache->y[pid + 1]);
    piz2.v = vec_set1(cell_cache->z[pid + 1]);
    v_hi2.v = vec_set1(hi2);
    v_vix2.v = vec_set1(cell_cache->vx[pid + 1]);
    v_viy2.v = vec_set1(cell_cache->vy[pid + 1]);
    v_viz2.v = vec_set1(cell_cache->vz[pid + 1]);

    const float hig2 = hi * hi * kernel_gamma2;
    const float hig2_2 = hi2 * hi2 * kernel_gamma2;
    v_hig2.v = vec_set1(hig2);
    v_hig2_2.v = vec_set1(hig2_2);

    vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum,
        curlvySum, curlvzSum;
    vector rhoSum2, rho_dhSum2, wcountSum2, wcount_dhSum2, div_vSum2,
        curlvxSum2, curlvySum2, curlvzSum2;

    vector v_hi_inv, v_hi_inv2;

    v_hi_inv = vec_reciprocal(v_hi);
    v_hi_inv2 = vec_reciprocal(v_hi2);

    rhoSum.v = vec_setzero();
    rho_dhSum.v = vec_setzero();
    wcountSum.v = vec_setzero();
    wcount_dhSum.v = vec_setzero();
    div_vSum.v = vec_setzero();
    curlvxSum.v = vec_setzero();
    curlvySum.v = vec_setzero();
    curlvzSum.v = vec_setzero();

    rhoSum2.v = vec_setzero();
    rho_dhSum2.v = vec_setzero();
    wcountSum2.v = vec_setzero();
    wcount_dhSum2.v = vec_setzero();
    div_vSum2.v = vec_setzero();
    curlvxSum2.v = vec_setzero();
    curlvySum2.v = vec_setzero();
    curlvzSum2.v = vec_setzero();

    /* Pad cache if there is a serial remainder. */
    count_align = count;
    int rem = count % (NUM_VEC_PROC * VEC_SIZE);
    if (rem != 0) {
      int pad = (NUM_VEC_PROC * VEC_SIZE) - rem;

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

    /* Find all of particle pi's interacions and store needed values in
     * secondary cache.*/
    for (int pjd = 0; pjd < count_align; pjd += (NUM_VEC_PROC * VEC_SIZE)) {

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
      vector v_dx2_tmp, v_dy2_tmp, v_dz2_tmp;
      vector v_dx2_tmp2, v_dy2_tmp2, v_dz2_tmp2, v2_r2_2;

      v_dx_tmp.v = vec_sub(pix.v, pjx.v);
      v_dy_tmp.v = vec_sub(piy.v, pjy.v);
      v_dz_tmp.v = vec_sub(piz.v, pjz.v);
      v_dx_tmp2.v = vec_sub(pix.v, pjx2.v);
      v_dy_tmp2.v = vec_sub(piy.v, pjy2.v);
      v_dz_tmp2.v = vec_sub(piz.v, pjz2.v);

      v_dx2_tmp.v = vec_sub(pix2.v, pjx.v);
      v_dy2_tmp.v = vec_sub(piy2.v, pjy.v);
      v_dz2_tmp.v = vec_sub(piz2.v, pjz.v);
      v_dx2_tmp2.v = vec_sub(pix2.v, pjx2.v);
      v_dy2_tmp2.v = vec_sub(piy2.v, pjy2.v);
      v_dz2_tmp2.v = vec_sub(piz2.v, pjz2.v);

      v_r2.v = vec_mul(v_dx_tmp.v, v_dx_tmp.v);
      v_r2.v = vec_fma(v_dy_tmp.v, v_dy_tmp.v, v_r2.v);
      v_r2.v = vec_fma(v_dz_tmp.v, v_dz_tmp.v, v_r2.v);
      v_r2_2.v = vec_mul(v_dx_tmp2.v, v_dx_tmp2.v);
      v_r2_2.v = vec_fma(v_dy_tmp2.v, v_dy_tmp2.v, v_r2_2.v);
      v_r2_2.v = vec_fma(v_dz_tmp2.v, v_dz_tmp2.v, v_r2_2.v);

      v2_r2.v = vec_mul(v_dx2_tmp.v, v_dx2_tmp.v);
      v2_r2.v = vec_fma(v_dy2_tmp.v, v_dy2_tmp.v, v2_r2.v);
      v2_r2.v = vec_fma(v_dz2_tmp.v, v_dz2_tmp.v, v2_r2.v);
      v2_r2_2.v = vec_mul(v_dx2_tmp2.v, v_dx2_tmp2.v);
      v2_r2_2.v = vec_fma(v_dy2_tmp2.v, v_dy2_tmp2.v, v2_r2_2.v);
      v2_r2_2.v = vec_fma(v_dz2_tmp2.v, v_dz2_tmp2.v, v2_r2_2.v);

/* Form a mask from r2 < hig2 and r2 > 0.*/
#ifdef HAVE_AVX512_F
      // KNL_MASK_16 doi_mask, doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi2_mask_check, doi2_mask2, doi2_mask2_check;

      doi_mask_check = vec_cmp_gt(v_r2.v, vec_setzero());
      doi_mask = vec_cmp_lt(v_r2.v, v_hig2.v);

      doi2_mask_check = vec_cmp_gt(v2_r2.v, vec_setzero());
      doi2_mask = vec_cmp_lt(v2_r2.v, v_hig2_2.v);

      doi_mask2_check = vec_cmp_gt(v_r2_2.v, vec_setzero());
      doi_mask2 = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      doi2_mask2_check = vec_cmp_gt(v2_r2_2.v, vec_setzero());
      doi2_mask2 = vec_cmp_lt(v2_r2_2.v, v_hig2_2.v);

      doi_mask = doi_mask & doi_mask_check;
      doi_mask2 = doi_mask2 & doi_mask2_check;

      doi2_mask = doi2_mask & doi2_mask_check;
      doi2_mask2 = doi2_mask2 & doi2_mask2_check;
#else
      vector v_doi_mask, v_doi_mask_check, v_doi_mask2, v_doi_mask2_check;
      int doi_mask2;

      vector v_doi2_mask, v_doi2_mask_check, v_doi2_mask2, v_doi2_mask2_check;
      int doi2_mask2;

      v_doi_mask_check.v = vec_cmp_gt(v_r2.v, vec_setzero());
      v_doi_mask.v = vec_cmp_lt(v_r2.v, v_hig2.v);

      v_doi2_mask_check.v = vec_cmp_gt(v2_r2.v, vec_setzero());
      v_doi2_mask.v = vec_cmp_lt(v2_r2.v, v_hig2_2.v);

      v_doi_mask2_check.v = vec_cmp_gt(v_r2_2.v, vec_setzero());
      v_doi_mask2.v = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      v_doi2_mask2_check.v = vec_cmp_gt(v2_r2_2.v, vec_setzero());
      v_doi2_mask2.v = vec_cmp_lt(v2_r2_2.v, v_hig2_2.v);

      doi_mask = vec_cmp_result(vec_and(v_doi_mask.v, v_doi_mask_check.v));
      doi_mask2 = vec_cmp_result(vec_and(v_doi_mask2.v, v_doi_mask2_check.v));
      doi2_mask = vec_cmp_result(vec_and(v_doi2_mask.v, v_doi2_mask_check.v));
      doi2_mask2 =
          vec_cmp_result(vec_and(v_doi2_mask2.v, v_doi2_mask2_check.v));
#endif /* HAVE_AVX512_F */

      /* Hit or miss? */
      // if (doi_mask) {
      storeInteractions(doi_mask, pjd, &v_r2, &v_dx_tmp, &v_dy_tmp, &v_dz_tmp,
                        &mj, &pjvx, &pjvy, &pjvz, cell_cache, &int_cache,
                        &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum,
                        &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv,
                        v_vix, v_viy, v_viz);
      //}
      // if (doi2_mask) {
      storeInteractions(
          doi2_mask, pjd, &v2_r2, &v_dx2_tmp, &v_dy2_tmp, &v_dz2_tmp, &mj,
          &pjvx, &pjvy, &pjvz, cell_cache, &int_cache2, &icount2, &rhoSum2,
          &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2,
          &curlvySum2, &curlvzSum2, v_hi_inv2, v_vix2, v_viy2, v_viz2);
      //}
      /* Hit or miss? */
      // if (doi_mask2) {
      storeInteractions(doi_mask2, pjd + VEC_SIZE, &v_r2_2, &v_dx_tmp2,
                        &v_dy_tmp2, &v_dz_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2,
                        cell_cache, &int_cache, &icount, &rhoSum, &rho_dhSum,
                        &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum,
                        &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      //}
      // if (doi2_mask2) {
      storeInteractions(doi2_mask2, pjd + VEC_SIZE, &v2_r2_2, &v_dx2_tmp2,
                        &v_dy2_tmp2, &v_dz2_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2,
                        cell_cache, &int_cache2, &icount2, &rhoSum2,
                        &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2,
                        &curlvxSum2, &curlvySum2, &curlvzSum2, v_hi_inv2,
                        v_vix2, v_viy2, v_viz2);
      //}
    }

    /* Perform padded vector remainder interactions if any are present. */
    calcRemInteractions(&int_cache, icount, &rhoSum, &rho_dhSum, &wcountSum,
                        &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum,
                        &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz,
                        &icount_align);

    calcRemInteractions(&int_cache2, icount2, &rhoSum2, &rho_dhSum2,
                        &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2,
                        &curlvySum2, &curlvzSum2, v_hi_inv2, v_vix2, v_viy2,
                        v_viz2, &icount_align2);

    /* Initialise masks to true incase remainder interactions have been
     * performed. */
    vector int_mask, int_mask2;
    vector int2_mask, int2_mask2;
#ifdef HAVE_AVX512_F
    KNL_MASK_16 knl_mask = 0xFFFF;
    KNL_MASK_16 knl_mask2 = 0xFFFF;
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
    int2_mask.m = vec_setint1(0xFFFFFFFF);
    int2_mask2.m = vec_setint1(0xFFFFFFFF);
#else
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);

    int2_mask.m = vec_setint1(0xFFFFFFFF);
    int2_mask2.m = vec_setint1(0xFFFFFFFF);
#endif

    /* Perform interaction with 2 vectors. */
    for (int pjd = 0; pjd < icount_align; pjd += (NUM_VEC_PROC * VEC_SIZE)) {
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

    for (int pjd = 0; pjd < icount_align2; pjd += (NUM_VEC_PROC * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(
          &int_cache2.r2q[pjd], &int_cache2.dxq[pjd], &int_cache2.dyq[pjd],
          &int_cache2.dzq[pjd], v_hi_inv2, v_vix2, v_viy2, v_viz2,
          &int_cache2.vxq[pjd], &int_cache2.vyq[pjd], &int_cache2.vzq[pjd],
          &int_cache2.mq[pjd], &rhoSum2, &rho_dhSum2, &wcountSum2,
          &wcount_dhSum2, &div_vSum2, &curlvxSum2, &curlvySum2, &curlvzSum2,
          int2_mask, int2_mask2,
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

    VEC_HADD(rhoSum2, pi2->rho);
    VEC_HADD(rho_dhSum2, pi2->density.rho_dh);
    VEC_HADD(wcountSum2, pi2->density.wcount);
    VEC_HADD(wcount_dhSum2, pi2->density.wcount_dh);
    VEC_HADD(div_vSum2, pi2->density.div_v);
    VEC_HADD(curlvxSum2, pi2->density.rot_v[0]);
    VEC_HADD(curlvySum2, pi2->density.rot_v[1]);
    VEC_HADD(curlvzSum2, pi2->density.rot_v[2]);

    /* Reset interaction count. */
    icount = 0;
    icount2 = 0;
  } /* loop over all particles. */

  TIMER_TOC(timer_doself_density);
#endif /* WITH_VECTORIZATION */
}
