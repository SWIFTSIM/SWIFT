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

#define array_align sizeof(float) * VEC_SIZE
#define NUM_VEC_PROC 2
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)

#ifdef WITH_VECTORIZATION
__attribute__((always_inline)) INLINE static void calcRemInteractions(const struct cache *const cell_cache, float *r2q, float *dxq, float *dyq, float *dzq, float *mq, float *vxq, float *vyq, float *vzq, const int icount, vector *rhoSum, vector *rho_dhSum, vector *wcountSum, vector *wcount_dhSum, vector *div_vSum, vector *curlvxSum,vector *curlvySum, vector *curlvzSum, vector v_hi_inv, vector v_vix, vector v_viy, vector v_viz, int *icount_align) {

#ifdef HAVE_AVX512_F
  KNL_MASK_16 knl_mask, knl_mask2;
#endif
  vector int_mask, int_mask2;
  
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
    /* Pad secondary cache  */
    for(int i=icount; i<*icount_align; i++) {
      mq[i] = 0.f;
      r2q[i] = 1.f;
      dxq[i] = 0.f;
      dyq[i] = 0.f;
      dzq[i] = 0.f;
      vxq[i] = 0.f;
      vyq[i] = 0.f;
      vzq[i] = 0.f;
    }

    /* Zero parts of mask that represent the padded values.*/
    if (pad < VEC_SIZE) {
#ifdef HAVE_AVX512_F
      knl_mask2 = knl_mask2 >> pad;
#else
      for(int i=VEC_SIZE - pad; i<VEC_SIZE; i++) int_mask2.i[i] = 0;
#endif
    }
    else {
#ifdef HAVE_AVX512_F
      knl_mask = knl_mask >> (VEC_SIZE - rem);
      knl_mask2 = 0;
#else
      for(int i=rem; i<VEC_SIZE; i++) int_mask.i[i] = 0;
      int_mask2.v = vec_setzero();
#endif
    }

    /* Perform remainder interaction and remove remainder from aligned interaction count. */
    *icount_align = icount - rem;
    runner_iact_nonsym_2_vec_density(&r2q[*icount_align], &dxq[*icount_align], &dyq[*icount_align], &dzq[*icount_align], v_hi_inv, v_vix, v_viy, v_viz, &vxq[*icount_align], &vyq[*icount_align], &vzq[*icount_align], &mq[*icount_align], rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum, int_mask, int_mask2,
#ifdef HAVE_AVX512_F
    knl_mask, knl_mask2);
#else
    0,0);
#endif
  }
}

__attribute__((always_inline)) INLINE static void storeInteractions(const int mask, const int pjd, vector *v_r2, vector *v_dx, vector *v_dy, vector *v_dz, vector *v_mj, vector *v_vjx, vector *v_vjy, vector *v_vjz, const struct cache *const cell_cache, float *r2q, float *dxq, float *dyq, float *dzq, float *mq, float *vxq, float *vyq, float *vzq, int *icount, vector *rhoSum, vector *rho_dhSum, vector *wcountSum, vector *wcount_dhSum, vector *div_vSum, vector *curlvxSum,vector *curlvySum, vector *curlvzSum, vector v_hi_inv, vector v_vix, vector v_viy, vector v_viz) {

#if defined(HAVE_AVX2) || defined(HAVE_AVX512_F)
  int pack = 0;

#ifdef HAVE_AVX512_F
  pack += __builtin_popcount(mask);
  VEC_LEFT_PACK(v_r2->v,mask,&r2q[*icount]);
  VEC_LEFT_PACK(v_dx->v,mask,&dxq[*icount]);
  VEC_LEFT_PACK(v_dy->v,mask,&dyq[*icount]);
  VEC_LEFT_PACK(v_dz->v,mask,&dzq[*icount]);
  VEC_LEFT_PACK(v_mj->v,mask,&mq[*icount]);
  VEC_LEFT_PACK(v_vjx->v,mask,&vxq[*icount]);
  VEC_LEFT_PACK(v_vjy->v,mask,&vyq[*icount]);
  VEC_LEFT_PACK(v_vjz->v,mask,&vzq[*icount]);
#else
  vector v_mask;
  VEC_FORM_PACKED_MASK(mask,v_mask.m,pack);
  
  VEC_LEFT_PACK(v_r2->v,v_mask.m,&r2q[*icount]);
  VEC_LEFT_PACK(v_dx->v,v_mask.m,&dxq[*icount]);
  VEC_LEFT_PACK(v_dy->v,v_mask.m,&dyq[*icount]);
  VEC_LEFT_PACK(v_dz->v,v_mask.m,&dzq[*icount]);
  VEC_LEFT_PACK(v_mj->v,v_mask.m,&mq[*icount]);
  VEC_LEFT_PACK(v_vjx->v,v_mask.m,&vxq[*icount]);
  VEC_LEFT_PACK(v_vjy->v,v_mask.m,&vyq[*icount]);
  VEC_LEFT_PACK(v_vjz->v,v_mask.m,&vzq[*icount]);
#endif

  (*icount) += pack;
#else
  for(int bit_index = 0; bit_index<VEC_SIZE; bit_index++) {
    if (mask & (1 << bit_index)) {
      /* Add this interaction to the queue. */
      r2q[*icount] = v_r2->f[bit_index];
      dxq[*icount] = v_dx->f[bit_index];
      dyq[*icount] = v_dy->f[bit_index];
      dzq[*icount] = v_dz->f[bit_index];
      mq[*icount] = cell_cache->m[pjd + bit_index];
      vxq[*icount] = cell_cache->vx[pjd + bit_index];
      vyq[*icount] = cell_cache->vy[pjd + bit_index];
      vzq[*icount] = cell_cache->vz[pjd + bit_index];

      (*icount)++;
    }
  }
  if(*icount >= (C2_CACHE_SIZE - (NUM_VEC_PROC * VEC_SIZE))) {

    int icount_align = *icount;
    calcRemInteractions(cell_cache, r2q, dxq, dyq, dzq, mq, vxq, vyq, vzq, *icount, rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum, v_hi_inv, v_vix, v_viy, v_viz, &icount_align);

    vector int_mask, int_mask2;
    int_mask.m = vec_setint1(0xFFFFFFFF);
    int_mask2.m = vec_setint1(0xFFFFFFFF);
    for (int pjd = 0; pjd < icount_align; pjd+=(NUM_VEC_PROC * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(&r2q[pjd], &dxq[pjd], &dyq[pjd], &dzq[pjd], v_hi_inv, v_vix, v_viy, v_viz, &vxq[pjd], &vyq[pjd], &vzq[pjd], &mq[pjd], rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum, int_mask, int_mask2, 0, 0);
    }
    *icount = 0;
  }

#endif
}
#endif

/**
 * @brief Compute the cell self-interaction (non-symmetric) vec.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void runner_doself1_density_vec(struct runner *r, struct cell *restrict c) {

#ifdef WITH_VECTORIZATION
  const int ti_current = r->e->ti_current;
  int doi_mask;
  struct part *restrict pi;
  int count_align;
  int num_vec_proc = NUM_VEC_PROC;

  struct part *restrict parts = c->parts;
  const int count = c->count;
  struct cache *restrict cell_cache = &r->par_cache;
  
  int icount = 0, icount_align = 0;
  float r2q[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dxq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dyq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dzq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float mq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vxq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vyq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vzq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  
  vector v_hi, v_vix, v_viy, v_viz, v_hig2, v_r2;

  //TIMER_TIC

  if (c->ti_end_min > ti_current) return;
  if (c->ti_end_max < ti_current) error("Cell in an impossible time-zone");

  if(cell_cache->count < count) {
    cache_init(cell_cache,count);
  }

  cache_read_particles(c,cell_cache);

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    pi = &parts[pid];

    /* Is the ith particle active? */
    if (pi->ti_end > ti_current) continue;

    vector pix, piy, piz;

    const float hi = cell_cache->h[pid];

    /* Fill pi position vector. */
    pix.v = vec_set1(cell_cache->x[pid]);
    piy.v = vec_set1(cell_cache->y[pid]);
    piz.v = vec_set1(cell_cache->z[pid]);
    v_hi.v = vec_set1(hi);
    v_vix.v = vec_set1(cell_cache->vx[pid]);
    v_viy.v = vec_set1(cell_cache->vy[pid]);
    v_viz.v = vec_set1(cell_cache->vz[pid]);

    const float hig2 = hi * hi * kernel_gamma2;
    v_hig2.v = vec_set1(hig2);

    vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum;
    
    vector v_hi_inv;
    
    VEC_RECIPROCAL(v_hi.v, v_hi_inv.v);
    
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
      /* Set positions to the same as particle pi so when the r2 > 0 mask is applied these extra contributions are masked out.*/
      for(int i=count; i<count_align; i++) {
        cell_cache->x[i] = pix.f[0];
        cell_cache->y[i] = piy.f[0];
        cell_cache->z[i] = piz.f[0];
      }
    }

    vector pjx, pjy, pjz;
    vector pjvx, pjvy, pjvz, mj;
    vector pjx2, pjy2, pjz2;
    vector pjvx2, pjvy2, pjvz2, mj2;

    /* Find all of particle pi's interacions and store needed values in secondary cache.*/
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
      mj2.v = vec_load(&cell_cache->m[pjd +VEC_SIZE]);

      /* Compute the pairwise distance. */
      vector v_dx_tmp, v_dy_tmp, v_dz_tmp;
      vector v_dx_tmp2, v_dy_tmp2, v_dz_tmp2, v_r2_2;

      v_dx_tmp.v = vec_sub(pix.v,pjx.v);
      v_dy_tmp.v = vec_sub(piy.v,pjy.v);
      v_dz_tmp.v = vec_sub(piz.v,pjz.v);
      v_dx_tmp2.v = vec_sub(pix.v,pjx2.v);
      v_dy_tmp2.v = vec_sub(piy.v,pjy2.v);
      v_dz_tmp2.v = vec_sub(piz.v,pjz2.v);
      
      v_r2.v = vec_mul(v_dx_tmp.v,v_dx_tmp.v);
      v_r2.v = vec_fma(v_dy_tmp.v,v_dy_tmp.v,v_r2.v);
      v_r2.v = vec_fma(v_dz_tmp.v,v_dz_tmp.v,v_r2.v);
      v_r2_2.v = vec_mul(v_dx_tmp2.v,v_dx_tmp2.v);
      v_r2_2.v = vec_fma(v_dy_tmp2.v,v_dy_tmp2.v,v_r2_2.v);
      v_r2_2.v = vec_fma(v_dz_tmp2.v,v_dz_tmp2.v,v_r2_2.v);
      
      /* Form a mask from r2 < hig2 and r2 > 0.*/
#ifdef HAVE_AVX512_F
      //KNL_MASK_16 doi_mask, doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi_mask_check, doi_mask2, doi_mask2_check;

      doi_mask_check = vec_cmp_gt(v_r2.v,vec_setzero());
      doi_mask = vec_cmp_lt(v_r2.v, v_hig2.v);

      doi_mask2_check = vec_cmp_gt(v_r2_2.v,vec_setzero());
      doi_mask2 = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      doi_mask = doi_mask & doi_mask_check;
      doi_mask2 = doi_mask2 & doi_mask2_check;

#else
      vector v_doi_mask, v_doi_mask_check, v_doi_mask2, v_doi_mask2_check;
      int doi_mask2;

      v_doi_mask_check.v = vec_cmp_gt(v_r2.v,vec_setzero());
      v_doi_mask.v = vec_cmp_lt(v_r2.v, v_hig2.v);

      v_doi_mask2_check.v = vec_cmp_gt(v_r2_2.v,vec_setzero());
      v_doi_mask2.v = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      doi_mask = vec_cmp_result(vec_and(v_doi_mask.v, v_doi_mask_check.v));
      doi_mask2 = vec_cmp_result(vec_and(v_doi_mask2.v, v_doi_mask2_check.v));
#endif

      /* Hit or miss? */
      if (doi_mask) {
        storeInteractions(doi_mask,pjd, &v_r2, &v_dx_tmp,&v_dy_tmp, &v_dz_tmp, &mj, &pjvx, &pjvy, &pjvz, cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      }
      /* Hit or miss? */
      if (doi_mask2) {
        storeInteractions(doi_mask2,pjd + VEC_SIZE, &v_r2_2, &v_dx_tmp2,&v_dy_tmp2, &v_dz_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2, cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      }
    }

    /* Perform padded vector remainder interactions if any are present. */    
    calcRemInteractions(cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz, &icount_align);
    
    /* Initialise masks to true incase remainder interactions have been performed. */
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
    for (int pjd = 0; pjd < icount_align; pjd+=(num_vec_proc * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(&r2q[pjd], &dxq[pjd], &dyq[pjd], &dzq[pjd], v_hi_inv, v_vix, v_viy, v_viz, &vxq[pjd], &vyq[pjd], &vzq[pjd], &mq[pjd], &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, int_mask, int_mask2,
#ifdef HAVE_AVX512_F
      knl_mask, knl_mask2);
#else
      0, 0);      
#endif
    }

    /* Perform horizontal adds on vector sums and store result in particle pi. */
    VEC_HADD(rhoSum,pi->rho);
    VEC_HADD(rho_dhSum,pi->density.rho_dh);
    VEC_HADD(wcountSum,pi->density.wcount);
    VEC_HADD(wcount_dhSum,pi->density.wcount_dh);
    VEC_HADD(div_vSum,pi->density.div_v);
    VEC_HADD(curlvxSum,pi->density.rot_v[0]);
    VEC_HADD(curlvySum,pi->density.rot_v[1]);
    VEC_HADD(curlvzSum,pi->density.rot_v[2]);

    /* Reset interaction count. */
    icount = 0;
  } /* loop over all particles. */

  //TIMER_TOC(TIMER_DOSELF);
#endif
}

/**
 * @brief Compute the cell self-interaction (non-symmetric) vec.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void runner_doself1_density_vec_2(struct runner *r, struct cell *restrict c) {

#ifdef WITH_VECTORIZATION
  const int ti_current = r->e->ti_current;
  int doi_mask;
  int doi2_mask;
  struct part *restrict pi;
  struct part *restrict pi2;
  int count_align;

  int icount = 0, icount_align = 0;
  float r2q[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dxq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dyq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dzq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float mq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vxq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vyq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vzq[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  
  int icount2 = 0, icount_align2 = 0;
  float r2q2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dxq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dyq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float dzq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float mq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vxq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vyq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));
  float vzq2[C2_CACHE_SIZE] __attribute__((aligned(array_align)));

  vector v_hi, v_vix, v_viy, v_viz, v_hig2, v_r2;
  vector v_hi2, v_vix2, v_viy2, v_viz2, v_hig2_2, v2_r2;

  //TIMER_TIC

  if (c->ti_end_min > ti_current) return;
  if (c->ti_end_max < ti_current) error("Cell in an impossible time-zone");

  struct part *restrict parts = c->parts;
  const int count = c->count;
  struct cache *restrict cell_cache = &r->par_cache;

  if(cell_cache->count < count) {
    cache_init(cell_cache,count);
  }

  cache_read_particles(c,&r->par_cache);

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid+=2) {

    /* Get a pointer to the ith particle. */
    pi = &parts[pid];
    pi2 = &parts[pid + 1];

    /* Is the ith particle active? */
    if (pi->ti_end > ti_current) continue;

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

    vector rhoSum, rho_dhSum, wcountSum, wcount_dhSum, div_vSum, curlvxSum, curlvySum, curlvzSum;
    vector rhoSum2, rho_dhSum2, wcountSum2, wcount_dhSum2, div_vSum2, curlvxSum2, curlvySum2, curlvzSum2;
    
    vector v_hi_inv, v_hi_inv2;
    
    VEC_RECIPROCAL(v_hi.v, v_hi_inv.v);
    VEC_RECIPROCAL(v_hi2.v, v_hi_inv2.v);
    
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
      /* Set positions to the same as particle pi so when the r2 > 0 mask is applied these extra contributions are masked out.*/
      for(int i=count; i<count_align; i++) {
        cell_cache->x[i] = pix.f[0];
        cell_cache->y[i] = piy.f[0];
        cell_cache->z[i] = piz.f[0];
      }
    }

    vector pjx, pjy, pjz;
    vector pjvx, pjvy, pjvz, mj;
    vector pjx2, pjy2, pjz2;
    vector pjvx2, pjvy2, pjvz2, mj2;

    /* Find all of particle pi's interacions and store needed values in secondary cache.*/
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
      mj2.v = vec_load(&cell_cache->m[pjd +VEC_SIZE]);

      /* Compute the pairwise distance. */
      vector v_dx_tmp, v_dy_tmp, v_dz_tmp;
      vector v_dx_tmp2, v_dy_tmp2, v_dz_tmp2, v_r2_2;
      vector v_dx2_tmp, v_dy2_tmp, v_dz2_tmp;
      vector v_dx2_tmp2, v_dy2_tmp2, v_dz2_tmp2, v2_r2_2;

      v_dx_tmp.v = vec_sub(pix.v,pjx.v);
      v_dy_tmp.v = vec_sub(piy.v,pjy.v);
      v_dz_tmp.v = vec_sub(piz.v,pjz.v);
      v_dx_tmp2.v = vec_sub(pix.v,pjx2.v);
      v_dy_tmp2.v = vec_sub(piy.v,pjy2.v);
      v_dz_tmp2.v = vec_sub(piz.v,pjz2.v);
      
      v_dx2_tmp.v = vec_sub(pix2.v,pjx.v);
      v_dy2_tmp.v = vec_sub(piy2.v,pjy.v);
      v_dz2_tmp.v = vec_sub(piz2.v,pjz.v);
      v_dx2_tmp2.v = vec_sub(pix2.v,pjx2.v);
      v_dy2_tmp2.v = vec_sub(piy2.v,pjy2.v);
      v_dz2_tmp2.v = vec_sub(piz2.v,pjz2.v);

      v_r2.v = vec_mul(v_dx_tmp.v,v_dx_tmp.v);
      v_r2.v = vec_fma(v_dy_tmp.v,v_dy_tmp.v,v_r2.v);
      v_r2.v = vec_fma(v_dz_tmp.v,v_dz_tmp.v,v_r2.v);
      v_r2_2.v = vec_mul(v_dx_tmp2.v,v_dx_tmp2.v);
      v_r2_2.v = vec_fma(v_dy_tmp2.v,v_dy_tmp2.v,v_r2_2.v);
      v_r2_2.v = vec_fma(v_dz_tmp2.v,v_dz_tmp2.v,v_r2_2.v);
      
      v2_r2.v = vec_mul(v_dx2_tmp.v,v_dx2_tmp.v);
      v2_r2.v = vec_fma(v_dy2_tmp.v,v_dy2_tmp.v,v2_r2.v);
      v2_r2.v = vec_fma(v_dz2_tmp.v,v_dz2_tmp.v,v2_r2.v);
      v2_r2_2.v = vec_mul(v_dx2_tmp2.v,v_dx2_tmp2.v);
      v2_r2_2.v = vec_fma(v_dy2_tmp2.v,v_dy2_tmp2.v,v2_r2_2.v);
      v2_r2_2.v = vec_fma(v_dz2_tmp2.v,v_dz2_tmp2.v,v2_r2_2.v);

      /* Form a mask from r2 < hig2 and r2 > 0.*/
#ifdef HAVE_AVX512_F
      //KNL_MASK_16 doi_mask, doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi_mask_check, doi_mask2, doi_mask2_check;
      KNL_MASK_16 doi2_mask_check, doi2_mask2, doi2_mask2_check;

      doi_mask_check = vec_cmp_gt(v_r2.v,vec_setzero());
      doi_mask = vec_cmp_lt(v_r2.v, v_hig2.v);

      doi2_mask_check = vec_cmp_gt(v2_r2.v,vec_setzero());
      doi2_mask = vec_cmp_lt(v2_r2.v, v_hig2_2.v);

      doi_mask2_check = vec_cmp_gt(v_r2_2.v,vec_setzero());
      doi_mask2 = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      doi2_mask2_check = vec_cmp_gt(v2_r2_2.v,vec_setzero());
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

      v_doi_mask_check.v = vec_cmp_gt(v_r2.v,vec_setzero());
      v_doi_mask.v = vec_cmp_lt(v_r2.v, v_hig2.v);

      v_doi2_mask_check.v = vec_cmp_gt(v2_r2.v,vec_setzero());
      v_doi2_mask.v = vec_cmp_lt(v2_r2.v, v_hig2_2.v);

      v_doi_mask2_check.v = vec_cmp_gt(v_r2_2.v,vec_setzero());
      v_doi_mask2.v = vec_cmp_lt(v_r2_2.v, v_hig2.v);

      v_doi2_mask2_check.v = vec_cmp_gt(v2_r2_2.v,vec_setzero());
      v_doi2_mask2.v = vec_cmp_lt(v2_r2_2.v, v_hig2_2.v);

      doi_mask = vec_cmp_result(vec_and(v_doi_mask.v, v_doi_mask_check.v));
      doi_mask2 = vec_cmp_result(vec_and(v_doi_mask2.v, v_doi_mask2_check.v));
      doi2_mask = vec_cmp_result(vec_and(v_doi2_mask.v, v_doi2_mask_check.v));
      doi2_mask2 = vec_cmp_result(vec_and(v_doi2_mask2.v, v_doi2_mask2_check.v));
#endif

      /* Hit or miss? */
      //if (doi_mask) {
        storeInteractions(doi_mask,pjd, &v_r2, &v_dx_tmp,&v_dy_tmp, &v_dz_tmp, &mj, &pjvx, &pjvy, &pjvz, cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      //}
      //if (doi2_mask) {
        storeInteractions(doi2_mask,pjd, &v2_r2, &v_dx2_tmp,&v_dy2_tmp, &v_dz2_tmp, &mj, &pjvx, &pjvy, &pjvz, cell_cache, &r2q2[0], &dxq2[0], &dyq2[0], &dzq2[0], &mq2[0], &vxq2[0], &vyq2[0], &vzq2[0], &icount2, &rhoSum2, &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2, &curlvySum2, &curlvzSum2, v_hi_inv2, v_vix2, v_viy2, v_viz2);
      //}       
      /* Hit or miss? */
      //if (doi_mask2) {
        storeInteractions(doi_mask2,pjd + VEC_SIZE, &v_r2_2, &v_dx_tmp2,&v_dy_tmp2, &v_dz_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2, cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], &icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz);
      //}
      //if (doi2_mask2) {
        storeInteractions(doi2_mask2,pjd + VEC_SIZE, &v2_r2_2, &v_dx2_tmp2,&v_dy2_tmp2, &v_dz2_tmp2, &mj2, &pjvx2, &pjvy2, &pjvz2, cell_cache, &r2q2[0], &dxq2[0], &dyq2[0], &dzq2[0], &mq2[0], &vxq2[0], &vyq2[0], &vzq2[0], &icount2, &rhoSum2, &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2, &curlvySum2, &curlvzSum2, v_hi_inv2, v_vix2, v_viy2, v_viz2);
      //}
    }

    /* Perform padded vector remainder interactions if any are present. */    
    calcRemInteractions(cell_cache, &r2q[0], &dxq[0], &dyq[0], &dzq[0], &mq[0], &vxq[0], &vyq[0], &vzq[0], icount, &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, v_hi_inv, v_vix, v_viy, v_viz, &icount_align);
    
    calcRemInteractions(cell_cache, &r2q2[0], &dxq2[0], &dyq2[0], &dzq2[0], &mq2[0], &vxq2[0], &vyq2[0], &vzq2[0], icount2, &rhoSum2, &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2, &curlvySum2, &curlvzSum2, v_hi_inv2, v_vix2, v_viy2, v_viz2, &icount_align2);

    /* Initialise masks to true incase remainder interactions have been performed. */
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
    for (int pjd = 0; pjd < icount_align; pjd+=(NUM_VEC_PROC * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(&r2q[pjd], &dxq[pjd], &dyq[pjd], &dzq[pjd], v_hi_inv, v_vix, v_viy, v_viz, &vxq[pjd], &vyq[pjd], &vzq[pjd], &mq[pjd], &rhoSum, &rho_dhSum, &wcountSum, &wcount_dhSum, &div_vSum, &curlvxSum, &curlvySum, &curlvzSum, int_mask, int_mask2,
#ifdef HAVE_AVX512_F
      knl_mask, knl_mask2);
#else
      0, 0);      
#endif
    }

    for (int pjd = 0; pjd < icount_align2; pjd+=(NUM_VEC_PROC * VEC_SIZE)) {
      runner_iact_nonsym_2_vec_density(&r2q2[pjd], &dxq2[pjd], &dyq2[pjd], &dzq2[pjd], v_hi_inv2, v_vix2, v_viy2, v_viz2, &vxq2[pjd], &vyq2[pjd], &vzq2[pjd], &mq2[pjd], &rhoSum2, &rho_dhSum2, &wcountSum2, &wcount_dhSum2, &div_vSum2, &curlvxSum2, &curlvySum2, &curlvzSum2, int2_mask, int2_mask2,
#ifdef HAVE_AVX512_F
      knl_mask, knl_mask2);
#else
      0, 0);      
#endif
    }
    /* Perform horizontal adds on vector sums and store result in particle pi. */
    VEC_HADD(rhoSum,pi->rho);
    VEC_HADD(rho_dhSum,pi->density.rho_dh);
    VEC_HADD(wcountSum,pi->density.wcount);
    VEC_HADD(wcount_dhSum,pi->density.wcount_dh);
    VEC_HADD(div_vSum,pi->density.div_v);
    VEC_HADD(curlvxSum,pi->density.rot_v[0]);
    VEC_HADD(curlvySum,pi->density.rot_v[1]);
    VEC_HADD(curlvzSum,pi->density.rot_v[2]);

    VEC_HADD(rhoSum2,pi2->rho);
    VEC_HADD(rho_dhSum2,pi2->density.rho_dh);
    VEC_HADD(wcountSum2,pi2->density.wcount);
    VEC_HADD(wcount_dhSum2,pi2->density.wcount_dh);
    VEC_HADD(div_vSum2,pi2->density.div_v);
    VEC_HADD(curlvxSum2,pi2->density.rot_v[0]);
    VEC_HADD(curlvySum2,pi2->density.rot_v[1]);
    VEC_HADD(curlvzSum2,pi2->density.rot_v[2]);

    /* Reset interaction count. */
    icount = 0;
    icount2 = 0;
  } /* loop over all particles. */

  //TIMER_TOC(TIMER_DOSELF);
#endif

}
