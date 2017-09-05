/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

// NEED TO USE
//   gpuErrchk( cudaPeekAtLastError() );
//  gpuErrchk( cudaDeviceSynchronize() );

#ifndef static
#define static
#endif
#ifndef restrict
#define restrict __restrict__
#endif
#define FULL_COMP
extern "C" {
#include "testcuda.h"
/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

}

#include <cuda.h>


/* Host function to check cuda functions don't return errors */
__host__ inline void cudaErrCheck(cudaError_t status) {
  if (status != cudaSuccess) {
    printf("%s\n", cudaGetErrorString(status));
  }
}

/* Task function to perform a self force task. No symmetry */
__device__ void doself_force(struct cell_cuda *ci) {
  /* Is the cell active? */
  if (!cuda_cell_is_active(ci)) return;

  const int count_i = ci->part_count;
  int part_i = ci->first_part;

  float3 a_hydro;
  float h_dt, v_sig_stor, entropy_dt;

  /* Loop over the particles */
  for (int pid = part_i + threadIdx.x; pid < part_i + count_i;
       pid += blockDim.x) {

    const float hi = cuda_parts.h[pid];
    if (!cuda_part_is_active(pid)) continue;
    /* Reset the values. */
    a_hydro.x = 0.0f;
    a_hydro.y = 0.0f;
    a_hydro.z = 0.0f;
    h_dt = 0.0f;
    v_sig_stor = cuda_parts.v_sig[pid];
    entropy_dt = 0.0f;

    double pix[3];
    pix[0] = cuda_parts.x_x[pid];
    pix[1] = cuda_parts.x_y[pid];
    pix[2] = cuda_parts.x_z[pid];

    const float hig2 = hi * hi * kernel_gamma2;

    /* Loop over the particles in cj. */
    for (int pjd = part_i; pjd < part_i + count_i; pjd++) {

      if (pid == pjd) continue;
      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      dx[0] = pix[0] - cuda_parts.x_x[pjd];
      r2 += dx[0] * dx[0];
      dx[1] = pix[1] - cuda_parts.x_y[pjd];
      r2 += dx[1] * dx[1];
      dx[2] = pix[2] - cuda_parts.x_z[pjd];
      r2 += dx[2] * dx[2];
      const float hj = cuda_parts.h[pjd];
      if (r2 < hig2 || r2 < hj * hj * kernel_gamma2) {
        float wi, wj, wi_dx, wj_dx;
        const float fac_mu = 1.f;

        const float r = sqrtf(r2);
        const float r_inv = 1.0f / r;

        /* Load some data.*/
        const float mj = cuda_parts.mass[pjd];
        const float rhoi = cuda_parts.rho[pid];
        const float rhoj = cuda_parts.rho[pjd];

        /* Get the kernel for hi. */
        const float hi_inv = 1.0f / hi;
        const float hid_inv = cuda_pow_dimension_plus_one(hi_inv);
        const float ui = r * hi_inv;
        cuda_kernel_deval(ui, &wi, &wi_dx);
        const float wi_dr = hid_inv * wi_dx;

        /* Get the kernel for hj. */
        const float hj_inv = 1.0f / hj;
        const float hjd_inv = cuda_pow_dimension_plus_one(hj_inv);
        const float xj = r * hj_inv;
        cuda_kernel_deval(xj, &wj, &wj_dx);
        const float wj_dr = hjd_inv * wj_dx;

        /* Compute h-gradient terms */
        const float f_i = cuda_parts.f[pid];
        const float f_j = cuda_parts.f[pjd];

        /* Compute pressure terms */
        const float P_over_rho2_i = cuda_parts.P_over_rho2[pid];
        const float P_over_rho2_j = cuda_parts.P_over_rho2[pjd];

        /* Compute sound speeds*/
        const float ci = cuda_parts.soundspeed[pid];
        const float cj = cuda_parts.soundspeed[pjd];

        /* Compute dv dot r. */
        const float dvdr = (cuda_parts.v[pid].x - cuda_parts.v[pjd].x) * dx[0] +
                           (cuda_parts.v[pid].y - cuda_parts.v[pjd].y) * dx[1] +
                           (cuda_parts.v[pid].z - cuda_parts.v[pjd].z) * dx[2];

        /* Balsara term */
        const float balsara_i = cuda_parts.balsara[pid];
        const float balsara_j = cuda_parts.balsara[pjd];

        /* Are the particles moving towards each other? */
        const float omega_ij = (dvdr < 0.f) ? dvdr : 0.f;
        const float mu_ij = fac_mu * r_inv * omega_ij;

        /* Signal velocity */
        const float v_sig = ci + cj - 3.f * mu_ij;
        const float rho_ij = 0.5f * (rhoi + rhoj);
        const float visc = -0.25 * const_viscosity_alpha * v_sig * mu_ij *
                           (balsara_i + balsara_j) / rho_ij;

        /* Now convolce with the kernel */
        const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
        const float sph_term =
            (f_i * P_over_rho2_i * wi_dr + f_j * P_over_rho2_j * wj_dr) * r_inv;

        /* Compute the acceleration */
        const float acc = visc_term + sph_term;

        /* Compute the force */
        a_hydro.x -= mj * acc * dx[0];
        a_hydro.y -= mj * acc * dx[1];
        a_hydro.z -= mj * acc * dx[2];


        /* Get the time derivative for h. */
        h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

        /* Update the signal velocity. */
        v_sig_stor = (v_sig_stor > v_sig) ? v_sig_stor : v_sig;

        /* Change in entropy */
        entropy_dt += mj * visc_term * dvdr;
      }

    }  // Inner loop

    /* Flush to global stores.*/
    atomicAdd(&cuda_parts.a_hydro[pid].x, a_hydro.x);
    atomicAdd(&cuda_parts.a_hydro[pid].y, a_hydro.y);
    atomicAdd(&cuda_parts.a_hydro[pid].z, a_hydro.z);
    atomicAdd(&cuda_parts.h_dt[pid], h_dt);
    atomicAdd(&cuda_parts.entropy_dt[pid], entropy_dt);

    /* Update the signal velocity */
    float global_vsig = cuda_parts.v_sig[pid];
    int *address_as_int = (int *)&cuda_parts.v_sig[pid];
    int old = *address_as_int;
    int assumed;
    do {
      global_vsig = cuda_parts.v_sig[pid];  // Scary line.
      assumed = old;
      if (v_sig_stor > global_vsig)
        old = atomicCAS(address_as_int, assumed, __float_as_int(v_sig_stor));
    } while (assumed != old && v_sig_stor > global_vsig);

  }  // Outer loop
}

/* Task function to execute a self-density task. */
__device__ void doself_density(struct cell_cuda *ci) {

  /* Is the cell active? */
  if (!cuda_cell_is_active(ci)) {
    printf(
        "Cell isn't active..., ti_end_min=%i, ti_current=%i, "
        "max_active_bin=%i\n",
        ci->ti_end_min, ti_current, max_active_bin);
    return;
  }

  const int count_i = ci->part_count;
  int part_i = ci->first_part;
  float rho, rho_dh, div_v, wcount, wcount_dh;
  float3 rot_v;

  for (int pid = part_i + threadIdx.x; pid < part_i + count_i;
       pid += blockDim.x) {
    double pix[3];
    pix[0] = cuda_parts.x_x[pid];
    pix[1] = cuda_parts.x_y[pid];
    pix[2] = cuda_parts.x_z[pid];
    const float hi = cuda_parts.h[pid];
    const float hig2 = hi * hi * kernel_gamma2;

    /* Reset local values. */
    rho = 0.0f;
    rho_dh = 0.0f;
    div_v = 0.0f;
    wcount = 0.0f;
    wcount_dh = 0.0f;
    rot_v.x = 0.0f;
    rot_v.y = 0.0f;
    rot_v.z = 0.0f;

    /* If the particle isn't active skip it. */
    if (!cuda_part_is_active(pid)) {
      continue;
    }

    /* Search for the neighbours! */
    for (int pjd = part_i; pjd < part_i + count_i; pjd++) {
      /* Particles don't interact with themselves */
      if (pid == pjd) continue;
      float dx[3], r2 = 0.0f;
      dx[0] = pix[0] - cuda_parts.x_x[pjd];
      r2 += dx[0] * dx[0];
      dx[1] = pix[1] - cuda_parts.x_y[pjd];
      r2 += dx[1] * dx[1];
      dx[2] = pix[2] - cuda_parts.x_z[pjd];
      r2 += dx[2] * dx[2];

      /* If in range then interact. */
      if (r2 < hig2) {
        float w, dw_dx;
        float dv[3], curlvr[3];
        /* Load mass on particle pj. */
        const float mj = cuda_parts.mass[pjd];

        /* Get r and 1/r */
        const float r = sqrtf(r2);
        const float ri = 1.0f / r;

        /* Compute the kernel function */
        const float hi_inv = 1.0f / hi;
        const float ui = r * hi_inv;

        cuda_kernel_deval(ui, &w, &dw_dx);
        /* Compute contribution to the density. */
        rho += mj * w;
        rho_dh -= mj * (hydro_dimension * w + ui * dw_dx);

        /* Compute contribution to the number of neighbours */
        wcount += w;
        wcount_dh -= (hydro_dimension * w + ui * dw_dx);

        const float fac = mj * dw_dx * ri;

        /* Compute dv dot r */
        float3 piv, pjv;
        piv = cuda_parts.v[pid];
        pjv = cuda_parts.v[pjd];
        dv[0] = piv.x - pjv.x;
        dv[1] = piv.y - pjv.y;
        dv[2] = piv.z - pjv.z;
        const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

        div_v -= fac * dvdr;

        curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
        curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
        curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

        rot_v.x += fac * curlvr[0];
        rot_v.y += fac * curlvr[1];
        rot_v.z += fac * curlvr[2];
      }
    }  // Loop over cj.
    /* Write data for particle pid back to global stores. */
    atomicAdd(&cuda_parts.rho[pid], rho);
    atomicAdd(&cuda_parts.rho_dh[pid], rho_dh);
    atomicAdd(&cuda_parts.wcount[pid], wcount);
    atomicAdd(&cuda_parts.wcount_dh[pid], wcount_dh);
    atomicAdd(&cuda_parts.div_v[pid], div_v);
    atomicAdd(&cuda_parts.rot_v[pid].x, rot_v.x);
    atomicAdd(&cuda_parts.rot_v[pid].y, rot_v.y);
    atomicAdd(&cuda_parts.rot_v[pid].z, rot_v.z);
  }
}


__global__ void do_test(struct cell_cuda *ci, struct cell_cuda *cj) {

  doself_density(ci);
  doself_force(ci);
  //dopair_density(ci, cj);
}

void allocate_parts(struct particle_arrays *p, size_t num_part) { 
  p->id = (long long int*) malloc (sizeof(long long int) * num_part);
  p->x_x = (double*) malloc(sizeof(double) * num_part);
  p->x_y = (double*) malloc(sizeof(double) * num_part);
  p->x_z = (double*) malloc(sizeof(double) * num_part);
  p->v = (float3*) malloc(sizeof(float3) * num_part);
  p->a_hydro = (float3*) malloc(sizeof(float3) * num_part);
  p->h = (float*) malloc(sizeof(float) * num_part);
  p->mass = (float*) malloc(sizeof(float) * num_part);
  p->rho = (float*) malloc(sizeof(float) * num_part);
  p->entropy = (float*) malloc(sizeof(float) * num_part);
  p->entropy_dt = (float*) malloc(sizeof(float) * num_part);

  // density
  p->wcount = (float*) malloc(sizeof(float) * num_part);
  p->wcount_dh = (float*) malloc(sizeof(float) * num_part);
  p->rho_dh = (float*) malloc(sizeof(float) * num_part);
  p->rot_v = (float3*) malloc(sizeof(float3) * num_part);
  p->div_v = (float*) malloc(sizeof(float) * num_part);

  // force
  p->balsara = (float*) malloc(sizeof(float) * num_part);
  p->f = (float*) malloc(sizeof(float) * num_part);
  p->P_over_rho2 = (float*) malloc(sizeof(float) * num_part);
  p->soundspeed = (float*) malloc(sizeof(float) * num_part);
  p->v_sig = (float*) malloc(sizeof(float) * num_part);
  p->h_dt = (float*) malloc(sizeof(float) * num_part);
  
  p->time_bin = (timebin_t*) malloc(sizeof(timebin_t) * num_part);
}

void allocate_device_parts(size_t num_part) {
  struct particle_arrays c_parts;
  
  cudaErrCheck(cudaMalloc(&c_parts.id, sizeof(long long int) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.x_x, sizeof(double) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.x_y, sizeof(double) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.x_z, sizeof(double) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.v, sizeof(float3) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.a_hydro, sizeof(float3) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.h, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.mass, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.rho, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.entropy, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.entropy_dt, sizeof(float) * num_part));

  // density
  cudaErrCheck(cudaMalloc(&c_parts.wcount, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.wcount_dh, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.rho_dh, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.rot_v, sizeof(float3) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.div_v, sizeof(float) * num_part));

  // force
  cudaErrCheck(cudaMalloc(&c_parts.balsara, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.f, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.P_over_rho2, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.soundspeed, sizeof(float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.v_sig, sizeof(volatile float) * num_part));
  cudaErrCheck(cudaMalloc(&c_parts.h_dt, sizeof(float) * num_part));

  cudaErrCheck(cudaMalloc(&c_parts.time_bin, sizeof(timebin_t) * num_part));

  cudaErrCheck(cudaMemcpyToSymbol(cuda_parts, &c_parts, sizeof(struct particle_arrays)));
}

void free_parts(struct particle_arrays p) {
  free(p.id); p.id = NULL;
  free(p.x_x); p.x_x = NULL;
  free(p.x_y); p.x_y = NULL;
  free(p.x_z); p.x_z = NULL;
  free(p.v); p.v = NULL;
  free(p.a_hydro); p.a_hydro = NULL;
  free(p.h); p.h = NULL;
  free(p.mass); p.mass = NULL;
  free(p.rho); p.rho = NULL;
  free(p.entropy); p.entropy = NULL;
  free(p.entropy_dt); p.entropy_dt = NULL;

  // density
  free(p.wcount); p.wcount = NULL;
  free(p.wcount_dh); p.wcount_dh = NULL;
  free(p.rho_dh); p.rho_dh = NULL;
  free(p.rot_v); p.rot_v = NULL;
  free(p.div_v); p.div_v = NULL;

  // force
  free(p.balsara); p.balsara = NULL;
  free(p.f); p.f = NULL;
  free(p.P_over_rho2); p.P_over_rho2 = NULL;
  free(p.soundspeed); p.soundspeed = NULL;
  free((void *)p.v_sig); p.v_sig = NULL;
  free(p.h_dt); p.h_dt = NULL;
  
  free(p.time_bin); p.time_bin = NULL;
}

void free_device_parts() {
  struct particle_arrays parts;
  cudaErrCheck(cudaMemcpyFromSymbol(&parts, cuda_parts,
				    sizeof(struct particle_arrays)));

  cudaErrCheck(cudaFree(parts.id));
  cudaErrCheck(cudaFree(parts.x_x));
  cudaErrCheck(cudaFree(parts.x_y));
  cudaErrCheck(cudaFree(parts.x_z));
  cudaErrCheck(cudaFree(parts.v));
  cudaErrCheck(cudaFree(parts.a_hydro));
  cudaErrCheck(cudaFree(parts.h));
  cudaErrCheck(cudaFree(parts.mass));
  cudaErrCheck(cudaFree(parts.rho));
  cudaErrCheck(cudaFree(parts.entropy));
  cudaErrCheck(cudaFree(parts.entropy_dt));

  // density
  cudaErrCheck(cudaFree(parts.wcount));
  cudaErrCheck(cudaFree(parts.wcount_dh));
  cudaErrCheck(cudaFree(parts.rho_dh));
  cudaErrCheck(cudaFree(parts.rot_v));
  cudaErrCheck(cudaFree(parts.div_v));

  // force
  cudaErrCheck(cudaFree(parts.balsara));
  cudaErrCheck(cudaFree(parts.f));
  cudaErrCheck(cudaFree(parts.P_over_rho2));
  cudaErrCheck(cudaFree(parts.soundspeed));
  cudaErrCheck(cudaFree((void*)parts.v_sig));
  cudaErrCheck(cudaFree(parts.h_dt));

  
  cudaErrCheck(cudaFree(parts.time_bin));
}


void copy_to_device_array(struct cell *ci, int offset) {

  struct particle_arrays h_p;
    
  int num_part = ci->count;
  
  allocate_parts(&h_p, num_part);
  
  //copy particles data
  for(int i=0;i<num_part;i++) {
    int p_i = i + offset;
    struct part p = ci->parts[i];
    h_p.id[p_i] = p.id;
    h_p.x_x[p_i] = p.x[0];
    h_p.x_y[p_i] = p.x[1];
    h_p.x_z[p_i] = p.x[2];
    h_p.v[p_i].x = p.v[0];
    h_p.v[p_i].y = p.v[1];
    h_p.v[p_i].z = p.v[2];
    h_p.a_hydro[p_i].x = p.a_hydro[0];
    h_p.a_hydro[p_i].y = p.a_hydro[1];
    h_p.a_hydro[p_i].z = p.a_hydro[2];
    h_p.h[p_i] = p.h;
    h_p.mass[p_i] = p.mass;
    h_p.rho[p_i] = p.rho;
    h_p.entropy[p_i] = p.entropy;
    h_p.entropy_dt[p_i] = p.entropy_dt;

    // density
    h_p.wcount[p_i] = p.density.wcount;
    h_p.wcount_dh[p_i] = p.density.wcount_dh;
    h_p.rho_dh[p_i] = p.density.rho_dh;
    h_p.rot_v[p_i].x = p.density.rot_v[0];
    h_p.rot_v[p_i].y = p.density.rot_v[1];
    h_p.rot_v[p_i].z = p.density.rot_v[2];
    h_p.div_v[p_i] = p.density.div_v;

    // force
    h_p.balsara[p_i] = p.force.balsara;
    h_p.f[p_i] = p.force.f;
    h_p.P_over_rho2[p_i] = p.force.P_over_rho2;
    h_p.soundspeed[p_i] = p.force.soundspeed;
    h_p.v_sig[p_i] = p.force.v_sig;
    h_p.h_dt[p_i] = p.force.h_dt;
    
    h_p.time_bin[p_i] = p.time_bin;
  }

  struct particle_arrays c_parts;
  cudaErrCheck(cudaMemcpyFromSymbol(&c_parts, cuda_parts,
                                    sizeof(struct particle_arrays)));

  void *p_data = c_parts.id + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.id, sizeof(long long int) * num_part,cudaMemcpyHostToDevice));

  p_data = c_parts.x_x + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.x_x, sizeof(double) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.x_y + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.x_y, sizeof(double) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.x_z + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.x_z, sizeof(double) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.v + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.v, sizeof(float3) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.a_hydro + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.a_hydro,sizeof(float3) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.h + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.h,sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.mass + offset;
  cudaErrCheck(cudaMemcpy(p_data,h_p.mass, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.rho + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.rho,sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.entropy + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.entropy,sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.entropy_dt + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.entropy_dt,sizeof(float) * num_part, cudaMemcpyHostToDevice));

  // density
  p_data = c_parts.wcount + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.wcount, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.wcount_dh + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.wcount_dh, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.rho_dh + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.rho_dh, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.rot_v + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.rot_v, sizeof(float3) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.div_v + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.div_v, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  // force
  p_data = c_parts.balsara + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.balsara, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.f + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.f, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.P_over_rho2 + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.P_over_rho2, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.soundspeed + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.soundspeed, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  p_data = (void*)(c_parts.v_sig + offset);
  cudaErrCheck(cudaMemcpy(p_data, (void*)h_p.v_sig, sizeof(volatile float) * num_part, cudaMemcpyHostToDevice));

  p_data = c_parts.h_dt + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.h_dt, sizeof(float) * num_part, cudaMemcpyHostToDevice));

  
  p_data = c_parts.time_bin + offset;
  cudaErrCheck(cudaMemcpy(p_data, h_p.time_bin,sizeof(timebin_t) * num_part, cudaMemcpyHostToDevice));

  free_parts(h_p);
}

void copy_from_device_array(struct particle_arrays *h_p, int offset, size_t num_part) {

  struct particle_arrays c_parts;
  cudaErrCheck(cudaMemcpyFromSymbol(&c_parts, cuda_parts,
                                    sizeof(struct particle_arrays)));

  void *p_data = c_parts.id + offset;
  cudaErrCheck(cudaMemcpy(h_p->id, p_data, sizeof(long long int) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.x_x + offset;
  cudaErrCheck(cudaMemcpy(h_p->x_x, p_data, sizeof(double) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.x_y + offset;
  cudaErrCheck(cudaMemcpy(h_p->x_y, p_data, sizeof(double) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.x_z + offset;
  cudaErrCheck(cudaMemcpy(h_p->x_z, p_data, sizeof(double) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.v + offset;
  cudaErrCheck(cudaMemcpy(h_p->v, p_data, sizeof(float3) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.a_hydro + offset;
  cudaErrCheck(cudaMemcpy(h_p->a_hydro, p_data, sizeof(float3) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.h + offset;
  cudaErrCheck(cudaMemcpy(h_p->h, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.mass + offset;
  cudaErrCheck(cudaMemcpy(h_p->mass, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.rho + offset;
  cudaErrCheck(cudaMemcpy(h_p->rho, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.entropy + offset;
  cudaErrCheck(cudaMemcpy(h_p->entropy, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));

  p_data = c_parts.entropy_dt + offset;
  cudaErrCheck(cudaMemcpy(h_p->entropy_dt, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));

  // density
  p_data = c_parts.wcount + offset;
  cudaErrCheck(cudaMemcpy(h_p->wcount, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.wcount_dh + offset;
  cudaErrCheck(cudaMemcpy(h_p->wcount_dh, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.rho_dh + offset;
  cudaErrCheck(cudaMemcpy(h_p->rho_dh, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.rot_v + offset;
  cudaErrCheck(cudaMemcpy(h_p->rot_v, p_data, sizeof(float3) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.div_v + offset;
  cudaErrCheck(cudaMemcpy(h_p->div_v, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  

  // force
  p_data = c_parts.balsara + offset;
  cudaErrCheck(cudaMemcpy(h_p->balsara, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.f + offset;
  cudaErrCheck(cudaMemcpy(h_p->f, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.P_over_rho2 + offset;
  cudaErrCheck(cudaMemcpy(h_p->P_over_rho2, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.soundspeed + offset;
  cudaErrCheck(cudaMemcpy(h_p->soundspeed, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = (void*)(c_parts.v_sig + offset);
  cudaErrCheck(cudaMemcpy((void*)h_p->v_sig, p_data, sizeof(volatile float) * num_part, cudaMemcpyDeviceToHost));
  
  p_data = c_parts.h_dt + offset;
  cudaErrCheck(cudaMemcpy(h_p->h_dt, p_data, sizeof(float) * num_part, cudaMemcpyDeviceToHost));
  

  p_data = c_parts.time_bin + offset;
  cudaErrCheck(cudaMemcpy(h_p->time_bin, p_data, sizeof(timebin_t) * num_part, cudaMemcpyDeviceToHost));
}

struct cell_cuda* copy_from_host(struct cell *ci, int offset) {
  copy_to_device_array(ci, offset);
  /* Set the host pointer. */
  struct cell_cuda c2;

  for (int k=0; k<DIM; k++) {
    c2.loc[k] = ci->loc[k];
    c2.width[k] = ci->width[k];
  
  }
  c2.h_max = ci->h_max;
  c2.first_part = offset;
  c2.part_count = ci->count;
  if (ci->parent != NULL) {
    c2.parent = ci->parent->cuda_ID;
  } else {
    c2.parent = -1;
  }
  if (ci->super != NULL) {
    c2.super = ci->super->cuda_ID;
  } else {
    c2.super = -1;
  }
  c2.ti_end_min = ci->ti_end_min;
  c2.ti_end_max = ci->ti_end_max;
  c2.dmin = ci->dmin;
  c2.nr_links = 0;
  c2.split = ci->split;

  struct cell_cuda *ret;
  cudaErrCheck( cudaMalloc(&ret, sizeof(struct cell_cuda)) );
  cudaErrCheck( cudaMemcpy(ret, &c2, sizeof(struct cell_cuda), cudaMemcpyHostToDevice) );
  return ret;
}

void copy_to_host(struct cell_cuda *cuda_c, struct cell *c) {
  struct cell_cuda *h_cuda_cell = (struct cell_cuda*) malloc(sizeof(struct cell_cuda));

  cudaMemcpy(h_cuda_cell, cuda_c, sizeof(struct cell_cuda), cudaMemcpyDeviceToHost);
  int N = h_cuda_cell->part_count;
  int p_i = h_cuda_cell->first_part;

  struct particle_arrays h_cuda_part;
  allocate_parts(&h_cuda_part, N);
  copy_from_device_array(&h_cuda_part, p_i, N);

  for (int i=0; i<N; i++) {
    struct part *p = &c->parts[i];

    p->id = h_cuda_part.id[p_i];
    p->h = h_cuda_part.h[p_i];
    p->rho = h_cuda_part.rho[p_i];
    p->entropy = h_cuda_part.entropy[p_i];
    p->entropy_dt = h_cuda_part.entropy_dt[p_i];
    //density
    p->density.wcount = h_cuda_part.wcount[p_i];
    p->density.wcount_dh = h_cuda_part.wcount_dh[p_i];
    p->density.rho_dh = h_cuda_part.rho_dh[p_i];
    p->density.div_v = h_cuda_part.div_v[p_i];

    //force
    p->force.balsara = h_cuda_part.balsara[p_i];
    p->force.f = h_cuda_part.f[p_i];
    p->force.P_over_rho2 = h_cuda_part.P_over_rho2[p_i];
    p->force.soundspeed = h_cuda_part.soundspeed[p_i];
    p->force.v_sig = h_cuda_part.v_sig[p_i];
    p->force.h_dt = h_cuda_part.h_dt[p_i];
    
    p->x[0] = h_cuda_part.x_x[p_i];
    p->v[0] = h_cuda_part.v[p_i].x;
    p->a_hydro[0] = h_cuda_part.a_hydro[p_i].x;
    p->density.rot_v[0] = h_cuda_part.rot_v[p_i].x;
#if DIM > 1
    p->x[1] = h_cuda_part.x_y[p_i];
    p->v[1] = h_cuda_part.v[p_i].y;
    p->a_hydro[1] = h_cuda_part.a_hydro[p_i].y;
    p->density.rot_v[1] = h_cuda_part.rot_v[p_i].y;
#if DIM > 2
    p->x[2] = h_cuda_part.x_z[p_i];
    p->v[2] = h_cuda_part.v[p_i].z;
    p->a_hydro[2] = h_cuda_part.a_hydro[p_i].z;
    p->density.rot_v[2] = h_cuda_part.rot_v[p_i].z;
#endif
#endif

    p_i++;
  }
}


/* n is both particles per axis and box size:
 * particles are generated on a mesh with unit spacing
 */
struct cell *make_cell(size_t n, double *offset, double size, double h,
                       double density, unsigned long long *partId,
                       double pert) {
  const size_t count = n * n * n;
  const double volume = size * size * size;
  struct cell *cell = (struct cell*)malloc(sizeof(struct cell));

  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->parts, count * sizeof(struct part));

  /* Construct the parts */
  struct part *part = cell->parts;
  for (size_t x = 0; x < n; ++x) {
    for (size_t y = 0; y < n; ++y) {
      for (size_t z = 0; z < n; ++z) {
        part->x[0] = offset[0] +
	  size * (x + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        part->x[1] = offset[1] +
	  size * (y + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;
        part->x[2] = offset[2] +
	  size * (z + 0.5 + random_uniform(-0.5, 0.5) * pert) / (float)n;

	part->v[0] = random_uniform(-0.05, 0.05);
        part->v[1] = random_uniform(-0.05, 0.05);
        part->v[2] = random_uniform(-0.05, 0.05);
        part->h = size * h / (float)n;
        part->id = ++(*partId);
        part->mass = density * volume / count;
        part->time_bin = 1;

	part->density.wcount = 0.;
	part->density.wcount_dh = 0.;
	part->density.rho_dh = 0.;
	part->density.rot_v[0] = 0.;
	part->density.rot_v[1] = 0.;
	part->density.rot_v[2] = 0.;
	part->density.div_v = 0.;

	part->force.balsara = 0.;
	part->force.f = 0.;
	part->force.P_over_rho2 = 0.;
	part->force.soundspeed = 0.;
	part->force.v_sig = 0.;
	part->force.h_dt = 0.;

        ++part;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->h_max = h;
  cell->count = count;
  cell->dx_max_part = 0.;
  cell->dx_max_sort = 0.;
  cell->width[0] = n;
  cell->width[1] = n;
  cell->width[2] = n;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->ti_old_part = 8;
  cell->ti_end_min = 8;
  cell->ti_end_max = 8;

  //shuffle_particles(cell->parts, cell->count);

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->parts);
  //for (int k = 0; k < 13; k++)
  //if (ci->sort[k] != NULL) free(ci->sort[k]);
  free(ci);
}

/**
 * @brief Initializes all particles field to be ready for a density calculation
 */
void zero_particle_fields(struct cell *c) {
  for (int pid = 0; pid < c->count; pid++) {
    struct part* p = &c->parts[pid];
    p->rho = 0.f;

    // density
    p->density.wcount = 0.f;
    p->density.wcount_dh = 0.f;
    p->density.rho_dh = 0.f;
    p->density.div_v = 0.f;
    p->density.rot_v[0] = 0.f;
    p->density.rot_v[1] = 0.f;
    p->density.rot_v[2] = 0.f;

    // force
    p->force.balsara = 0.5f;
    p->force.f = 0.9f;
    p->force.P_over_rho2 = 1.2f;
    p->force.soundspeed = 0.f;
    p->force.v_sig = 0.f;
    p->force.h_dt = 0.f;
  }
}

/**
 * @brief Dump all the particles to a file
 */
void dump_particle_fields(char *fileName, struct cell *ci, struct cell *cj) {
  FILE *file = fopen(fileName, "w");

  /* Write header */
#ifdef FULL_COMP
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %13s %13s %13s %13s %13s "
          "%13s %13s %13s %13s %13s %13s %13s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho", "rho_dh",
          "wcount", "wcount_dh", "div_v", "curl_vx", "curl_vy", "curl_vz",
	  "a_hydro.x","a_hydro.y","a_hydro.z","h_dt","entropy_dt");
#else
  fprintf(file,
          "# %4s %10s %10s %10s %10s %10s %10s %13s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "rho");
#endif

  fprintf(file, "# ci --------------------------------------------\n");

  for (int pid = 0; pid < ci->count; pid++) {
#ifdef FULL_COMP
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e %13e %13e %13e %13e %13e %13e\n",
            ci->parts[pid].id, ci->parts[pid].x[0], ci->parts[pid].x[1],
            ci->parts[pid].x[2], ci->parts[pid].v[0], ci->parts[pid].v[1],
            ci->parts[pid].v[2], ci->parts[pid].rho,
            ci->parts[pid].density.rho_dh,
            ci->parts[pid].density.wcount, ci->parts[pid].density.wcount_dh,
            ci->parts[pid].density.div_v, ci->parts[pid].density.rot_v[0],
            ci->parts[pid].density.rot_v[1], ci->parts[pid].density.rot_v[2],
            ci->parts[pid].a_hydro[0], ci->parts[pid].a_hydro[1], ci->parts[pid].a_hydro[2],
            ci->parts[pid].force.h_dt, ci->parts[pid].entropy_dt
            );
#else
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %13e %13e\n",
            ci->parts[pid].id, ci->parts[pid].x[0], ci->parts[pid].x[1],
            ci->parts[pid].x[2], ci->parts[pid].v[0], ci->parts[pid].v[1],
            ci->parts[pid].v[2], ci->parts[pid].rho
            );
    
#endif
  }

  fprintf(file, "# cj --------------------------------------------\n");

  for (int pjd = 0; pjd < cj->count; pjd++) {
#ifdef FULL_COMP
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %13e %13e %13e %13e %13e "
            "%13e %13e %13e %13e %13e %13e %13e %13e %13e\n",
            cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
            cj->parts[pjd].x[2], cj->parts[pjd].v[0], cj->parts[pjd].v[1],
            cj->parts[pjd].v[2], cj->parts[pjd].rho,
            cj->parts[pjd].density.rho_dh,
            cj->parts[pjd].density.wcount, cj->parts[pjd].density.wcount_dh,
            cj->parts[pjd].density.div_v, cj->parts[pjd].density.rot_v[0],
            cj->parts[pjd].density.rot_v[1], cj->parts[pjd].density.rot_v[2],
            cj->parts[pjd].a_hydro[0], cj->parts[pjd].a_hydro[1], cj->parts[pjd].a_hydro[2],
            cj->parts[pjd].force.h_dt, cj->parts[pjd].entropy_dt
            );
#else
    fprintf(file,
            "%6llu %10f %10f %10f %10f %10f %13e %13e\n",
            cj->parts[pjd].id, cj->parts[pjd].x[0], cj->parts[pjd].x[1],
            cj->parts[pjd].x[2], cj->parts[pjd].v[0], cj->parts[pjd].v[1],
            cj->parts[pjd].v[2], cj->parts[pjd].rho
            );
#endif
  }
  fclose(file);
}


int main(int argc, char *argv[]) {

  int tmp = 8;
  cudaErrCheck(
      cudaMemcpyToSymbol(ti_current, &tmp, sizeof(int)));
  tmp = 1000;
  cudaErrCheck(
      cudaMemcpyToSymbol(max_active_bin, &tmp, sizeof(int)));

  float host_kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)];
  for (int a = 0; a < (kernel_degree + 1) * (kernel_ivals + 1); a++) {
    host_kernel_coeffs[a] = kernel_coeffs[a];
  }
  cudaErrCheck(cudaMemcpyToSymbol(
    cuda_kernel_coeffs, &host_kernel_coeffs,
    sizeof(float) * (kernel_degree + 1) * (kernel_ivals + 1)));

  size_t particles = 0, runs = 0, volume, type = 0;
  double offset[3] = {0, 0, 0}, h = 1.1255, size = 1., rho = 1.;
  double perturbation = 0.;
  struct cell *ci, *cj;

  char c;
  static unsigned long long partId = 0;
  char outputFileNameExtension[200] = "";
  char outputFileName[200] = "";

  //  ticks tic, toc, time;

  /* Initialize CPU frequency, this also starts time. */
  //unsigned long long cpufreq = 0;
  //clocks_set_cpufreq(cpufreq);

  /* Choke on FP-exceptions */
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  srand(0);

  while ((c = getopt(argc, argv, "h:p:r:t:d:f:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 'p':
        sscanf(optarg, "%zu", &particles);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case 't':
        sscanf(optarg, "%zu", &type);
        break;
      case 'd':
        sscanf(optarg, "%lf", &perturbation);
        break;
      case 'f':
        strcpy(outputFileNameExtension, optarg);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || particles == 0 || runs == 0 || type > 2) {
    printf(
        "\nUsage: %s -p PARTICLES_PER_AXIS -r NUMBER_OF_RUNS [OPTIONS...]\n"
        "\nGenerates a cell pair, filled with particles on a Cartesian grid."
        "\nThese are then interacted using runner_dopair1_density."
        "\n\nOptions:"
        "\n-t TYPE=0          - cells share face (0), edge (1) or corner (2)"
        "\n-h DISTANCE=1.1255 - smoothing length"
        "\n-d pert            - perturbation to apply to the particles [0,1["
        "\n-f fileName        - part of the file name used to save the dumps\n",
        argv[0]);
    exit(1);
  }

  volume = particles * particles * particles;
  message("particles: %zu B\npositions: 0 B", 2 * volume * sizeof(struct part));

  ci = make_cell(particles, offset, size, h, rho, &partId, perturbation);
  for (size_t i = 0; i < type + 1; ++i) offset[i] = 1.;
  cj = make_cell(particles, offset, size, h, rho, &partId, perturbation);

  sprintf(outputFileName, "swift_gpu_init_%s.dat", outputFileNameExtension);
  dump_particle_fields(outputFileName, ci, cj);

  allocate_device_parts(2*volume);

  struct cell_cuda *cuda_ci;
  struct cell_cuda *cuda_cj;
  
  //  time = 0;
  for (size_t i = 0; i < runs; ++i) {
    /* Zero the fields */
    zero_particle_fields(ci);
    //zero_particle_fields(cj);

    cuda_ci = copy_from_host(ci, 0);
    //cuda_cj = copy_from_host(cj, ci->count);
    //exit(1);
    //tic = getticks();

    /* Run the test */
    do_test<<<volume/CUDA_THREADS + 1,CUDA_THREADS>>>(cuda_ci, cuda_cj);
    cudaErrCheck( cudaPeekAtLastError() );
    cudaErrCheck( cudaDeviceSynchronize() );

    //toc = getticks();
    //time += toc - tic;

    copy_to_host(cuda_ci, ci);
    //copy_to_host(cuda_cj, cj);
    /* Dump if necessary */
    if (i % 50 == 0) {
      sprintf(outputFileName, "swift_gpu_%s.dat", outputFileNameExtension);
      dump_particle_fields(outputFileName, ci, cj);
    }
  }

  printf("rho %g\n", ci->parts[0].rho);
  /* Output timing */
//  message("SWIFT calculation took       %lli ticks.", time / runs);

  /* Now perform a brute-force version for accuracy tests */

  /* Zero the fields */
  //zero_particle_fields(ci);
  //zero_particle_fields(cj);

//  tic = getticks();

  /* Run the brute-force test */
  //pairs_all_density(&runner, ci, cj);

//  toc = getticks();

  /* Dump */
  sprintf(outputFileName, "swift_gpu_end_%s.dat", outputFileNameExtension);
  dump_particle_fields(outputFileName, ci, cj);

  /* Output timing */
//  message("Brute force calculation took %lli ticks.", toc - tic);

  /* Clean things to make the sanitizer happy ... */
  clean_up(ci);
  clean_up(cj);

  return 0;
}
