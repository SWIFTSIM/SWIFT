#include <stdio.h>
#include <string.h>
#include "hydro.h"
#include "error.h"
#include "engine.h"
#include "todo_temporary_globals.h"



#define DO_MLADEN_DEBUG_OUTPUT
/* #define DO_MLADEN_TRACKING */

#define MAXDEBUGDUMPS 10


struct mladen_globals mladen_globs;

/* ========================================= */
void mladen_setup(struct engine* e){
/* ========================================= */
  /* sets up the necessary stuff. To be called in main.
   * Just drop in everything you need to be done first
   * in here.
   * written to be called at the beginning of main.c
   */
#ifdef DO_MLADEN_DEBUG_OUTPUT

  /* mladen_globs.outfilep = fopen("mladen_outputfile_all.txt", "w");   */
  /* mladen_globs.oneTimeFlagFilep = fopen("mladen_flags_all.txt", "w");   */

  mladen_globs.called_fluxes = 1;
  mladen_globs.e = e;

  message("MESSAGE TO MLADEN: calling init your own global temporary stuff.");


  mladen_globs.dump_nr = 0;

#endif


}






/* ========================================= */
void mladen_cleanup(void){
/* ========================================= */
  /* clean up after yourself */
  /* to be called at the end of the run in main.c. */

#ifdef DO_MLADEN_DEBUG_OUTPUT

  /* fclose(mladen_globs.outfilep); */
  /* fclose(mladen_globs.oneTimeFlagFilep); */

#endif

}








/* ========================================= */
void mladen_setup_data_dump(long long npart){
/* ========================================= */

  /* allocate data dump array, set up necessary stuff.
   * call in main after you established all the necessary
   * parameters necessary.
   * Call somewhere in main.c before sim loop starts*/

#ifdef DO_MLADEN_DEBUG_OUTPUT
  mladen_globs.npart = npart+1;
  mladen_globs.data = (struct gizmo_debug_dump*) swift_malloc("debug_data", (npart+1)*sizeof(struct gizmo_debug_dump));

  /* initialize array to zero */
  mladen_reset_dump_data();

#endif


}








/* ======================================= */
void mladen_reset_dump_data(){
/* ======================================= */
  /* reset dump data values
   * called after dump is written, or when dump data is initialized, both in mladen_* functions*/
#ifdef DO_MLADEN_DEBUG_OUTPUT
  struct mladen_globals m = mladen_globs;
  for (int i = 0; i < m.npart; i++){
    struct gizmo_debug_dump *d = &m.data[i];
    d->nneigh = 0;
    d->nneigh_grads = 0;
    d->volume_store = 0;
    d->omega = 0;
    d->wgrads_store[0] = 0;
    d->wgrads_store[1] = 0;
    d->wgrads_store[2] = 0;
    d->pos[0] = 0;
    d->pos[1] = 0;
    d->pos[2] = 0;
    d->h = 0;
    d->id = 0;

    for (int j = 0; j < MLADENASSN; j++){
      d->neighbour_ids[j] = 0;
      d->neighbour_ids_grad[j] = 0;
      d->dwdr[j] = 0;
      d->wjxi[j] = 0;
      d->r[j] = 0;
    }
    for (int j = 0; j < 2*MLADENASSN; j++){
      d->grads_sum_contrib[j] = 0;
      d->Aij[j] = 0;
      d->grads_sum_dx[j] = 0;
      d->grads_final[j] = 0;
    }
  }

#endif

}







/* ========================================= */
void mladen_dump_after_timestep(void){
/* ========================================= */

  /* written to be called from engine_step in engine.c
   * right after
    if (!e->restarting)
      fprintf(
          e->file_timesteps,
          "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %12lld "
          "%21.3f %6d\n",
          e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
          e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
          e->s_updates, e->b_updates, e->wallclock_time, e->step_props);
    */

#ifdef DO_MLADEN_DEBUG_OUTPUT

#ifdef MAXDEBUGDUMPS
  if (mladen_globs.dump_nr == MAXDEBUGDUMPS) return;
#endif

  if (mladen_globs.e->updates==mladen_globs.npart-1){ /* npart = npart+1 in setup*/

    printf("Dumping Mladen debugging data %d after timestep\n", mladen_globs.dump_nr);

    /* get filename */
    char filename[200] = "swift-gizmo-debug-dump_";
    char dumpstring[5];
    sprintf(dumpstring, "%04d", mladen_globs.dump_nr);
    strcat(filename, dumpstring);
    strcat(filename, ".dat");


    FILE *fp;
    fp = fopen(filename, "wb");

    long long np = mladen_globs.npart;

    /* first dump nparts */
    fwrite(&np, sizeof(long long), 1, fp);
    int maxneigh_guess = MLADENASSN;
    fwrite(&maxneigh_guess, sizeof(int), 1, fp);

    /* now dump particle after particle */

    for (int p=1; p<np; p++){
      struct gizmo_debug_dump * d = &(mladen_globs.data[p]);
      fwrite(&(d->id), sizeof(long long), 1, fp);
      fwrite(&(d->h), sizeof(float), 1, fp);
      fwrite(&(d->omega), sizeof(float), 1, fp);
      fwrite(&(d->volume_store), sizeof(float), 1, fp);
      fwrite(&(d->wgrads_store), sizeof(float), 3, fp);
      fwrite(&(d->pos), sizeof(float), 3, fp);


      fwrite(&(d->nneigh_grads), sizeof(int), 1, fp);
      fwrite(&(d->neighbour_ids_grad), sizeof(long long), MLADENASSN, fp);
      fwrite(&(d->grads_sum_contrib), sizeof(float), 2*MLADENASSN, fp);
      fwrite(&(d->dwdr), sizeof(float), MLADENASSN, fp);
      fwrite(&(d->wjxi), sizeof(float), MLADENASSN, fp);
      fwrite(&(d->grads_sum_dx), sizeof(float), 2*MLADENASSN, fp);
      fwrite(&(d->r), sizeof(float), MLADENASSN, fp);

      fwrite(&(d->nneigh), sizeof(int), 1, fp);
      fwrite(&(d->neighbour_ids), sizeof(long long), MLADENASSN, fp);
      fwrite(&(d->Aij), sizeof(float), 2*MLADENASSN, fp);
      fwrite(&(d->grads_final), sizeof(float), 2*MLADENASSN, fp);


      /* to check whether you're writing the correct stuff: you need to have read in 'teststring' as well */
      /* char teststring[11] = "teststring"; */
      /* fwrite(&teststring, sizeof(char), 11, fp); */

    }

    fclose(fp);

    /* increase dump index */
    mladen_globs.dump_nr += 1;
  }

  /* reset values */
  mladen_reset_dump_data();

#endif

}







/* ======================================================== */
void mladen_store_particle_data(struct part *p, float h){
/* ======================================================== */

  /* store default particle data
   * Gets called whenever a mladen_store* routine is called
   * */

#ifdef DO_MLADEN_DEBUG_OUTPUT

  int ind = (int) p->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);
  dumploc->id = p->id;
  dumploc->pos[0] = p->x[0];
  dumploc->pos[1] = p->x[1];
  dumploc->pos[2] = p->x[2];
  dumploc->h = h;

#endif

}








/* ======================================================== */
void mladen_store_neighbour_data(struct part *restrict pi,
    long long pjid, float wi, float GSCX, float GSCY,
    float GSDX, float GSDY, float dwdr, const float r,
    float hi){
/* ======================================================== */
  /* pi: particle i
   * pjid: id of particle j
   * wi:  Wj(xi)
   * GSCX, GSCY: Gradient Sum Contribution in x, y direction
   * GSDX, GSDY: Gradient sum dx, dy
   * dwdr: dwdr
   * r: r
   * hi: h_i
   *
   *
   * written to be called from runner_iact_density and
   * runner_iact_nonsym_density
   */


  /* printf("   Storing data for particle %lld\n", pi->id); */
#ifdef DO_MLADEN_DEBUG_OUTPUT

#ifdef WITH_IVANOVA
  /* first store particle data to be safe */
  mladen_store_particle_data(pi, hi);

  /* quit if you're outside the compact support radius */
  if (r > hi*1.778002) return;

  int ind = (int) pi->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);

  int ng = dumploc->nneigh_grads;
  int new =  1;

  /* first check whether the neighbour is already in there */
  for (int i = 0; i<ng; i++){
    if (dumploc->neighbour_ids_grad[i] == pjid){
      ng = i;
      new = 0;
      break;
    }
  }

  /* if new, add neighbours */
  if (new){
    ng = dumploc->nneigh_grads;
    dumploc->nneigh_grads += 1;
  }

  /* check to not make segfaults */
  if (dumploc->nneigh_grads == MLADENASSN) error("Particle %lld has > %d neighbours\n", pi->id, MLADENASSN);

  /* store stuff */
  dumploc->neighbour_ids_grad[ng] = pjid;
  dumploc->grads_sum_contrib[2*ng] = GSCX;
  dumploc->grads_sum_contrib[2*ng+1] = GSCY;
  dumploc->grads_sum_dx[2*ng] = GSDX;
  dumploc->grads_sum_dx[2*ng+1] = GSDY;
  dumploc->dwdr[ng] = dwdr;
  dumploc->wjxi[ng] = wi;
  dumploc->r[ng] = r;
#endif
#endif

}








/* ======================================================== */
void mladen_store_density_data(struct part *restrict pi,
    float hi, float Vi){
/* ======================================================== */

  /* written to be called from runner_iact_fluxes_common */

#ifdef DO_MLADEN_DEBUG_OUTPUT
#ifdef WITH_IVANOVA
  mladen_store_particle_data(pi, hi);
  int ind = (int) pi->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);

  dumploc->wgrads_store[0] = pi->density.wgrads[0];
  dumploc->wgrads_store[1] = pi->density.wgrads[1];
  dumploc->wgrads_store[2] = pi->density.wgrads[2];

  dumploc->volume_store = Vi;
  dumploc->omega = pi->density.wcount;
#endif
#endif


}







/* ========================================================================================== */
void mladen_store_Aij(struct part *restrict pi, struct part *restrict pj, float r, float hi,
  float* A, float grad_final_x, float grad_final_y, int negative){
/* ========================================================================================== */

  /* written to be called from runner_iact_fluxes_common */
  /* grad_final_x: \del psi_j (x_i) / \del x: The gradient of psi which will
   * be multiplied by V_i in the A_ij */

#ifdef DO_MLADEN_DEBUG_OUTPUT
  mladen_store_particle_data(pi, hi);
#ifdef WITH_IVANOVA

  /* quit if you're outside the compact support radius */
  if (r > hi*1.778002) return;

  int ind = (int) pi->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);

  int n = dumploc->nneigh;
  int new =  1;

  /* first check whether the neighbour is already in there */
  for (int i = 0; i<n; i++){
    if (dumploc->neighbour_ids[i] == pj->id){
      n = i;
      new = 0;
      break;
    }
  }

  /* if new, add neighbours */
  if (new){
    n = dumploc->nneigh;
    dumploc->nneigh+= 1;
  }

  /* now store stuff */
  if (negative){
    dumploc->Aij[2*n] = -A[0];
    dumploc->Aij[2*n+1] = -A[1];
  }
  else{
    dumploc->Aij[2*n] = A[0];
    dumploc->Aij[2*n+1] = A[1];
  }
  dumploc->grads_final[2*n] = grad_final_x;
  dumploc->grads_final[2*n+1] = grad_final_y;
  dumploc->neighbour_ids[n] = pj->id;

#endif
#endif

}




/*============================================================================================================*/
void mladen_track_volume(const struct part *restrict pi, const struct part *restrict pj){
/*============================================================================================================*/
  /* Track what particle adds what contribution to the volume and print information */
  /* to be called at the start of runner_iact_density and runner_iact_nonsym density*/
  /* call it twice with pi, pj inversed in symmetric one                            */

#ifdef DO_MLADEN_DEBUG_OUTPUT
/* which particle to look for */
#define MLADEN_PART_I 2066
#ifdef WITH_IVANOVA

  if (pi->id != MLADEN_PART_I) return;


  float wi, wi_dx;

  double dx[3];

  for (int i = 0; i<3; i++){
    dx[i] = pi->x[i] - pj->x[i];
  }

  float r2 = (float) dx[0]*dx[0] + dx[1]*dx[1];
  float hi = pi->h;

  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  wi = wi * hi_inv * hi_inv;

  float gamma = 1.778002;

  printf("Adding to i=%6lld from j=%6lld \n", pi->id, pj->id);
  printf("wcount is now: %14.7e || adding %14.7e \n", pi->density.wcount*hi_inv*hi_inv, wi);
  printf("r:   %14.7e \nH:   %14.7e \nr/H: %14.7e \n", r, hi*gamma, r/hi/gamma);
  printf("\n");

#endif
#endif
}





/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------------------------------------------------ */


/* ========================================================================================== */
void mladen_track_particle_stdout(struct part* restrict pi, int condition){
/* ========================================================================================== */

  /* track particles that satisfy some condition given by condition int
   * 1: 
   */

#ifdef DO_MLADEN_TRACKING
  
  if (condition == 1){
    float r = sqrtf((float) pi->x[0]*pi->x[0] + (float) pi->x[1]*pi->x[1] + (float) pi->x[2]+pi->x[2]);
    float u = hydro_get_comoving_internal_energy(pi);
    if (r>=0.2 && u >=5){
      printf("Got u=%14.7e for r=%14.7f\n", u, r);
    }
  }
  


#endif

}
