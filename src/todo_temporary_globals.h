#ifndef TODO_TEMPORARY_GLOBALS_H_
#define TODO_TEMPORARY_GLOBALS_H_

#define MLADENASSN 200 /*assume max number of neighbours*/




/* ================================ */
struct mladen_globals {
/* ================================ */

  FILE *outfilep;
  FILE *oneTimeFlagFilep;
  int called_fluxes;
  long long npart;

  int dump_nr;              /* index of dump */

  struct gizmo_debug_dump * data; /* array for data dump */

  struct engine* e;

} ;




extern struct mladen_globals mladen_globs;





/* ================================ */
struct gizmo_debug_dump {
/* ================================ */

  long long id;                     /* particle ID */
  float h;                          /* smoothing length */

  float volume_store;               /* particle volume */
  float omega;                      /* normalization for psi */
  float wgrads_store[3];            /* store sum of individual cartesian gradients contributions */
  float pos[3];                     /* particle positions */

  int nneigh;                             /* number of neighbours this particle interacts with */
  long long neighbour_ids[MLADENASSN];    /* IDs of each neighbour this particle interacts with */
  float Aij[2*MLADENASSN];                /* effective surface towards each neighbour */

  int nneigh_grads;
  long long neighbour_ids_grad[MLADENASSN];     /* IDs of neighbour for individual gradient contributions */
  float grads_sum_contrib[2*MLADENASSN];        /* contributions to the gradient sum from each neighbour */
  float dwdr[MLADENASSN];                       /* radial derivative of the kernel */
  float wjxi[MLADENASSN];                       /* Wj(xi) */
  float grads_final[2*MLADENASSN];              /* Actual analytical gradient */

  float grads_sum_dx[2*MLADENASSN];        /* pi.x - pj.x*/
  float r[MLADENASSN];                     /* |pi.x - pj.x | */

};



/* general functions */
void mladen_setup(struct engine* e);
void mladen_cleanup(void);


/* Aij debug dumps */
void mladen_dump_after_timestep(void);
void mladen_setup_data_dump(long long npart);
void mladen_reset_dump_data(void);
void mladen_store_particle_data(struct part *p, float h);

void mladen_store_neighbour_data(struct part *restrict pi,
    long long pjid, float wi, float GSCX, float GSCY, float GSDX, float GSDY,
    float dwdr, const float r, float hi);

void mladen_store_density_data(struct part *restrict pi,
    float hi, float Vi);

void mladen_store_Aij(struct part *restrict pi, struct part *restrict pj, float r, float hi,
    float* A, float grad_final_x, float grad_final_y, int negative);

void mladen_track_volume(const struct part *restrict pi, const struct part *restrict pj);



/* particle tracking */

void mladen_track_particle_stdout(struct part* restrict pi, int condition);


#endif /* todo_temporary_globals.h */
