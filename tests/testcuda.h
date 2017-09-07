typedef long long integertime_t;
typedef char timebin_t;

extern "C" {
#include <stdio.h>
#include <float.h>
}

#define part_align 128
#define CUDA_MAX_LINKS 27
#define CUDA_THREADS 256
#define DIM 3
#define NBRE_DIR 13
#define atomic_or(v, i) __sync_fetch_and_or(v, i)

#define kernel_degree 3
#define kernel_ivals 2  /*!< Number of branches */
#define kernel_gamma ((float)(1.825742))
#define kernel_gamma2 ((float)(kernel_gamma * kernel_gamma))
#define kernel_constant ((float)(16. * M_1_PI))
#define kernel_gamma_inv ((float)(1. / kernel_gamma))
#define kernel_gamma_dim ((float)(kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_ivals_f ((float)(kernel_ivals))
#define hydro_dimension 3.f
#define kernel_gamma_dim_plus_one \
  ((float)(kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_inv_dim \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma)))
#define kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma)))

__constant__ float cuda_kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)];

static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

__device__ __inline__ void cuda_kernel_deval(float u, float *restrict W,
                                             float *restrict dW_dx) {
  /* Go to the range [0,1] from [0,H] */
  const float x = u * kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *coeffs = &cuda_kernel_coeffs[ind * (kernel_degree + 1)];

  /* First two terms of the polynomial*/
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];
  /* and the rest of them*/
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  /* Return the results */
  *W = w * kernel_constant * kernel_gamma_inv_dim;
  *dW_dx = dw_dx * kernel_constant * kernel_gamma_inv_dim_plus_one;
}



#define error(s, ...)                                                      \
  ({                                                                       \
    fprintf(stderr, "%s:%s():%i: " s "\n", \
            __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);              \
    abort();                                                               \
  })

#define message(s, ...)                                                 \
  ({                                                                    \
    printf("%s: " s "\n", __FUNCTION__, \
           ##__VA_ARGS__);                                              \
  })


/**
 * Returns a random number (uniformly distributed) in [a,b[
 */
double random_uniform(double a, double b) {
  return (rand() / (double)RAND_MAX) * (b - a) + a;
}

/**
 * @brief The default struct alignment in SWIFT.
 */
#define SWIFT_STRUCT_ALIGNMENT 32

/**
 * @brief Defines alignment of structures
 */
#define SWIFT_STRUCT_ALIGN __attribute__((aligned(SWIFT_STRUCT_ALIGNMENT)))


/* Data of a single particle. */
struct part {

  /* Particle ID. */
  long long id;

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle mass. */
  float mass;

  /* Particle density. */
  float rho;

  /* Particle entropy. */
  float entropy;

  /* Entropy time derivative */
  float entropy_dt;

  //union {

    struct {

      /* Number of neighbours. */
      float wcount;

      /* Number of neighbours spatial derivative. */
      float wcount_dh;

      /* Derivative of the density with respect to h. */
      float rho_dh;

      /* Particle velocity curl. */
      float rot_v[3];

      /* Particle velocity divergence. */
      float div_v;

    } density;

    struct {

      /* Balsara switch */
      float balsara;

      /*! "Grad h" term */
      float f;

      /* Pressure over density squared  */
      float P_over_rho2;

      /* Particle sound speed. */
      float soundspeed;

      /* Signal velocity. */
      float v_sig;

      /* Time derivative of the smoothing length */
      float h_dt;

    } force;
  //};

  /* Time-step length */
  timebin_t time_bin;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

/* sorting stuff */
struct entry {

  /*! Distance on the axis */
  float d;

  /*! Particle index */
  int i;
};

/**
 * @brief Cell within the tree structure.
 *
 * Contains particles, links to tasks, a multipole object and counters.
 */
struct cell {

  /*! The cell location on the grid. */
  double loc[3];

  /*! The cell dimensions. */
  double width[3];

  /*! Max smoothing length in this cell. */
  double h_max;

  /*! Linking pointer for "memory management". */
  //struct cell *next;

  /*! Pointer to the #part data. */
  struct part *parts;

  /*! Pointer for the sorted indices. */
  struct entry *sort[NBRE_DIR];

  /*! Pointers to the next level of cells. */
  struct cell *progeny[8];

  /*! Parent cell. */
  struct cell *parent;

  /*! Super cell, i.e. the highest-level parent cell that has pair/self tasks */
  struct cell *super;

  /*! The task computing this cell's sorts. */
  //struct task *sorts;

  /*! Linked list of the tasks computing this cell's hydro density. */
  struct link *density;

  /* Linked list of the tasks computing this cell's hydro gradients. */
  //struct link *gradient;

  /*! Linked list of the tasks computing this cell's hydro forces. */
  //struct link *force;

  /*! Linked list of the tasks computing this cell's gravity forces. */
  //struct link *grav;

    /*! The ghost tasks */
  //struct task *ghost_in;
  //struct task *ghost_out;
  //struct task *ghost;

  /*! The extra ghost task for complex hydro schemes */
  //struct task *extra_ghost;

  /*! The drift task for parts */
  //struct task *drift_part;
  
  /*! The drift task for gparts */
  //struct task *drift_gpart;

  /*! The first kick task */
  //struct task *kick1;

  /*! The second kick task */
  //struct task *kick2;

  /*! The task to compute time-steps */
  //struct task *timestep;

  /*! Task linking the FFT mesh to the rest of gravity tasks */
  //struct task *grav_ghost[2];

  /*! Task computing long range non-periodic gravity interactions */
  //struct task *grav_long_range;

  /*! Task propagating the multipole to the particles */
  //struct task *grav_down;

  /*! Task for cooling */
  //struct task *cooling;

  /*! Task for source terms */
  //struct task *sourceterms;

  /*! Minimum end of (integer) time step in this cell. */
  integertime_t ti_end_min;

  /*! Maximum end of (integer) time step in this cell. */
  integertime_t ti_end_max;

  /*! Maximum beginning of (integer) time step in this cell. */
  //integertime_t ti_beg_max;

  /*! Last (integer) time the cell's part were drifted forward in time. */
  integertime_t ti_old_part;

  /*! Last (integer) time the cell's gpart were drifted forward in time. */
  //integertime_t ti_old_gpart;

  /*! Last (integer) time the cell's multipole was drifted forward in time. */
  //integertime_t ti_old_multipole;

  /*! Minimum dimension, i.e. smallest edge of this cell (min(width)). */
  float dmin;

  /*! Maximum particle movement in this cell since the last sort. */
  float dx_max_sort;

  /*! Maximum part movement in this cell since last construction. */
  float dx_max_part;

  /*! Maximum gpart movement in this cell since last construction. */
  float dx_max_gpart;

  /*! Nr of #part in this cell. */
  int count;

  /*! Nr of #gpart in this cell. */
  //int gcount;

  /*! Nr of #spart in this cell. */
  //int scount;

  /*! Bit-mask indicating the sorted directions */
  unsigned int sorted;

  /*! Spin lock for various uses (#part case). */
  //swift_lock_type lock;

  /*! Spin lock for various uses (#gpart case). */
  //swift_lock_type glock;

  /*! Spin lock for various uses (#multipole case). */
  //swift_lock_type mlock;

  /*! Spin lock for various uses (#spart case). */
  //swift_lock_type slock;

  /*! ID of the previous owner, e.g. runner. */
  //int owner;

  /*! Number of #part updated in this cell. */
  int updated;

  /*! Number of #gpart updated in this cell. */
  //int g_updated;

  /*! Number of #spart updated in this cell. */
  //int s_updated;

  /*! ID of the node this cell lives on. */
  int nodeID;

  /*! Is the #part data of this cell being used in a sub-cell? */
  //int hold;

  /*! Is the #gpart data of this cell being used in a sub-cell? */
  //int ghold;

  /*! Is the #multipole data of this cell being used in a sub-cell? */
  //int mhold;

  /*! Is the #spart data of this cell being used in a sub-cell? */
  //int shold;

  /*! Number of tasks that are associated with this cell. */
  //short int nr_tasks;

  /*! The depth of this cell in the tree. */
  //char depth;

  /*! Is this cell split ? */
  char split;

  /*! The maximal depth of this cell and its progenies */
  //char maxdepth;

  /*! Values of dx_max and h_max before the drifts, used for sub-cell tasks. */
  //float dx_max_old;
  //float h_max_old;
  //float dx_max_sort_old;

  /* Bit mask of sort directions that will be needed in the next timestep. */
  //unsigned int requires_sorts;

  /*! Does this cell need to be drifted? */
  //char do_drift;

  /*! Do any of this cell's sub-cells need to be drifted? */
  //char do_sub_drift;

  /*! Bit mask of sorts that need to be computed for this cell. */
  //unsigned int do_sort;

  /*! Do any of this cell's sub-cells need to be sorted? */
  //char do_sub_sort;

  /* Gives the cells a unique indentifier (slash index). */
  int cuda_ID;

  /* ID of the cells' load and unload task used for dependencies. */
  //int load_task, unload_task;

} SWIFT_STRUCT_ALIGN;


struct cell_cuda {
  /* The cell location on the grid. */
  double loc[3];

  /* The cell dimensions. */
  double width[3];

  /* Max smoothing length in this cell. */
  double h_max;

  /* Index of the particle data. */
  int first_part;
  
  /* Number of particles in the cell. */
  int part_count;

  /* Indices of the next level of cells. */
  int progeny[8];

  /* Index of the parent cell. */
  int parent;

  /* Index of the super cell.*/
  int super;

  /* Minimum end of time step in this cell. */
  integertime_t ti_end_min;

  /* Maximum end of time step in this cell. */
  integertime_t ti_end_max;
  
  /* Minimum dimension of this cell */
  float dmin;

  /* Need the density links only, maximum of 27 */
  int links[CUDA_MAX_LINKS];

  /* Number of links */  
  int nr_links;

  /* IS split? */
  int split;

};


/*Definition of particle structures (SoA on GPU) */
struct particle_arrays {
  /* Device array containing the particle ID */
  long long int *id;
  /*Device array containing the particle positions. */
  double *x_x, *x_y, *x_z;
  /*Device array containing the particle predicted velocity */
  float3 *v;
  /*Device array containing the particle acceleration. */
  float3 *a_hydro;
  /* Device array contining the particle cutoff radius */
  float *h;
  /* Device array containing particle masses */
  float *mass;
  /* Device array containing particle densities. */
  float *rho;
  /* Device array containing the particle entropies. */
  float *entropy;
  /* Device array containing the particle entropy_dt */
  float *entropy_dt;

  /* Device array containing sort */
  struct entry *sort[NBRE_DIR];
  /* Unclear if unions work on the GPU, so unrolling the union for now */
  /* DENSITY */
  /* Device array containing the number of neighbours. */
  float *wcount;
  /* Device array containing the number of neighbours spatial derivative */
  float *wcount_dh;
  /* Device array containing the derative of the density w.r.t. h */
  float *rho_dh;
  /* Device array containing the particle velocity curl. */
  float3 *rot_v;
  /* Device array containing the particle velocity divergence */
  float *div_v;

  /* FLOAT */
  /* Device array containing the balsara switch */
  float *balsara;
  /* Device array containing the "Grad h" term */
  float *f;
  /* Device array containing the pressure over density squared */
  float *P_over_rho2;
  /* Device array containing the particle sound speeds */
  float *soundspeed;
  /* Device array containing the signal velocities */
  volatile float *v_sig;
  /* Device array containing the time derivative of the smoothing lengths */
  float *h_dt;

  /* Device array containing time step length */
  timebin_t *time_bin;
};

__device__ struct particle_arrays cuda_parts;


__device__ __constant__ int ti_current = 8;

__device__ __constant__ int max_active_bin = 1;

__device__ __constant__ double dim[3];

__device__ __inline__ int cuda_cell_is_active(struct cell_cuda *c) {
  return (c->ti_end_min == ti_current);
}

__device__ __inline__ int cuda_part_is_active(int pid) {

  return (cuda_parts.time_bin[pid] <= max_active_bin);
}


__device__ __constant__ float const_viscosity_alpha=0.8;

#define HYDRO_DIMENSION_3D

__device__ float cuda_pow_dimension_plus_one(float x) {
#if defined(HYDRO_DIMENSION_3D)
  const float x2 = x * x;
  return x2 * x2;
#elif defined(HYDRO_DIMENSION_2D)
  return x * x * x;
#elif defined(HYDRO_DIMENSION_1D)
  return x * x;
#else
  printf("The dimension is not defined!");
  return 0.f;
#endif
}



/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */
void runner_do_sort_ascending(struct entry *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[10];
  int qpos, i, j, lo, hi, imin;
  struct entry temp;
  float pivot;

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    if (hi - lo < 15) {
      for (i = lo; i < hi; i++) {
        imin = i;
        for (j = i + 1; j <= hi; j++)
          if (sort[j].d < sort[imin].d) imin = j;
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      while (i <= j) {
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        if (i <= j) {
          if (i < j) {
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

/* Orientation of the cell pairs */
static const double runner_shift[13][3] = {
    {5.773502691896258e-01, 5.773502691896258e-01, 5.773502691896258e-01},
    {7.071067811865475e-01, 7.071067811865475e-01, 0.0},
    {5.773502691896258e-01, 5.773502691896258e-01, -5.773502691896258e-01},
    {7.071067811865475e-01, 0.0, 7.071067811865475e-01},
    {1.0, 0.0, 0.0},
    {7.071067811865475e-01, 0.0, -7.071067811865475e-01},
    {5.773502691896258e-01, -5.773502691896258e-01, 5.773502691896258e-01},
    {7.071067811865475e-01, -7.071067811865475e-01, 0.0},
    {5.773502691896258e-01, -5.773502691896258e-01, -5.773502691896258e-01},
    {0.0, 7.071067811865475e-01, 7.071067811865475e-01},
    {0.0, 1.0, 0.0},
    {0.0, 7.071067811865475e-01, -7.071067811865475e-01},
    {0.0, 0.0, 1.0},
};

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param c The #cell.
 */
void do_sort(struct cell *c) {

  struct part *parts = c->parts;
  const int count = c->count;

  /* start by allocating the entry arrays in the requested dimensions. */
  for (int j = 0; j < NBRE_DIR; j++) {
    if (c->sort[j] == NULL) {
      if ((c->sort[j] = (struct entry *)malloc(sizeof(struct entry) *
                                               (count + 1))) == NULL)
        error("Failed to allocate sort memory.");
    }
  }
  
  /* Fill the sort array. */

  for (int k = 0; k < count; k++) {
    const double px[3] = {parts[k].x[0], parts[k].x[1], parts[k].x[2]};
    for (int j = 0; j < NBRE_DIR; j++)
      {
	c->sort[j][k].i = k;
	c->sort[j][k].d = px[0] * runner_shift[j][0] +
	  px[1] * runner_shift[j][1] +
	  px[2] * runner_shift[j][2];
      }
  }

  /* Add the sentinel and sort. */
  for (int j = 0; j < NBRE_DIR; j++) {
      c->sort[j][count].d = FLT_MAX;
      c->sort[j][count].i = 0;
      runner_do_sort_ascending(c->sort[j], count);
      atomic_or(&c->sorted, 1 << j);
    }
  
}
