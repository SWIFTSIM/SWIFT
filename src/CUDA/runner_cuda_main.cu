
#include "queue_cuda.h"
#include "task_cuda.h"

/*Definition of particle structures (SoA on GPU) */
struct particle_arrays{
/* Device array containing the particle ID */
__device__ long long int *id;
/*Device array containing the particle positions. */
__device__ double *x_x, *x_y, *x_z;
/*Device array containing the particle predicted velocity */
__device__ float3 *v;
/*Device array containing the particle acceleration. */
__device__ float3 *a_hydro;
/* Device array contining the particle cutoff radius */
__device__ float *h;
/* Device array containing particle masses */
__device__ float *mass;
/* Device array containing particle densities. */
__device__ float *rho;
/* Device array containing the particle entropies. */
__device__ float *entropy;

/* Unclear if unions work on the GPU, so unrolling the union for now */
/* DENSITY */
/* Device array containing the number of neighbours. */
__device__ float *wcount;
/* Device array containing the number of neighbours spatial derivative */
__device__ float *wcount_dh;
/* Device array containing the derative of the density w.r.t. h */
__device__ float *rho_dh;
/* Device array containing the particle velocity curl. */
__device__ float3 *rot_v;
/* Device array containing the particle velocity divergence */
__device__ float *div_v;

/* FLOAT */
/* Device array containing the balsara switch */
__device__ float *balsara;
/* Device array containing the "Grad h" term */
__device__ float f;
/* Device array containing the pressure over density squared */
__device__ float *P_over_rho2;
/* Device array containing the particle sound speeds */
__device__ float *soundspeed;
/* Device array containing the signal velocities */
__device__ float *v_sig;
/* Device array containing the time derivative of the smoothing lengths */
__device__ float *h_dt;

/* Device array containing time step length */
__device__ timebin_t *time_bin;
} cuda_parts;

/* Cell array and array of CPU cell pointers.*/
__device__ struct cell_cuda *cells;
__device__ struct cell *cpu_cells;

/* Queue variables*/
__device__ struct queue_cuda cuda_queues[ cuda_numqueues ];
__device__ struct load_queue;
__device__ struct unload_queue;

/* Array of tasks on the device */
__device__ struct task_cuda *tasks;



__device__ void load_cell(int cell_index)
{
  /* Get the pointer to the relevant cells. */
  struct cell *cpu_cell = &cpu_cells[cell_index];
  struct cell *cell = &cells[cell_index];
  struct part *parts = cpu_cell->parts;
  int i;
  /* Get the index to copy the data for the 0th particle in this cell to.*/
  int start = cell->first_part;

  for( i = threadIdx.x; i < cell->part_count; i+= blockDim.x ){
    struct part *current = &parts[i];
    /* Load the ID and position*/
    long long local_id = current->id;
    double local_x[3];
    local_x[0] = current->x[0];
    local_x[1] = current->x[1];
    local_x[2] = current->x[2];

    /* Unpack the ID and position and put into the device arrays. */
    cuda_parts.id[start+i] = local_id;
    cuda_parts.x_x[start+i] = local_x[0];
    cuda_parts.x_y[start+i] = local_x[1];
    cuda_parts.x_z[start+i] = local_x[2];

    /* Load the predicted velocity */ 
    float local_tempv[3];
    local_tempv[0] = current->v[0];
    local_tempv[1] = current->v[1];
    local_tempv[2] = current->v[2];

    /* Repack the predicted velocity */
    float3 local_v;
    local_v.x = local_tempv[0];
    local_v.y = local_tempv[1];
    local_v.z = local_tempv[2];

    /* Copy the predicted velocity to the particle array. */
    cuda_parts.v[start+i] = local_v;
   
    /* Load the acceleration */ 
    float local_tempa[3];
    local_tempa[0] = current->a_hydro[0];
    local_tempa[1] = current->a_hydro[1];
    local_tempa[2] = current->a_hydro[2];

    /* Repack the acceleration */
    float3 local_a;
    local_a.x = local_tempa[0];
    local_a.y = local_tempa[1];
    local_a.z = local_tempa[2];

    /* Copy the predicted velocity to the particle array. */
    cuda_parts.a_hydro[start+i] = local_a;

    /* Copy the cutoff, mass, density, entropy and entropy derivative*/
    float local_h, local_mass, local_rho, local_entropy, local_entropy_dt;
    local_h = current->h;
    local_mass = current->mass;
    local_rho = current->rho;
    local_entropy = current->entropy;
    local_entropy_dt = current->entropy_dt;

    cuda_parts.h[start+i] = local_h;
    cuda_parts.mass[start+i] = local_mass;
    cuda_parts.rho[start+i] = local_rho;
    cuda_parts.entropy[start+i] = local_entropy;
    cuda_parts.entropy_dt[start+i] = local_entropy_dt;

    /* Copy the density union to the device */
    float local_wcount, local_wcount_dh, local_rho_dh, local_div_v;
    float3 local_rot_v;
    local_wcount = current->density.wcount;
    local_wcount_dh = current->density.wcount_dh;
    local_rho_dh = current->density.rho_dh;
    local_rot_v.x = current->density.rot_v[0];
    local_rot_v.y = current->density.rot_v[1];
    local_rot_v.z = current->density.rot_v[2];
    local_div_v = current->density.div_v;

    cuda_parts.wcount[start+i] = local_wcount;
    cuda_parts.wcount_dh[start+i] = local_wcount_dh;
    cuda_parts.rho_dh[start+i] = local_rho_dh;
    cuda_parts.rot_v[start+i] = local_rot_v;
    cuda_parts.div_v = local_div_v;

    /* Copy the force union to the device. TODO Do we need this?*/
    float local_balsara, local_f, local_P_over_rho2, local_soundspeed, local_vsig, local_h_dt;

    local_balsara = current->force.balsara;
    local_f = current->force.f;
    local_P_over_rho2 = current->force.P_over_rho2;
    local_soundspeed = current->force.soundspeed;
    local_v_sig = current->force.v_sig;
    local_h_dt = current->force.h_dt;

    cuda_parts.balsara[start+i] = local_balsara;
    cuda_parts.f[start+i] = local_f;
    cuda_parts.P_over_rho2[start+i] = local_P_over_rho2;
    cuda_parts.soundspeed[start+i] = local_soundspeed;
    cuda_parts.v_sig[start+i] = local_v_sig;
    cuda_parts.h_dt[start+i] = local_h_dt;

    cuda_parts.time_bin[start+i] = current->time_bin;
  }
}

/* Task function to unload a specific cell. */
/* TODO: Note that this copies back certain data that I expect is not modified.*/
__device__ void unload_cell( int cell_index ) {

  /* Get the pointer to the relevant cells. */
  struct cell *cpu_cell = &cpu_cells[cell_index];
  struct cell *cell = &cells[cell_index];
  struct part *parts = cpu_cell->parts;
  int i;
  /* Get the index to copy the data for the 0th particle in this cell to.*/
  int start = cell->first_part;

  for(i = threadIdx.x; i < cell->part_count; i+= blockDim.x )
  { 
    struct part *current = &parts[i];

    /* Copy back the ID and position.*/
    current->id = cuda_parts.id[start+i]
    current->x[0] = cuda_parts.x_x[start+i];
    current->x[1] = cuda_parts.x_y[start+i];
    current->x[2] = cuda_parts.x_z[start+i];

    /* Copy back the velocity*/
    float3 local_v = cuda_parts.v[start+i];
    current->v[0] = local_v.x;
    current->v[1] = local_v.y;
    current->v[2] = local_v.z;

    /* Copy back the acceleration */
    float3 local_a_hydro = cuda_parts.a_hydro[start+i];
    current->a_hydro[0] = local_a_hydro.x;
    current->a_hydro[1] = local_a_hydro.y;
    current->a_hydro[2] = local_a_hydro.z;

    /* Copy back the cutoff, mass, density, entropy and entropy_dt*/
    current->h = cuda_parts.h[start+i];
    current->mass = cuda_parts.mass[start+i];
    current->rho = cuda_parts.rho[start+i];
    current->entropy = cuda_parts.entropy[start+i];
    current->entropy_dt = cuda_parts.entropy_dt[start+i];

    /* Copy back the force union.*/
    current->force.balsara = cuda_parts.balsara[start+i];
    current->force.f = cuda_parts.f[start+i];
    current->force.P_over_rho2 = cuda_parts.P_over_rho2[start+i];
    current->force.soundspeed = cuda_parts.soundspeed[start+i];
    current->force.v_sig = cuda_parts.v_sig[start+i];
    current->force.h_dt = cuda_parts.h_dt[start+i];

    /* Copy back the timebin. */
    current->time_bin = cuda_parts.time_bin[start+i];
  }

}

/* Runner function to retrieve a task index from a queue. */
__device__ int runner_cuda_gettask( struct queue_cuda *q ) {

  int tid = -1;
  if( atomicAdd((int*)&q->nr_avail_tasks, -1) <= 0)
  {
    atomicAdd((int*)&q->nr_avail_tasks, 1);
    return -1;
  }

  /* Main loop */
  while( ( tid = cuda_queue_gettask( q ) ) >= 0 ){
    //TODO Do we need to lock anything here? Probably no.
    //TODO Does this need to be a while?
    break;
  }
  if( tid >= 0 ) {
    q->rec_data[atomicAdd( (int *)&q->rec_count, 1 ) ] = tid; //TODO Are we keeping rec_data
  }

  return tid;
}


/* The main kernel. */
__global__ void swift_device_kernel( )
{
  __shared__ volatile int tid;
  __shared__ volatile int done;
  int *src, *dest;
  int i;

  /* Main loop */
  while( 1 ){
    __syncthreads();
    /* Get a task from the queue. */
    if(threadIdx.x == 0){
      tid = -1;
      /* Highest Priority queue, unload tasks. */
      if(unload_queue.nr_avail_tasks > 0 ) {
        tid = runner_cuda_gettask( &unload_queue);
      }

      /* Next highest priority queue, load tasks. Only some blocks look in here */
      if(tid < 0 && load_queue.nr_avail_tasks > 0 && blockIdx.x < cuda_numloaders ) {
        tid = runner_cuda_gettask( &load_queue );
      }

      /* Finally loop through work queues in priority order. queue 0 is highest priority*/
      for( i = 0; i < cuda_numqueues && tid < 0; i++ ) {
        if( cuda_queues[i].nr_avail_tasks > 0 ){
          tid = runner_cuda_gettask( &cuda_queues[i] );
        }
      }

    }//Get task from queue

    /* Threads need to wait until they have work to do */
    __syncthreads();

    /* If the tasks are all complete and we don't have anything to do exit*/
    if(tid < 0 && tot_num_tasks == 0)
      break;

    if( tid < 0 )
      continue;

    //TODO Do work here!


  } // End of main loop

}

