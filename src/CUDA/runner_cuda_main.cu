
#include "queue_cuda.h"
#include "task_cuda.h"

/*Definition of particle structures (SoA on GPU) */

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

/* Cell array and array of CPU cell pointers.*/
__device__ struct cell_cuda *cells;
__device__ struct cell *cpu_cells;

/* Queue variables*/
__device__ struct queue_cuda cuda_queues[ cuda_numqueues ];
__device__ struct load_queue;
__device__ struct unload_queue;

/* Array of tasks on the device */
__device__ struct task_cuda *tasks;



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

