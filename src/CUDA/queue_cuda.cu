#ifdef WITH_CUDA
#ifndef static
#define static
#endif
#ifndef restrict
#define restrict __restrict__
#endif
#endif


#include "queue_cuda.h"

/* Queue function to search for a task index from a specific queue. */
__device__ int cuda_queue_gettask( struct queue_cuda *q ) {

  int ind, tid = -1;

  /* Don't try if queue is empty. */
  if( q->rec_count == q->count )
    return -1;

  /* Get index of the next task. */
  ind = atomicAdd(&q->first, 1);

  ind %= cuda_queue_size;
  /* Loop until there is a valid task at that index. */
  while( q->rec_count < q->count && (tid = q->data[ind] ) < 0);

  /* Remove the task from the queue. */
  if( tid >= 0 )
  {
    q->data[ind] = -1;
    atomicAdd((int*) &tot_num_tasks, -1);
  }

  /* Return the acquired task ID */
  return tid;
}

/* Queue function to add task index tid to the queue*/
__device__ void cuda_queue_puttask( struct queue_cuda *q, int tid ) {

     int ind;

    /* Get the index of the next task. */
    ind = atomicAdd( &q->last , 1 ) % cuda_queue_size;

    /* Wait for the slot in the queue to be empty. */
    while ( q->data[ind] != -1 );

    /* Write the task back to the queue. */
    q->data[ind] = tid;

    atomicAdd((int*)&q->nr_avail_tasks, 1);
}

#ifdef WITH_CUDA
#undef static
#undef restrict
#endif
