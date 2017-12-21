/*******************************************************************************
 * This file is part of QuickSched.
 * Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk), Aidan Chalk (aidan.chalk@durham.ac.uk)
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
#ifndef SWIFT_CUDA_QUEUE_H
#define SWIFT_CUDA_QUEUE_H

extern "C" {
#include "task.h"
}

/* The number of priority queues. Only 2 is supported at the moment. */
#define cuda_numqueues 2

/* The number of blocks to perform load tasks. */
#define cuda_numloaders 12

/* GPU types, references the types declared in task.h.*/
#define type_load task_type_load
#define type_unload task_type_unload
#define type_implicit_load task_type_implicit_load
#define type_implicit_unload task_type_implicit_unload
#define type_recv_load task_type_recv_load
#define type_send_unload task_type_send_unload

const int num_gpu_types = 5;
const int gpu_work_task_array[num_gpu_types] = {task_type_self, task_type_pair, task_type_sub_self,
                                                task_type_sub_pair, task_type_ghost};

/* Queue constant data required to be set on the CPU */
__constant__ int cuda_queue_size;
__constant__ int cuda_numtasks;
__constant__ int median_cost; //Used for priority queues.
__device__ int tot_num_tasks;

/* Barrier variable to exit nicely */
__device__ int cuda_barrier = 0;

/** Struct for a task queue. */
struct queue_cuda {

    /* Indices to the first and last elements. */
    int first, last;
    
    /* Number of elements in this queue. */
    volatile int count;
    
    /* Number of elements in the recycled list. */
    volatile int rec_count;
    
    /* The queue data. */
    volatile int *data;
    
    /* The recycling list. */
    volatile int *rec_data;
    
    volatile int nr_avail_tasks;

    };


__device__ int cuda_queue_gettask( struct queue_cuda *q);
__device__ void cuda_queue_puttask( struct queue_cuda *q, int tid );

#endif /* SWIFT_CUDA_QUEUE_H */
