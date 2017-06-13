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

/* The number of priority queues. Only 2 is supported at the moment. */
#define cuda_numqueues 2

/* The number of blocks to perform load tasks. */
#define cuda_numloaders 12

/* Task types for reserved tasks */
#define type_load -101
#define type_unload -102
#define type_implicit -103

/* Queue constant data required to be set on the CPU */
__constant__ int cuda_queue_size;
__constant__ int cuda_nrqueues;
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


int cuda_queue_gettask( struct queue_cuda *q)

#endif /* SWIFT_CUDA_QUEUE_H */
