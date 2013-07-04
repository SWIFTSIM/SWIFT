/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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


/* Some constants. */
#define queue_maxsuper           50
#define queue_sizeinit           100
#define queue_sizegrow           2


/* Counters. */
enum {
    queue_counter_swap = 0,
    queue_counter_count,
    };
extern int queue_counter[ queue_counter_count ];


/** The queue struct. */
struct queue {

    /* The lock to access this queue. */
    lock_type lock;

    /* Size, count and next element. */
    int size, count;
    
    /* The actual tasks to which the indices refer. */
    struct task *tasks;
    
    /* The task indices. */
    int *tid;

    } __attribute__((aligned (64)));
    

/* Function prototypes. */
struct task *queue_gettask ( struct queue *q , int qid , struct cell *super , int blocking );
void queue_init ( struct queue *q , struct task *tasks );
void queue_insert ( struct queue *q , struct task *t );
