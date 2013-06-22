/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#define scheduler_maxwait                    3
#define scheduler_maxunlock                  40
#define scheduler_dosub                      1
#define scheduler_maxsubdepth                3

/* Flags . */
#define scheduler_flag_none                  0
#define scheduler_flag_steal                 1
#define scheduler_flag_maxsteal              2


/* Data of a scheduler. */
struct scheduler {

    /* Scheduler flags. */
    unsigned int flags;

    /* Number of queues in this scheduler. */
    int nr_queues;
    
    /* Array of queues. */
    struct queue *queues;
    
    /* Total number of tasks. */
    int nr_tasks, size, tasks_next;
    
    /* Total number of waiting tasks. */
    int waiting;
    
    /* The task array. */
    struct task *tasks;
    
    /* The task indices. */
    int *tasks_ind;
    
    /* Lock for this scheduler. */
    lock_type lock;
    
    /* Waiting queue. */
    pthread_mutex_t sleep_mutex;
    pthread_cond_t sleep_cond;
    
    /* The space associated with this scheduler. */
    struct space *space;

    };


/* Function prototypes. */
void scheduler_init ( struct scheduler *s , struct space *space , int nr_queues , unsigned int flags );
struct task *scheduler_gettask ( struct scheduler *s , int qid );
void scheduler_enqueue ( struct scheduler *s , struct task *t );
void scheduler_start ( struct scheduler *s );
void scheduler_reset ( struct scheduler *s , int nr_tasks );
void scheduler_ranktasks ( struct scheduler *s );
struct task *scheduler_addtask ( struct scheduler *s , int type , int subtype , int flags , int wait , struct cell *ci , struct cell *cj , int tight );
void scheduler_splittasks ( struct scheduler *s );
void scheduler_map_mkghosts ( struct cell *c , void *data );
void scheduler_done ( struct scheduler *s , struct task *t );
