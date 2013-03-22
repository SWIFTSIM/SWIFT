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
#define engine_policy_none          0
#define engine_policy_rand          1
#define engine_policy_steal         2
#define engine_policy_keep          4
#define engine_policy_block         8

#define engine_queue_scale          1.2


/* Data structure for the engine. */
struct engine {

    /* Number of threads on which to run. */
    int nr_threads;
    
    /* The space with which the runner is associated. */
    struct space *s;
    
    /* The runner's threads. */
    struct runner *runners;
    
    /* The running policy. */
    int policy;
    
    /* The number of queues. */
    int nr_queues;
    
    /* The queues. */
    struct queue *queues;
    
    /* The maximum dt to step (current). */
    float dt_step;
    
    /* The minimum dt over all particles in the system. */
    float dt_min;
    
    /* The system time step. */
    float dt, dt_orig;
    
    /* The current step number. */
    int step;
    
    /* The current system time. */
    float time;
    
    /* Data for the threads' barrier. */
    pthread_mutex_t barrier_mutex;
    pthread_cond_t barrier_cond;
    int barrier_count;
    
    };


/* Function prototypes. */
void engine_barrier( struct engine *e );
void engine_init ( struct engine *e , struct space *s , float dt , int nr_threads , int nr_queues , int policy );
void engine_prepare ( struct engine *e );
void engine_step ( struct engine *e , int sort_queues );
