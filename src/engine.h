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
#define engine_policy_fixdt         16
#define engine_policy_multistep     32
#define engine_policy_cputight      64
#define engine_policy_mpi           128
#define engine_policy_setaffinity   256

#define engine_queue_scale          1.2
#define engine_maxtaskspercell      32
#define engine_maxproxies           64
#define engine_tasksreweight        10


/* The rank of the engine as a global variable (for messages). */
extern int engine_rank;


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
    
    /* The task scheduler. */
    struct scheduler sched;
    
    /* The maximum dt to step (current). */
    float dt_step;
    
    /* The minimum dt over all particles in the system. */
    float dt_min, dt_max;
    
    /* The system time step. */
    float dt, dt_orig;
    
    /* The system energies from the previous step. */
    double ekin, epot;
    
    /* The current step number. */
    int step, nullstep;
    
    /* The number of particles updated in the previous step. */
    int count_step;
    
    /* The current system time. */
    float time;
    
    /* Data for the threads' barrier. */
    pthread_mutex_t barrier_mutex;
    pthread_cond_t barrier_cond;
    volatile int barrier_running, barrier_launch, barrier_launchcount;
    
    /* ID of the node this engine lives on. */
    int nr_nodes, nodeID;
    
    /* Proxies for the other nodes in this simulation. */
    struct proxy *proxies;
    int nr_proxies, *proxy_ind;
    
    /* Force the engine to rebuild? */
    int forcerebuild, forcerepart;
    
    /* How many steps have we done with the same set of tasks? */
    int tasks_age;
    
    };


/* Function prototypes. */
void engine_barrier( struct engine *e , int tid );
void engine_init ( struct engine *e , struct space *s , float dt , int nr_threads , int nr_queues , int nr_nodes , int nodeID , int policy );
void engine_prepare ( struct engine *e );
void engine_step ( struct engine *e );
void engine_maketasks ( struct engine *e );
void engine_split ( struct engine *e , int *grid );
int engine_exchange_strays ( struct engine *e , struct part *parts , struct xpart *xparts , int *ind , int N );
void engine_rebuild ( struct engine *e );
void engine_repartition ( struct engine *e );
void engine_makeproxies ( struct engine *e );
