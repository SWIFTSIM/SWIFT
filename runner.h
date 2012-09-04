/*******************************************************************************
 * This file is part of GadgetSMP.
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
#define runner_policy_none          0
#define runner_policy_rand          1
#define runner_policy_steal         2
#define runner_policy_keep          4
#define runner_policy_block         8

#define runner_queue_scale          1.2


/* The timers themselves. */
enum {
    runner_timer_none = 0,
    runner_timer_dosort,
    runner_timer_doself,
    runner_timer_dopair,
    runner_timer_dosub,
    runner_timer_getpair,
    runner_timer_queue,
    runner_timer_tree,
    runner_timer_bubble,
    runner_timer_steal,
    runner_timer_stalled,
    runner_timer_count,
    };
extern ticks runner_timer[ runner_timer_count ];


/* Define the timer macros. */
#ifdef TIMER_VERBOSE
    #define TIMER
#endif
#ifdef TIMER
    #define TIMER_TIC ticks tic = getticks();
    #define TIMER_TOC(t) timer_toc( t , tic )
    #define TIMER_TIC2 ticks tic2 = getticks();
    #define TIMER_TOC2(t) timer_toc( t , tic2 )
    #ifndef INLINE
    # if __GNUC__ && !__GNUC_STDC_INLINE__
    #  define INLINE extern inline
    # else
    #  define INLINE inline
    # endif
    #endif
    INLINE ticks timer_toc ( int t , ticks tic ) {
        ticks d = (getticks() - tic);
        __sync_add_and_fetch( &runner_timer[t] , d );
        return d;
        }
#else
    #define TIMER_TIC
    #define TIMER_TOC(t)
#endif


/* Counters. */
enum {
    runner_counter_swap = 0,
    runner_counter_stall,
    runner_counter_steal_stall,
    runner_counter_steal_empty,
    runner_counter_keep,
    runner_counter_count,
    };
extern int runner_counter[ runner_counter_count ];


/* Counter macros. */
#ifdef COUNTER
    #define COUNT(c) ( __sync_add_and_fetch( &runner_counter[ c ] , 1 ) )
#else
    #define COUNT(c)
#endif


/* Histogram functions. */
#define runner_hist_a 1.0
#define runner_hist_b 10.0
#define runner_hist_N 99
long long int runner_hist_bins[ runner_hist_N ];
#define runner_hist_hit( x ) __sync_add_and_fetch( &runner_hist_bins[ (int)fmax( 0.0 , fmin( runner_hist_N-1 , ((x) - runner_hist_a) / (runner_hist_b - runner_hist_a) * runner_hist_N ) ) ] , 1 )


/* Get the inlining right. */
#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif


/**
 * @brief Compute the 'interaction' between two particles.
 *
 * @param r2 the inter-particle radius squared.
 * @param hi the screening distance of the ith particle.
 * @param hj the screening distance of the jth particle.
 * @param io Pointer to where to store the interaction of the ith particle.
 * @param jo Pointer to where to store the interaction of the ith particle.
 */
 
__attribute__ ((always_inline)) INLINE void iact_nopart ( float r2 , float hi , float hj , float *force_i , float *force_j , int *count_i , int *count_j ) {

    #define  KERNEL_COEFF_1  2.546479089470f
    #define  KERNEL_COEFF_2  15.278874536822f
    #define  KERNEL_COEFF_3  45.836623610466f
    #define  KERNEL_COEFF_4  30.557749073644f
    #define  KERNEL_COEFF_5  5.092958178941f
    #define  KERNEL_COEFF_6  (-15.278874536822f)
    #define  NORM_COEFF      4.188790204786f

    float r = sqrtf( r2 );
    float ui, uj, wi, wj;
    
    if ( r2 < hi*hi && !( force_i == NULL && count_i == NULL ) ) {
        
        ui = r / hi;
        if ( ui < 0.5 )
            wi = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (ui - 1.0f) * ui * ui;
        else
            wi = KERNEL_COEFF_5 * (1.0f - ui) * (1.0f - ui) * (1.0 - ui);
        if ( force_i != NULL )
            *force_i += NORM_COEFF * wi;
        if ( count_i != NULL )
            *count_i += 1;
        
        }

    if ( r2 < hj*hj && !( force_j == NULL && count_j == NULL ) ) {
        
        uj = r / hj;
        if ( uj < 0.5 )
            wj = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (uj - 1.0f) * uj * uj;
        else
            wj = KERNEL_COEFF_5 * (1.0f - uj) * (1.0f - uj) * (1.0 - uj);
        if ( force_j != NULL )
            *force_j += NORM_COEFF * wj;
        if ( count_j != NULL )
            *count_j += 1;
            
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    


__attribute__ ((always_inline)) INLINE void iact ( float r2 , float hi , float hj , struct part *pi , struct part *pj ) {

    #define  KERNEL_COEFF_1  2.546479089470f
    #define  KERNEL_COEFF_2  15.278874536822f
    #define  KERNEL_COEFF_3  45.836623610466f
    #define  KERNEL_COEFF_4  30.557749073644f
    #define  KERNEL_COEFF_5  5.092958178941f
    #define  KERNEL_COEFF_6  (-15.278874536822f)
    #define  NORM_COEFF      4.188790204786f

    float r = sqrtf( r2 );
    float ui, uj, wi, wj;
    float ui_dh, uj_dh, wi_dh, wj_dh;
    
    if ( r2 < hi*hi && pi != NULL ) {
        
        ui = r / hi;
        ui_dh = -r / hi / hi;
        if ( ui < 0.5f ) {
            wi = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (ui - 1.0f) * ui * ui;
            wi_dh = KERNEL_COEFF_2 * ui_dh * ui * ui
                  + 2 * KERNEL_COEFF_2 * (ui - 1.0f) * ui_dh * ui;
            }
        else {
            wi = KERNEL_COEFF_5 * (1.0f - ui) * (1.0f - ui) * (1.0f - ui);
            wi_dh = -3 * KERNEL_COEFF_5 * ui_dh * (1.0f - ui) * (1.0f - ui);
            }
        pi->count += NORM_COEFF * wi;
        pi->count_dh += NORM_COEFF * wi_dh;
        pi->icount += 1;
        
        }

    if ( r2 < hj*hj && pj != NULL ) {
        
        uj = r / hj;
        uj_dh = -r / hj / hj;
        if ( uj < 0.5f ) {
            wj = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (uj - 1.0f) * uj * uj;
            wj_dh = KERNEL_COEFF_2 * uj_dh * uj * uj
                  + 2 * KERNEL_COEFF_2 * (uj - 1.0f) * uj_dh * uj;
            }
        else {
            wj = KERNEL_COEFF_5 * (1.0f - uj) * (1.0f - uj) * (1.0f - uj);
            wj_dh = -3 * KERNEL_COEFF_5 * uj_dh * (1.0f - uj) * (1.0f - uj);
            }
        pj->count += NORM_COEFF * wj;
        pj->count_dh += NORM_COEFF * wj_dh;
        pj->icount += 1;
            
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    


/* A task queue. */
struct queue {

    /* The lock to access this queue. */
    lock_type lock;

    /* Size, count and next element. */
    int size, count, next;
    
    /* The runner in which this queue lives. */
    struct runner *r;
    
    /* The actual tasks to which the indices refer. */
    struct task *tasks;
    
    /* The task indices. */
    int *tid;

    } __attribute__((aligned (64)));
    

/* A struct representing a runner's thread and its data. */
struct runner_thread {

    /* The id of this thread. */
    int id;

    /* The thread which it is running. */
    pthread_t thread;
    
    /* The underlying runner. */
    struct runner *r;
    
    };


/* Data structure for the runner. */
struct runner {

    /* Number of threads on which to run. */
    int nr_threads;
    
    /* The space with which the runner is associated. */
    struct space *s;
    
    /* The runner's threads. */
    struct runner_thread *threads;
    
    /* The running policy. */
    int policy;
    
    /* The number of queues. */
    int nr_queues;
    
    /* The queues. */
    struct queue *queues;
    
    /* Data for the threads' barrier. */
    pthread_mutex_t barrier_mutex;
    pthread_cond_t barrier_cond;
    int barrier_count;
    
    };


/* Function prototypes. */
void runner_run ( struct runner *r , int sort_queues );
void runner_dopair ( struct runner_thread *rt , struct cell *ci , struct cell *cj );
void runner_doself ( struct runner_thread *rt , struct cell *c );
void runner_dosort ( struct runner_thread *rt , struct cell *c , int flag );
void runner_init ( struct runner *r , struct space *s , int nr_threads , int nr_queues , int policy );
struct task *queue_gettask ( struct queue *q , int blocking , int keep );
