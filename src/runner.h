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



/* The timers themselves. */
enum {
    runner_timer_none = 0,
    runner_timer_dosort,
    runner_timer_doself_density,
    runner_timer_doself_force,
    runner_timer_dopair_density,
    runner_timer_dopair_force,
    runner_timer_dosub_density,
    runner_timer_dosub_force,
    runner_timer_dopair_subset,
    runner_timer_doghost,
    runner_timer_getpair,
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
    #define TIMER_TIC2
    #define TIMER_TOC2(t)
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
#define runner_hist_b 100.0
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


/* A struct representing a runner's thread and its data. */
struct runner {

    /* The id of this thread. */
    int id;

    /* The thread which it is running. */
    pthread_t thread;
    
    /* The underlying runner. */
    struct engine *e;
    
    };


/* Function prototypes. */
void runner_doghost ( struct runner *r , struct cell *c );
void runner_dopair_density ( struct runner *r , struct cell *ci , struct cell *cj );
void runner_doself_density ( struct runner *r , struct cell *c );
void runner_dosub_density ( struct runner *r , struct cell *ci , struct cell *cj , int flags );
void runner_dosort ( struct runner *r , struct cell *c , int flag );
void *runner_main ( void *data );
