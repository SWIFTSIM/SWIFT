
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
 
__attribute__ ((always_inline)) INLINE void iact ( float r2 , float hi , float hj , float *io , float *jo , int *ic , int *jc ) {

    #define  KERNEL_COEFF_1  2.546479089470f
    #define  KERNEL_COEFF_2  15.278874536822f
    #define  KERNEL_COEFF_3  45.836623610466f
    #define  KERNEL_COEFF_4  30.557749073644f
    #define  KERNEL_COEFF_5  5.092958178941f
    #define  KERNEL_COEFF_6  (-15.278874536822f)
    #define  NORM_COEFF      4.188790204786f

    float r = sqrtf( r2 );
    float ui, uj, wi, wj;
    
    if ( r2 < hi*hi ) {
        
        ui = r / hi;
        if ( ui < 0.5 )
            wi = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (ui - 1.0f) * ui * ui;
        else
            wi = KERNEL_COEFF_5 * (1.0f - ui) * (1.0f - ui) * (1.0 - ui);
        if ( io != NULL )
            *io += NORM_COEFF * wi;
        if ( ic != NULL )
            *ic += 1;
        
        }

    if ( r2 < hj*hj ) {
        
        uj = r / hj;
        if ( uj < 0.5 )
            wj = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (uj - 1.0f) * uj * uj;
        else
            wj = KERNEL_COEFF_5 * (1.0f - uj) * (1.0f - uj) * (1.0 - uj);
        if ( jo != NULL )
            *jo += NORM_COEFF * wj;
        if ( jc != NULL )
            *jc += 1;
            
        }
    
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
