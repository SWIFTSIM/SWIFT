
/* Get the inlining right. */
#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif
    
#ifdef PTHREAD_LOCK
    #define lock_type pthread_spinlock_t
    #define lock_init( l ) ( pthread_spin_init( l , PTHREAD_PROCESS_PRIVATE ) != 0 )
    #define lock_destroy( l ) ( pthread_spin_destroy( l ) != 0 )
    #define lock_lock( l ) ( pthread_spin_lock( l ) != 0 )
    #define lock_trylock( l ) ( pthread_spin_lock( l ) != 0 )
    #define lock_unlock( l ) ( pthread_spin_unlock( l ) != 0 )
#else
    #define lock_type volatile int
    #define lock_init( l ) ( *l = 0 )
    #define lock_destroy( l ) 0
    INLINE int lock_lock ( volatile int *l ) {
        while ( __sync_val_compare_and_swap( l , 0 , 1 ) != 0 )
            while( *l );
        return 0;
        }
    #define lock_trylock( l ) ( ( *(l) ) ? 1 : __sync_val_compare_and_swap( l , 0 , 1 ) )
    #define lock_unlock( l ) ( *l = 0 )
#endif
