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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <float.h>
#include <limits.h>
#include <math.h>

/* Local headers. */
#include "cycle.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "cell.h"


/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/* The timers. */
ticks cell_timer[ cell_timer_count ];


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
        __sync_add_and_fetch( &cell_timer[t] , d );
        return d;
        }
#else
    #define TIMER_TIC
    #define TIMER_TOC(t)
#endif

/**
 * @brief Lock a cell and hold its parents.
 *
 * @param c The #cell.
 */
 
int cell_locktree( struct cell *c ) {

    struct cell *finger, *finger2;
    TIMER_TIC

    /* First of all, try to lock this cell. */
    if ( lock_trylock( &c->lock ) != 0 ) {
        TIMER_TOC(cell_timer_tree);
        return 1;
        }
        
    /* Did somebody hold this cell in the meantime? */
    if ( c->hold ) {
        
        /* Unlock this cell. */
        if ( lock_unlock( &c->lock ) != 0 )
            error( "Failed to unlock cell." );
            
        /* Admit defeat. */
        TIMER_TOC(cell_timer_tree);
        return 1;
    
        }
        
    /* Climb up the tree and lock/hold/unlock. */
    for ( finger = c->parent ; finger != NULL ; finger = finger->parent ) {
    
        /* Lock this cell. */
        if ( lock_trylock( &finger->lock ) != 0 )
            break;
            
        /* Increment the hold. */
        __sync_fetch_and_add( &finger->hold , 1 );
        
        /* Unlock the cell. */
        if ( lock_unlock( &finger->lock ) != 0 )
            error( "Failed to unlock cell." );
    
        }
        
    /* If we reached the top of the tree, we're done. */
    if ( finger == NULL ) {
        TIMER_TOC(cell_timer_tree);
        return 0;
        }
        
    /* Otherwise, we hit a snag. */
    else {
    
        /* Undo the holds up to finger. */
        for ( finger2 = c->parent ; finger2 != finger ; finger2 = finger2->parent )
            __sync_fetch_and_sub( &finger2->hold , 1 );
            
        /* Unlock this cell. */
        if ( lock_unlock( &c->lock ) != 0 )
            error( "Failed to unlock cell." );
            
        /* Admit defeat. */
        TIMER_TOC(cell_timer_tree);
        return 1;
    
        }

    }
    
    
/**
 * @brief Unock a cell's parents.
 *
 * @param c The #cell.
 */
 
void cell_unlocktree( struct cell *c ) {

    struct cell *finger;
    TIMER_TIC

    /* First of all, try to unlock this cell. */
    if ( lock_unlock( &c->lock ) != 0 )
        error( "Failed to unlock cell." );
        
    /* Climb up the tree and unhold the parents. */
    for ( finger = c->parent ; finger != NULL ; finger = finger->parent )
        __sync_fetch_and_sub( &finger->hold , 1 );
        
    TIMER_TOC(cell_timer_tree);
        
    }
    
    
/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 */
 
void cell_split ( struct cell *c  ) {

    int i, j, k, kk;
    struct part temp, *parts = c->parts;
    int left[8], right[8];
    double pivot[3];
    
    /* Init the pivot. */
    for ( k = 0 ; k < 3 ; k++ )
        pivot[k] = c->loc[k] + c->h[k]/2;
    
    /* Split along the x-axis. */
    i = 0; j = c->count - 1;
    while ( i <= j ) {
        while ( i <= c->count-1 && parts[i].x[0] <= pivot[0] )
            i += 1;
        while ( j >= 0 && parts[j].x[0] > pivot[0] )
            j -= 1;
        if ( i < j ) {
            temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
            }
        }
    for ( k = 0 ; k <= j ; k++ )
        if ( parts[k].x[0] > pivot[0] )
            error( "cell_split: sorting failed." );
    for ( k = i ; k < c->count ; k++ )
        if ( parts[k].x[0] < pivot[0] )
            error( "cell_split: sorting failed." );
    left[1] = i; right[1] = c->count - 1;
    left[0] = 0; right[0] = j;
    
    /* Split along the y axis, twice. */
    for ( k = 1 ; k >= 0 ; k-- ) {
        i = left[k]; j = right[k];
        while ( i <= j ) {
            while ( i <= right[k] && parts[i].x[1] <= pivot[1] )
                i += 1;
            while ( j >= left[k] && parts[j].x[1] > pivot[1] )
                j -= 1;
            if ( i < j ) {
                temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                }
            }
        for ( kk = left[k] ; kk <= j ; kk++ )
            if ( parts[kk].x[1] > pivot[1] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (left)." );
                }
        for ( kk = i ; kk <= right[k] ; kk++ )
            if ( parts[kk].x[1] < pivot[1] )
                error( "sorting failed (right)." );
        left[2*k+1] = i; right[2*k+1] = right[k];
        left[2*k] = left[k]; right[2*k] = j;
        }

    /* Split along the z axis, four times. */
    for ( k = 3 ; k >= 0 ; k-- ) {
        i = left[k]; j = right[k];
        while ( i <= j ) {
            while ( i <= right[k] && parts[i].x[2] <= pivot[2] )
                i += 1;
            while ( j >= left[k] && parts[j].x[2] > pivot[2] )
                j -= 1;
            if ( i < j ) {
                temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                }
            }
        for ( kk = left[k] ; kk <= j ; kk++ )
            if ( parts[kk].x[2] > pivot[2] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (left)." );
                }
        for ( kk = i ; kk <= right[k] ; kk++ )
            if ( parts[kk].x[2] < pivot[2] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (right)." );
                }
        left[2*k+1] = i; right[2*k+1] = right[k];
        left[2*k] = left[k]; right[2*k] = j;
        }
        
    /* Store the counts and offsets. */
    for ( k = 0 ; k < 8 ; k++ ) {
        c->progeny[k]->count = right[k] - left[k] + 1;
        c->progeny[k]->parts = &c->parts[ left[k] ];
        }
        
    /* Verify a few sub-cells. */
    /* for ( k = 0 ; k < c->progeny[0]->count ; k++ )
        if ( c->progeny[0]->parts[k].x[0] > pivot[0] ||
             c->progeny[0]->parts[k].x[1] > pivot[1] ||
             c->progeny[0]->parts[k].x[2] > pivot[2] )
            error( "Sorting failed (progeny=0)." );
    for ( k = 0 ; k < c->progeny[1]->count ; k++ )
        if ( c->progeny[1]->parts[k].x[0] > pivot[0] ||
             c->progeny[1]->parts[k].x[1] > pivot[1] ||
             c->progeny[1]->parts[k].x[2] <= pivot[2] )
            error( "Sorting failed (progeny=1)." );
    for ( k = 0 ; k < c->progeny[2]->count ; k++ )
        if ( c->progeny[2]->parts[k].x[0] > pivot[0] ||
             c->progeny[2]->parts[k].x[1] <= pivot[1] ||
             c->progeny[2]->parts[k].x[2] > pivot[2] )
            error( "Sorting failed (progeny=2)." ); */

    }


