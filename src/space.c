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
#include "space.h"
#include "runner.h"

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/* Split size. */
int space_splitsize = space_splitsize_default;

/* Map shift vector to sortlist. */
const int sortlistID[27] = {
    /* ( -1 , -1 , -1 ) */   0 ,
    /* ( -1 , -1 ,  0 ) */   1 , 
    /* ( -1 , -1 ,  1 ) */   2 ,
    /* ( -1 ,  0 , -1 ) */   3 ,
    /* ( -1 ,  0 ,  0 ) */   4 , 
    /* ( -1 ,  0 ,  1 ) */   5 ,
    /* ( -1 ,  1 , -1 ) */   6 ,
    /* ( -1 ,  1 ,  0 ) */   7 , 
    /* ( -1 ,  1 ,  1 ) */   8 ,
    /* (  0 , -1 , -1 ) */   9 ,
    /* (  0 , -1 ,  0 ) */   10 , 
    /* (  0 , -1 ,  1 ) */   11 ,
    /* (  0 ,  0 , -1 ) */   12 ,
    /* (  0 ,  0 ,  0 ) */   0 , 
    /* (  0 ,  0 ,  1 ) */   12 ,
    /* (  0 ,  1 , -1 ) */   11 ,
    /* (  0 ,  1 ,  0 ) */   10 , 
    /* (  0 ,  1 ,  1 ) */   9 ,
    /* (  1 , -1 , -1 ) */   8 ,
    /* (  1 , -1 ,  0 ) */   7 , 
    /* (  1 , -1 ,  1 ) */   6 ,
    /* (  1 ,  0 , -1 ) */   5 ,
    /* (  1 ,  0 ,  0 ) */   4 , 
    /* (  1 ,  0 ,  1 ) */   3 ,
    /* (  1 ,  1 , -1 ) */   2 ,
    /* (  1 ,  1 ,  0 ) */   1 , 
    /* (  1 ,  1 ,  1 ) */   0 
    };
    
    
/**
 * @breif Recursively dismantle a cell tree.
 *
 */
 
void space_rebuild_recycle ( struct space *s , struct cell *c ) {
    
    int k;
    
    if ( c->split )
        for ( k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL ) {
                space_rebuild_recycle( s , c->progeny[k] );
                space_recycle( s , c->progeny[k] );
                c->progeny[k] = NULL;
                }
    
    }

/**
 * @breif Recursively rebuild a cell tree.
 *
 */
 
int space_rebuild_recurse ( struct space *s , struct cell *c ) {
    
    int k, count, changes = 0, wasmt[8];
    float h, h_limit, h_max = 0.0f;
    struct cell *temp;
    
    /* If the cell is already split, check that the split is still ok. */
    if ( c->split ) {
    
        /* Check the depth. */
        if ( c->depth > s->maxdepth )
            s->maxdepth = c->depth;

        /* Set the minimum cutoff. */
        h_limit = fmin( c->h[0] , fmin( c->h[1] , c->h[2] ) ) / 2;

        /* Count the particles below that. */
        for ( count = 0 , k = 0 ; k < c->count ; k++ ) {
            h = c->parts[k].h;
            if ( h <= h_limit )
                count += 1;
            if ( h > h_max )
                h_max = h;
            }
        c->h_max = h_max;
            
        /* Un-split? */
        if ( count < c->count*space_splitratio || c->count < space_splitsize ) {
        
            /* Get rid of the progeny. */
            space_rebuild_recycle( s , c );
            
            /* Re-set the split flag. */
            c->split = 0;
        
            }
        
        /* Otherwise, recurse on the kids. */
        else {
        
            /* Populate all progeny. */
            for ( k = 0 ; k < 8 ; k++ )
                if ( ( wasmt[k] = ( c->progeny[k] == NULL ) ) ) {
                    temp = space_getcell( s );
                    temp->count = 0;
                    temp->loc[0] = c->loc[0];
                    temp->loc[1] = c->loc[1];
                    temp->loc[2] = c->loc[2];
                    temp->h[0] = c->h[0]/2;
                    temp->h[1] = c->h[1]/2;
                    temp->h[2] = c->h[2]/2;
                    if ( k & 4 )
                        temp->loc[0] += temp->h[0];
                    if ( k & 2 )
                        temp->loc[1] += temp->h[1];
                    if ( k & 1 )
                        temp->loc[2] += temp->h[2];
                    temp->depth = c->depth + 1;
                    temp->split = 0;
                    temp->h_max = 0.0;
                    temp->parent = c;
                    c->progeny[k] = temp;
                    }
        
            /* Make sure each part is in its place. */
            cell_split( c );
            
            /* Remove empty progeny. */
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k]->count == 0 ) {
                    changes += !wasmt[k];
                    space_recycle( s , c->progeny[k] );
                    c->progeny[k] = NULL;
                    }
                else
                    changes += wasmt[k];
        
            /* Recurse. */
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    changes += space_rebuild_recurse( s , c->progeny[k] );
                    
            }
    
        }
        
    /* Otherwise, try to split it anyway. */
    else {
        space_split( s , c );
        changes += c->split;
        }
        
    /* Return the grand total. */
    return changes;
    
    }

/**
 * @breif Re-build the cells as well as the tasks.
 *
 * @param s The #space in which to update the cells.
 * @param force Flag to force re-building the cells and tasks.
 *
 * @return 1 if changes to the cells and/or tasks were made.
 */
 
int space_rebuild ( struct space *s , int force ) {

    float h_max = 0.0f;
    int i, j, k, cdim[3];
    struct cell *c;
    int changes = 0;
    
    /* Run through the parts and get the current h_max. */
    for ( k = 0 ; k < s->nr_parts ; k++ )
        if ( s->parts[k].h > h_max )
            h_max = s->parts[k].h;
    
    /* Get the new putative cell dimensions. */
    for ( k = 0 ; k < 3 ; k++ )
        cdim[k] = floor( s->dim[k] / h_max );
        
    /* Do we need to re-build the upper-level cells? */
    if ( force || cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] || cdim[2] < s->cdim[2] ) {
    
        /* Free the old cells, if they were allocated. */
        if ( s->cells != NULL ) {
            for ( k = 0 ; k < s->nr_cells ; k++ )
                space_rebuild_recycle( s , &s->cells[k] );
            free( s->cells );
            s->maxdepth = 0;
            }
            
        /* Set the new cell dimensions. */
        for ( k = 0 ; k < 3 ; k++ ) {
            s->cdim[k] = cdim[k];
            s->h[k] = s->dim[k] / cdim[k];
            s->ih[k] = 1.0 / s->h[k];
            }
    
        /* Allocate the highest level of cells. */
        s->nr_cells = cdim[0] * cdim[1] * cdim[2];
        if ( posix_memalign( (void *)&s->cells , 64 , s->nr_cells * sizeof(struct cell) ) != 0 )
            error( "Failed to allocate cells." );
        bzero( s->cells , s->nr_cells * sizeof(struct cell) );
        for ( k = 0 ; k < s->nr_cells ; k++ )
            if ( lock_init( &s->cells[k].lock ) != 0 )
                error( "Failed to init spinlock." );

        /* Set the cell location and sizes. */
        for ( i = 0 ; i < cdim[0] ; i++ )
            for ( j = 0 ; j < cdim[1] ; j++ )
                for ( k = 0 ; k < cdim[2] ; k++ ) {
                    c = &s->cells[ cell_getid( cdim , i , j , k ) ];
                    c->loc[0] = i*s->h[0]; c->loc[1] = j*s->h[1]; c->loc[2] = k*s->h[2];
                    c->h[0] = s->h[0]; c->h[1] = s->h[1]; c->h[2] = s->h[2];
                    c->depth = 0;
                    }
                    
        /* There were massive changes. */
        changes = 1;
        
        } /* re-build upper-level cells? */
        
    /* At this point, we have the upper-level cells, old or new. Now make
       sure that the parts in each cell are ok. */
    for ( k = 0 ; k < s->nr_cells ; k++ )
        changes += space_rebuild_recurse( s , &s->cells[k] );
        
    /* Now that we have the cell structre, re-build the tasks. */
    if ( changes )
        space_maketasks( s , 1 );
    
    /* Return the number of changes. */
    return changes;

    }


/**
 * @brief Sort the particles according to the given indices.
 *
 * @param parts The list of #part
 * @param ind The indices with respect to which the parts are sorted.
 * @param N The number of parts
 * @param min Lowest index.
 * @param max highest index.
 *
 * This function calls itself recursively.
 */
 
void parts_sort ( struct part *parts , int *ind , int N , int min , int max ) {

    int pivot = (min + max) / 2;
    int i = 0, j = N-1;
    int temp_i;
    struct part temp_p;
    
    /* One pass of quicksort. */
    while ( i < j ) {
        while ( i < N && ind[i] <= pivot )
            i++;
        while ( j >= 0 && ind[j] > pivot )
            j--;
        if ( i < j ) {
            temp_i = ind[i]; ind[i] = ind[j]; ind[j] = temp_i;
            temp_p = parts[i]; parts[i] = parts[j]; parts[j] = temp_p;
            }
        }
        
    /* Verify sort. */
    for ( int k = 0 ; k <= j ; k++ )
        if ( ind[k] > pivot ) {
            printf( "parts_sort: sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i, N=%i.\n" , k , ind[k] , pivot , i , j , N );
            error( "Sorting failed (<=pivot)." );
            }
    for ( int k = j+1 ; k < N ; k++ )
        if ( ind[k] <= pivot ) {
            printf( "parts_sort: sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%i, j=%i, N=%i.\n" , k , ind[k] , pivot , i , j , N );
            error( "Sorting failed (>pivot)." );
            }
        
    /* Recurse on the left? */
    if ( j > 0 && pivot > min )
        parts_sort( parts , ind , j+1 , min , pivot );
        
    /* Recurse on the right? */
    if ( i < N && pivot+1 < max )
        parts_sort( &parts[i], &ind[i], N-i , pivot+1 , max );

    }


/**
 * @brief Mapping function to free the sorted indices buffers.
 */

void space_map_clearsort ( struct cell *c , void *data ) {

    if ( c->sort != NULL ) {
        free( c->sort );
        c->sort = NULL;
        }

    }


/**
 * @brief Mapping function to append a ghost task to each cell.
 *
 * Looks for the super cell, e.g. the highest-level cell above each
 * cell for which a pair is defined. All ghosts below this cell will
 * depend on the ghost of their parents (sounds spooky, but it isn't).
 */

void space_map_mkghosts ( struct cell *c , void *data ) {

    struct space *s = (struct space *)data;
    struct cell *finger;

    /* Find the super cell, i.e. the highest cell hierarchically above
       this one to still have at least one task associated with it. */
    c->super = c;
    for ( finger = c->parent ; finger != NULL ; finger = finger->parent )
        if ( finger->nr_tasks > 0 )
            c->super = finger;
            
    /* Make the ghost task */
    if ( c->super != c || c->nr_tasks > 0 )
        c->ghost = space_addtask( s , task_type_ghost , task_subtype_none , 0 , 0 , c , NULL , NULL , 0 , NULL , 0 );

    /* If we are not the super cell ourselves, make our ghost depend
       on our parent cell. */
    if ( c->super != c )
        task_addunlock( c->parent->ghost , c->ghost );
    
    }


/**
 * @brief Mapping function to clear the number of tasks in each cell.
 */

void space_map_clearnrtasks ( struct cell *c , void *data ) {

    c->nr_tasks = 0;
    c->nr_density = 0;

    }


/**
 * @brief Get a task free of dependencies and conflicts.
 *
 * @param s The #space.
 */
 
struct task *space_gettask ( struct space *s ) {

    int k, tid = -1;
    struct task *res = NULL;
    struct cell *c;

    /* Main loop, while there are tasks... */
    while ( s->next_task < s->nr_tasks ) {
    
        /* Grab the task lock. */
        if ( lock_lock( &s->task_lock ) != 0 )
            error( "Locking the task_lock failed.\n" );
            
        /* Loop over the remaining task IDs. */
        for ( k = s->next_task ; k < s->nr_tasks ; k++ ) {
        
            /* Put a finger on the task. */
            res = &s->tasks[ s->tasks_ind[k] ];
            
            /* Is this task blocked? */
            if ( res->wait )
                continue;
            
            /* Different criteria for different types. */
            switch ( res->type ) {
                case task_type_self:
                    if ( res->ci->lock || res->ci->wait > 0 )
                        continue;
                    break;
                case task_type_pair:
                    if ( res->ci->lock || res->cj->lock || res->ci->wait || res->cj->wait )
                        continue;
                    break;
                case task_type_sort:
                    if ( res->ci->lock )
                        continue;
                    break;
                }
            
            /* If we made it this far, we're safe. */
            break;
        
            } /* loop over the task IDs. */
            
        /* Did we get a task? */
        if ( k < s->nr_tasks ) {
        
            // /* Swap to front. */
            // tid = s->tasks_ind[k];
            // s->tasks_ind[k] = s->tasks_ind[ s->next_task ];
            // s->tasks_ind[ s->next_task ] = tid;
            
            /* Bubble-down the task. */
            tid = s->tasks_ind[k];
            while ( k > s->next_task ) {
                s->tasks_ind[ k ] = s->tasks_ind[ k-1 ];
                k -= 1;
                }
            s->tasks_ind[ s->next_task ] = tid;
            
            /* Lock the cells, if needed. */
            if ( s->tasks[tid].type != task_type_sort ) {
                for ( c = res->ci ; c != NULL ; c = c->parent )
                    __sync_fetch_and_add( &c->lock , 1 );
                for ( c = res->cj ; c != NULL ; c = c->parent )
                    __sync_fetch_and_add( &c->lock , 1 );
                }
            
            /* up the counter. */
            s->next_task += 1;
        
            }
    
        /* Release the task lock. */
        if ( lock_unlock( &s->task_lock ) != 0 )
            error( "Locking the task_lock failed.\n" );
            
        /* Leave? */
        if ( tid >= 0 )
            return &s->tasks[tid];
    
        } /* while there are tasks. */
        
    /* No beef. */
    return NULL;

    }


/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 */
 
void space_map_parts ( struct space *s , void (*fun)( struct part *p , struct cell *c , void *data ) , void *data ) {

    int i;

    void rec_map ( struct cell *c ) {
    
        int k;
        
        /* No progeny? */
        if ( !c->split )
            for ( k = 0 ; k < c->count ; k++ )
                fun( &c->parts[k] , c , data );
                
        /* Otherwise, recurse. */
        else
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    rec_map( c->progeny[k] );
                
        }
        
    /* Call the recursive function on all higher-level cells. */
    for ( i = 0 ; i < s->nr_cells ; i++ )
        rec_map( &s->cells[i] );

    }


/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 * @param full Map to all cells, including cells with sub-cells.
 */
 
void space_map_cells ( struct space *s , int full , void (*fun)( struct cell *c , void *data ) , void *data ) {

    int i;

    void rec_map ( struct cell *c ) {
    
        int k;
        
        /* No progeny? */
        if ( full || !c->split )
            fun( c , data );
                
        /* Recurse. */
        if ( c->split )
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    rec_map( c->progeny[k] );
                
        }
        
    /* Call the recursive function on all higher-level cells. */
    for ( i = 0 ; i < s->nr_cells ; i++ )
        rec_map( &s->cells[i] );

    }


/**
 * @brief Add a #task to the #space.
 *
 * @param s The #space we are working in.
 */
 
struct task *space_addtask ( struct space *s , int type , int subtype , int flags , int wait , struct cell *ci , struct cell *cj , struct task *unlock_tasks[] , int nr_unlock_tasks , struct cell *unlock_cells[] , int nr_unlock_cells ) {

    struct task *t = &s->tasks[ s->nr_tasks ];
    
    /* Copy the data. */
    t->type = type;
    t->subtype = subtype;
    t->flags = flags;
    t->wait = wait;
    t->ci = ci;
    t->cj = cj;
    if ( unlock_tasks != NULL )
        memcpy( t->unlock_tasks , unlock_tasks , sizeof(struct task *) * nr_unlock_tasks );
    t->nr_unlock_tasks = nr_unlock_tasks;
    if ( unlock_cells != NULL )
        memcpy( t->unlock_cells , unlock_cells , sizeof(struct task *) * nr_unlock_cells );
    t->nr_unlock_cells = nr_unlock_cells;
    
    /* Increase the task counter. */
    s->nr_tasks += 1;
    
    /* Return a pointer to the new task. */
    return t;

    }



/**
 * @brief Split tasks that may be too large.
 *
 * @param s The #space we are working in.
 */
 
void space_splittasks ( struct space *s ) {

    int k, sid, tid;
    struct cell *ci, *cj;
    double hi, hj, shift[3];
    struct task *t;

    /* Loop through the tasks... */
    for ( tid = 0 ; tid < s->nr_tasks ; tid++ ) {
    
        /* Get a pointer on the task. */
        t = &s->tasks[tid];
    
        /* If this task isn't a pair, i'm not interested. */
        if ( t->type != task_type_pair )
            continue;
            
        /* Get a handle on the cells involved. */
        ci = t->ci;
        cj = t->cj;
        hi = fmax( ci->h[0] , fmax( ci->h[1] , ci->h[2] ) );
        hj = fmax( cj->h[0] , fmax( cj->h[1] , cj->h[2] ) );
            
        /* Should this task be split-up? */
        if ( ci->split && cj->split && ci->h_max*space_stretch < hi/2 && cj->h_max*space_stretch < hj/2 ) {
        
            /* Get the relative distance between the pairs, wrapping. */
            for ( k = 0 ; k < 3 ; k++ ) {
                if ( cj->loc[k] - ci->loc[k] < -s->dim[k]/2 )
                    shift[k] = s->dim[k];
                else if ( cj->loc[k] - ci->loc[k] > s->dim[k]/2 )
                    shift[k] = -s->dim[k];
                else
                    shift[k] = 0.0;
                }

            /* Get the sorting index. */
            for ( sid = 0 , k = 0 ; k < 3 ; k++ )
                sid = 3*sid + ( (cj->loc[k] - ci->loc[k] + shift[k] < 0) ? 0 : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1 );

            /* Flip? */
            if ( sid < 13 ) {
                cj = t->ci;
                ci = t->cj;
                t->ci = ci; t->cj = cj;
                }
            else
                sid = 26 - sid;
                
            /* Remove the dependency of this task on the sorts of ci and cj. */
            task_rmunlock( ci->sorts[sid] , t );
            task_rmunlock( cj->sorts[sid] , t );
            ci->nr_pairs -= 1;
            cj->nr_pairs -= 1;
            t->nr_unlock_cells = 0;
                
            /* For each different sorting type... */
            switch ( sid ) {
            
                case 0: /* (  1 ,  1 ,  1 ) */
                    t->ci = ci->progeny[7]; t->cj = cj->progeny[0];
                    task_addunlock( ci->progeny[7]->sorts[0] , t ); task_addunlock( cj->progeny[0]->sorts[0] , t );
                    ci->progeny[7]->nr_pairs += 1;
                    cj->progeny[0]->nr_pairs += 1;
                    break;
                    
                case 1: /* (  1 ,  1 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split ) {
                        t->type = task_type_sub; t->flags = 1;
                        task_addunlock( ci->progeny[6]->sorts[1] , t ); task_addunlock( cj->progeny[0]->sorts[1] , t );
                        task_addunlock( ci->progeny[7]->sorts[1] , t ); task_addunlock( cj->progeny[1]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[0] , t ); task_addunlock( cj->progeny[1]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        }
                    else {
                        t->ci = ci->progeny[6]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[6]->sorts[1] , t ); task_addunlock( cj->progeny[0]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[1] , t ); task_addunlock( cj->progeny[1]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[0] , t ); task_addunlock( cj->progeny[1]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        }
                    ci->progeny[6]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    break;
                    
                case 2: /* (  1 ,  1 , -1 ) */
                    t->ci = ci->progeny[6]; t->cj = cj->progeny[1];
                    task_addunlock( ci->progeny[6]->sorts[2] , t ); task_addunlock( cj->progeny[1]->sorts[2] , t );
                    ci->progeny[6]->nr_pairs += 1;
                    cj->progeny[1]->nr_pairs += 1;
                    break;
                    
                case 3: /* (  1 ,  0 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[5]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[2]->split ) {
                        t->type = task_type_sub; t->flags = 3;
                        task_addunlock( ci->progeny[5]->sorts[3] , t ); task_addunlock( cj->progeny[0]->sorts[3] , t );
                        task_addunlock( ci->progeny[7]->sorts[3] , t ); task_addunlock( cj->progeny[2]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[0] , t ); task_addunlock( cj->progeny[2]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        }
                    else {
                        t->ci = ci->progeny[5]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[5]->sorts[3] , t ); task_addunlock( cj->progeny[0]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[3] , t ); task_addunlock( cj->progeny[2]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[0] , t ); task_addunlock( cj->progeny[2]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        }
                    ci->progeny[5]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[2]->nr_pairs += 2;
                    break;
                    
                case 4: /* (  1 ,  0 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[5]->split && !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split && !cj->progeny[2]->split && !cj->progeny[3]->split ) {
                        t->type = task_type_sub; t->flags = 4;
                        task_addunlock( ci->progeny[4]->sorts[4] , t ); task_addunlock( cj->progeny[0]->sorts[4] , t );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        task_addunlock( ci->progeny[4]->sorts[3] , t ); task_addunlock( cj->progeny[1]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[4] , t ); task_addunlock( cj->progeny[1]->sorts[4] , t );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        task_addunlock( ci->progeny[4]->sorts[1] , t ); task_addunlock( cj->progeny[2]->sorts[1] , t );
                        task_addunlock( ci->progeny[5]->sorts[2] , t ); task_addunlock( cj->progeny[2]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[4] , t ); task_addunlock( cj->progeny[2]->sorts[4] , t );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        task_addunlock( ci->progeny[4]->sorts[0] , t ); task_addunlock( cj->progeny[3]->sorts[0] , t );
                        task_addunlock( ci->progeny[5]->sorts[1] , t ); task_addunlock( cj->progeny[3]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[3] , t ); task_addunlock( cj->progeny[3]->sorts[3] , t );
                        task_addunlock( ci->progeny[7]->sorts[4] , t ); task_addunlock( cj->progeny[3]->sorts[4] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[4]->sorts[4] , t ); task_addunlock( cj->progeny[0]->sorts[4] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[4] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[3] , t ); task_addunlock( cj->progeny[1]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[4] , t ); task_addunlock( cj->progeny[1]->sorts[4] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[1] , t ); task_addunlock( cj->progeny[2]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[2] , t ); task_addunlock( cj->progeny[2]->sorts[2] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[4] , t ); task_addunlock( cj->progeny[2]->sorts[4] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[0] , t ); task_addunlock( cj->progeny[3]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[1] , t ); task_addunlock( cj->progeny[3]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[3] , t ); task_addunlock( cj->progeny[3]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[4] , t ); task_addunlock( cj->progeny[3]->sorts[4] , t );
                        }
                    ci->progeny[4]->nr_pairs += 4;
                    ci->progeny[5]->nr_pairs += 4;
                    ci->progeny[6]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[1]->nr_pairs += 4;
                    cj->progeny[2]->nr_pairs += 4;
                    cj->progeny[3]->nr_pairs += 4;
                    break;
                    
                case 5: /* (  1 ,  0 , -1 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[6]->split &&
                         !cj->progeny[1]->split && !cj->progeny[3]->split ) {
                        t->type = task_type_sub; t->flags = 5;
                        task_addunlock( ci->progeny[4]->sorts[5] , t ); task_addunlock( cj->progeny[1]->sorts[5] , t );
                        task_addunlock( ci->progeny[6]->sorts[5] , t ); task_addunlock( cj->progeny[3]->sorts[5] , t );
                        task_addunlock( ci->progeny[4]->sorts[2] , t ); task_addunlock( cj->progeny[3]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[1];
                        task_addunlock( ci->progeny[4]->sorts[5] , t ); task_addunlock( cj->progeny[1]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[5] , t ); task_addunlock( cj->progeny[3]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[2] , t ); task_addunlock( cj->progeny[3]->sorts[2] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        }
                    ci->progeny[4]->nr_pairs += 2;
                    ci->progeny[6]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    cj->progeny[3]->nr_pairs += 2;
                    break;
                    
                case 6: /* (  1 , -1 ,  1 ) */
                    t->ci = ci->progeny[5]; t->cj = cj->progeny[2];
                    task_addunlock( ci->progeny[5]->sorts[6] , t ); task_addunlock( cj->progeny[2]->sorts[6] , t );
                    ci->progeny[5]->nr_pairs += 1;
                    cj->progeny[2]->nr_pairs += 1;
                    break;
                    
                case 7: /* (  1 , -1 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[5]->split &&
                         !cj->progeny[2]->split && !cj->progeny[3]->split ) {
                        t->type = task_type_sub; t->flags = 7;
                        task_addunlock( ci->progeny[4]->sorts[6] , t ); task_addunlock( cj->progeny[3]->sorts[6] , t );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        task_addunlock( ci->progeny[4]->sorts[7] , t ); task_addunlock( cj->progeny[2]->sorts[7] , t );
                        task_addunlock( ci->progeny[5]->sorts[7] , t ); task_addunlock( cj->progeny[3]->sorts[7] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[3];
                        task_addunlock( ci->progeny[4]->sorts[6] , t ); task_addunlock( cj->progeny[3]->sorts[6] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[7] , t ); task_addunlock( cj->progeny[2]->sorts[7] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[7] , t ); task_addunlock( cj->progeny[3]->sorts[7] , t );
                        }
                    ci->progeny[4]->nr_pairs += 2;
                    ci->progeny[5]->nr_pairs += 2;
                    cj->progeny[2]->nr_pairs += 2;
                    cj->progeny[3]->nr_pairs += 2;
                    break;
                    
                case 8: /* (  1 , -1 , -1 ) */
                    t->ci = ci->progeny[4]; t->cj = cj->progeny[3];
                    task_addunlock( ci->progeny[4]->sorts[8] , t ); task_addunlock( cj->progeny[3]->sorts[8] , t );
                    ci->progeny[4]->nr_pairs += 1;
                    cj->progeny[3]->nr_pairs += 1;
                    break;
                    
                case 9: /* (  0 ,  1 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[3]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[4]->split ) {
                        t->type = task_type_sub; t->flags = 9;
                        task_addunlock( ci->progeny[3]->sorts[9] , t ); task_addunlock( cj->progeny[0]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[9] , t ); task_addunlock( cj->progeny[4]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[0] , t ); task_addunlock( cj->progeny[4]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        }
                    else {
                        t->ci = ci->progeny[3]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[3]->sorts[9] , t ); task_addunlock( cj->progeny[0]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[9] , t ); task_addunlock( cj->progeny[4]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[0] , t ); task_addunlock( cj->progeny[4]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        }
                    ci->progeny[3]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[4]->nr_pairs += 2;
                    break;
                    
                case 10: /* (  0 ,  1 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[2]->split && !ci->progeny[3]->split && !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split && !cj->progeny[4]->split && !cj->progeny[5]->split ) {
                        t->type = task_type_sub; t->flags = 10;
                        task_addunlock( ci->progeny[2]->sorts[10] , t ); task_addunlock( cj->progeny[0]->sorts[10] , t );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        task_addunlock( ci->progeny[2]->sorts[9] , t ); task_addunlock( cj->progeny[1]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[10] , t ); task_addunlock( cj->progeny[1]->sorts[10] , t );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        task_addunlock( ci->progeny[2]->sorts[1] , t ); task_addunlock( cj->progeny[4]->sorts[1] , t );
                        task_addunlock( ci->progeny[3]->sorts[2] , t ); task_addunlock( cj->progeny[4]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[10] , t ); task_addunlock( cj->progeny[4]->sorts[10] , t );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        task_addunlock( ci->progeny[2]->sorts[0] , t ); task_addunlock( cj->progeny[5]->sorts[0] , t );
                        task_addunlock( ci->progeny[3]->sorts[1] , t ); task_addunlock( cj->progeny[5]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[9] , t ); task_addunlock( cj->progeny[5]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[10] , t ); task_addunlock( cj->progeny[5]->sorts[10] , t );
                        }
                    else {
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[2]->sorts[10] , t ); task_addunlock( cj->progeny[0]->sorts[10] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[2] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[9] , t ); task_addunlock( cj->progeny[1]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[10] , t ); task_addunlock( cj->progeny[1]->sorts[10] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[2] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[1] , t ); task_addunlock( cj->progeny[4]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[2] , t ); task_addunlock( cj->progeny[4]->sorts[2] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[10] , t ); task_addunlock( cj->progeny[4]->sorts[10] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[0] , t ); task_addunlock( cj->progeny[5]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[1] , t ); task_addunlock( cj->progeny[5]->sorts[1] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[9] , t ); task_addunlock( cj->progeny[5]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[10] , t ); task_addunlock( cj->progeny[5]->sorts[10] , t );
                        }
                    ci->progeny[2]->nr_pairs += 4;
                    ci->progeny[3]->nr_pairs += 4;
                    ci->progeny[6]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[1]->nr_pairs += 4;
                    cj->progeny[4]->nr_pairs += 4;
                    cj->progeny[5]->nr_pairs += 4;
                    break;
                    
                case 11: /* (  0 ,  1 , -1 ) */
                    if ( space_dosub &&
                         !ci->progeny[2]->split && !ci->progeny[6]->split &&
                         !cj->progeny[1]->split && !cj->progeny[5]->split ) {
                        t->type = task_type_sub; t->flags = 11;
                        task_addunlock( ci->progeny[2]->sorts[11] , t ); task_addunlock( cj->progeny[1]->sorts[11] , t );
                        task_addunlock( ci->progeny[6]->sorts[11] , t ); task_addunlock( cj->progeny[5]->sorts[11] , t );
                        task_addunlock( ci->progeny[2]->sorts[2] , t ); task_addunlock( cj->progeny[5]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        }
                    else {
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[1];
                        task_addunlock( ci->progeny[2]->sorts[11] , t ); task_addunlock( cj->progeny[1]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[11] , t ); task_addunlock( cj->progeny[5]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[2] , t ); task_addunlock( cj->progeny[5]->sorts[2] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        }
                    ci->progeny[2]->nr_pairs += 2;
                    ci->progeny[6]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    cj->progeny[5]->nr_pairs += 2;
                    break;
                    
                case 12: /* (  0 ,  0 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[1]->split && !ci->progeny[3]->split && !ci->progeny[5]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[2]->split && !cj->progeny[4]->split && !cj->progeny[6]->split ) {
                        t->type = task_type_sub; t->flags = 12;
                        task_addunlock( ci->progeny[1]->sorts[12] , t ); task_addunlock( cj->progeny[0]->sorts[12] , t );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        task_addunlock( ci->progeny[1]->sorts[9] , t ); task_addunlock( cj->progeny[2]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[12] , t ); task_addunlock( cj->progeny[2]->sorts[12] , t );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        task_addunlock( ci->progeny[1]->sorts[3] , t ); task_addunlock( cj->progeny[4]->sorts[3] , t );
                        task_addunlock( ci->progeny[3]->sorts[6] , t ); task_addunlock( cj->progeny[4]->sorts[6] , t );
                        task_addunlock( ci->progeny[5]->sorts[12] , t ); task_addunlock( cj->progeny[4]->sorts[12] , t );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        task_addunlock( ci->progeny[1]->sorts[0] , t ); task_addunlock( cj->progeny[6]->sorts[0] , t );
                        task_addunlock( ci->progeny[3]->sorts[3] , t ); task_addunlock( cj->progeny[6]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[9] , t ); task_addunlock( cj->progeny[6]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[12] , t ); task_addunlock( cj->progeny[6]->sorts[12] , t );
                        }
                    else {
                        t->ci = ci->progeny[1]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[1]->sorts[12] , t ); task_addunlock( cj->progeny[0]->sorts[12] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[1] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[9] , t ); task_addunlock( cj->progeny[2]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[12] , t ); task_addunlock( cj->progeny[2]->sorts[12] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[1] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[3] , t ); task_addunlock( cj->progeny[4]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[6] , t ); task_addunlock( cj->progeny[4]->sorts[6] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[12] , t ); task_addunlock( cj->progeny[4]->sorts[12] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[1] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[0] , t ); task_addunlock( cj->progeny[6]->sorts[0] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[3] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[3] , t ); task_addunlock( cj->progeny[6]->sorts[3] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[5] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[9] , t ); task_addunlock( cj->progeny[6]->sorts[9] , t );
                        t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , ci->progeny[7] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[12] , t ); task_addunlock( cj->progeny[6]->sorts[12] , t );
                        }
                    ci->progeny[1]->nr_pairs += 4;
                    ci->progeny[3]->nr_pairs += 4;
                    ci->progeny[5]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[2]->nr_pairs += 4;
                    cj->progeny[4]->nr_pairs += 4;
                    cj->progeny[6]->nr_pairs += 4;
                    break;
            
                }
                
            /* Take a step back... */
            tid -= 1;
            
            } /* split this task? */
    
        } /* loop over all tasks. */
        
    }
    
    
/**
 * @brief Fill the #space's task list.
 *
 * @param s The #space we are working in.
 * @param do_sort Flag to add sorting tasks to the list.
 */
 
void space_maketasks ( struct space *s , int do_sort ) {

    int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd;
    int *cdim = s->cdim;
    struct task *t , *t2;
    int pts[7][8] = { { -1 , 12 , 10 , 9 , 4 , 3 , 1 , 0 } ,
                      { -1 , -1 , 11 , 10 , 5 , 4 , 2 , 1 } ,
                      { -1 , -1 , -1 , 12 , 7 , 6 , 4 , 3 } , 
                      { -1 , -1 , -1 , -1 , 8 , 7 , 5 , 4 } ,
                      { -1 , -1 , -1 , -1 , -1 , 12 , 10 , 9 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , 11 , 10 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , -1 , 12 } };
    int counts[task_type_count];

    /* Recursive function to generate tasks in the cell tree. */
    void maketasks_rec ( struct cell *c , struct task *sort_up[] , int nr_sort_up , struct cell *parent ) {

        int j, k, nr_sort = 0;
        struct task *sort[7], *t;

        /* Clear the waits on this cell. */
        c->wait = 0;
        sort[0] = NULL;
        
        /* Start by generating the sort task. */
        if ( c->count > 0 ) {
        
            if ( do_sort ) {
                if ( c->count < 1000 ) {
                    sort[0] = space_addtask( s , task_type_sort , task_subtype_none , 0x1fff , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    for ( k = 0 ; k < 13 ; k++ )
                        c->sorts[k] = sort[0];
                    nr_sort = 1;
                    }
                else if ( c->count < 5000 ) {
                    sort[0] = space_addtask( s , task_type_sort , task_subtype_none , 0x7f , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    sort[1] = space_addtask( s , task_type_sort , task_subtype_none , 0x1f80 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    for ( k = 0 ; k < 7 ; k++ )
                        c->sorts[k] = sort[0];
                    for ( k = 7 ; k < 14 ; k++ )
                        c->sorts[k] = sort[1];
                    nr_sort = 2;
                    }
                else {
                    c->sorts[0] = c->sorts[1] = sort[0] = space_addtask( s , task_type_sort , task_subtype_none , 0x1 + 0x2 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[2] = c->sorts[3] = sort[1] = space_addtask( s , task_type_sort , task_subtype_none , 0x4 + 0x8 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[4] = c->sorts[5] = sort[2] = space_addtask( s , task_type_sort , task_subtype_none , 0x10 + 0x20 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[6] = c->sorts[7] = sort[3] = space_addtask( s , task_type_sort , task_subtype_none , 0x40 + 0x80 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[8] = c->sorts[9] = sort[4] = space_addtask( s , task_type_sort , task_subtype_none , 0x100 + 0x200 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[10] = c->sorts[11] = sort[5] = space_addtask( s , task_type_sort , task_subtype_none , 0x400 + 0x800 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[12] = c->sorts[13] = sort[6] = space_addtask( s , task_type_sort , task_subtype_none , 0x1000 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    nr_sort = 7;
                    }
                }

            /* Generate a self-interaction if not split. */
            if ( !c->split && c->count > 1 )
                space_addtask( s , task_type_self , task_subtype_density , 0 , 0 , c , NULL , NULL , 0 , NULL , 0 );
                
            }
            
        /* Otherwise, add the interactions between progeny. */
        if ( c->split ) {
        
            /* Recurse. */
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    maketasks_rec( c->progeny[k] , sort , nr_sort , c );
                        
            /* Worth splitting into several tasks? */
            if ( !space_dosub || c->count > 2*space_splitsize ) {
            
                /* Make a task for each pair of progeny. */
                for ( j = 0 ; j < 8 ; j++ )
                    if ( c->progeny[j] != NULL && c->progeny[j]->count > 0 )
                        for ( k = j + 1 ; k < 8 ; k++ )
                            if ( c->progeny[k] != NULL && c->progeny[k]->count > 0 ) {
                                t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , c->progeny[j] , c->progeny[k] , NULL , 0 , NULL , 0 );
                                task_addunlock( c->progeny[j]->sorts[ pts[j][k] ] , t );
                                task_addunlock( c->progeny[k]->sorts[ pts[j][k] ] , t );
                                c->progeny[k]->nr_pairs += 1;
                                c->progeny[j]->nr_pairs += 1;
                                }
                                
                }
                
            /* Otherwise, dispatch as one large task. */
            else {
            
                /* Add the task. */
                t = space_addtask( s , task_type_sub , task_subtype_density , 0 , 0 , c , NULL , NULL , 0 , NULL , 0 );
                
                /* Make it depend on all the sorts of its progeny. */
                for ( k = 0 ; k < 8 ; k++ )
                    for ( j = 0 ; j < 13 ; j++ )
                        if ( c->progeny[k] != NULL )
                            task_addunlock( c->progeny[k]->sorts[j] , t );
            
                }

            }

        }
        
    /* Allocate the task-list, if needed. */
    if ( s->tasks == NULL || s->tasks_size < s->tot_cells * 30 ) {
        if ( s->tasks != NULL )
            free( s->tasks );
        if ( s->tasks_ind != NULL )
            free( s->tasks_ind );
        s->tasks_size = s->tot_cells * 30;
        if ( posix_memalign( (void *)&s->tasks , 64 , sizeof(struct task) * s->tasks_size ) != 0 )
            error( "Failed to allocate task list." );
        if ( ( s->tasks_ind = (int *)malloc( sizeof(int) * s->tasks_size ) ) == NULL )
            error( "Failed to allocate task indices." );
        }
    s->nr_tasks = 0;
    
    /* Loop over the cells and get their sub-tasks. */
    for ( k = 0 ; k < s->nr_cells ; k++ )
        maketasks_rec( &s->cells[k] , NULL , 0 , NULL );

    /* Run through the highest level of cells and add pairs. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                cid = cell_getid( cdim , i , j , k );
                if ( s->cells[cid].count == 0 )
                    continue;
                for ( ii = -1 ; ii < 2 ; ii++ ) {
                    iii = i + ii;
                    if ( !s->periodic && ( iii < 0 || iii >= cdim[0] ) )
                        continue;
                    iii = ( iii + cdim[0] ) % cdim[0];
                    for ( jj = -1 ; jj < 2 ; jj++ ) {
                        jjj = j + jj;
                        if ( !s->periodic && ( jjj < 0 || jjj >= cdim[1] ) )
                            continue;
                        jjj = ( jjj + cdim[1] ) % cdim[1];
                        for ( kk = -1 ; kk < 2 ; kk++ ) {
                            kkk = k + kk;
                            if ( !s->periodic && ( kkk < 0 || kkk >= cdim[2] ) )
                                continue;
                            kkk = ( kkk + cdim[2] ) % cdim[2];
                            cjd = cell_getid( cdim , iii , jjj , kkk );
                            if ( s->cells[cjd].count == 0 )
                                continue;
                            if ( cid >= cjd )
                                continue;
                            t = space_addtask( s , task_type_pair , task_subtype_density , 0 , 0 , &s->cells[cid] , &s->cells[cjd] , NULL , 0 , NULL , 0 );
                            task_addunlock( s->cells[cid].sorts[ sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ] ] , t );
                            task_addunlock( s->cells[cjd].sorts[ sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ] ] , t );
                            s->cells[cid].nr_pairs += 1;
                            s->cells[cjd].nr_pairs += 1;
                            }
                        }
                    }
                }

    /* Split the tasks. */
    space_splittasks( s );
    
    /* Remove sort tasks with no dependencies. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( t->type == task_type_sort && t->nr_unlock_tasks == 0 ) {
            if ( t->ci->split )
                for ( i = 0 ; i < 13 ; i++ )
                    if ( t->flags & ( 1 << i ) ) {
                        for ( j = 0 ; j < 8 ; j++ )
                            if ( t->ci->progeny[j] != NULL )
                                task_rmunlock_blind( t->ci->progeny[j]->sorts[i] , t );
                        t->ci->sorts[i] = NULL;
                        }
            t->type = task_type_none;
            }
        }
        
    /* Count the number of tasks associated with each cell and
       store the density tasks in each cell. */
    space_map_cells( s , 1 , &space_map_clearnrtasks , NULL );
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( t->type == task_type_self ) {
            t->ci->nr_tasks += 1;
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ t->ci->nr_density ] = t;
                t->ci->nr_density += 1;
                }
            }
        else if ( t->type == task_type_pair ) {
            t->ci->nr_tasks += 1;
            t->cj->nr_tasks += 1;
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ t->ci->nr_density ] = t;
                t->ci->nr_density += 1;
                t->cj->density[ t->cj->nr_density ] = t;
                t->cj->nr_density += 1;
                }
            }
        else if ( t->type == task_type_sub ) {
            t->ci->nr_tasks += 1;
            if ( t->cj != NULL )
                t->cj->nr_tasks += 1;
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ t->ci->nr_density ] = t;
                t->ci->nr_density += 1;
                if ( t->cj != NULL ) {
                    t->cj->density[ t->cj->nr_density ] = t;
                    t->cj->nr_density += 1;
                    }
                }
            }
        }
        
    /* Append a ghost task to each cell. */
    space_map_cells( s , 1 , &space_map_mkghosts , s );
    
    /* Run through the tasks and make iacts for each density task. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
    
        /* Get a pointer to the task. */
        t = &s->tasks[k];
        
        /* Self-interaction? */
        if ( t->type == task_type_self && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            t2 = space_addtask( s , task_type_self , task_subtype_force , 0 , 0 , t->ci , NULL , NULL , 0 , NULL , 0 );
            task_addunlock( t->ci->ghost , t2 );
            }
            
        /* Otherwise, pair interaction? */
        else if ( t->type == task_type_pair && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            task_addunlock( t , t->cj->super->ghost );
            t2 = space_addtask( s , task_type_pair , task_subtype_force , 0 , 0 , t->ci , t->cj , NULL , 0 , NULL , 0 );
            task_addunlock( t->ci->ghost , t2 );
            task_addunlock( t->cj->ghost , t2 );
            }

        /* Otherwise, sub interaction? */
        else if ( t->type == task_type_sub && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            if ( t->cj != NULL )
                task_addunlock( t , t->cj->super->ghost );
            t2 = space_addtask( s , task_type_sub , task_subtype_force , t->flags , 0 , t->ci , t->cj , NULL , 0 , NULL , 0 );
            task_addunlock( t->ci->ghost , t2 );
            if ( t->cj != NULL )
                task_addunlock( t->cj->ghost , t2 );
            }
            
        }
        
    /* Re-set the indices. */
    for ( k = 0 ; k < s->nr_tasks ; k++ )
        s->tasks_ind[k] = k;
            
    /* Count the number of each task type. */
    for ( k = 0 ; k < task_type_count ; k++ )
        counts[k] = 0;
    for ( k = 0 ; k < s->nr_tasks ; k++ )
        counts[ s->tasks[k].type ] += 1;
    printf( "space_maketasks: task counts are [ %s=%i" , taskID_names[0] , counts[0] );
    for ( k = 1 ; k < task_type_count ; k++ )
        printf( " %s=%i" , taskID_names[k] , counts[k] );
    printf( " ]\n" );
        
    /* Re-set the next task pointer. */
    s->next_task = 0;
            
    }
    
    

/**
 * @brief Split cells that contain too many particles.
 *
 * @param s The #space we are working in.
 * @param c The #cell under consideration.
 */
 
void space_split ( struct space *s , struct cell *c ) {

    int k, count;
    double h, h_limit, h_max = 0.0;
    struct cell *temp;
    
    /* Check the depth. */
    if ( c->depth > s->maxdepth )
        s->maxdepth = c->depth;
    
    /* Set the minimum cutoff. */
    h_limit = fmin( c->h[0] , fmin( c->h[1] , c->h[2] ) ) / 2;
    
    /* Count the particles below that. */
    for ( count = 0 , k = 0 ; k < c->count ; k++ ) {
        h = c->parts[k].h;
        if ( h <= h_limit )
            count += 1;
        if ( h > h_max )
            h_max = h;
        }
    c->h_max = h_max;
            
    /* Split or let it be? */
    if ( count > c->count*space_splitratio && c->count > space_splitsize ) {
    
        /* No longer just a leaf. */
        c->split = 1;
        
        /* Create the cell's progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            temp = space_getcell( s );
            temp->count = 0;
            temp->loc[0] = c->loc[0];
            temp->loc[1] = c->loc[1];
            temp->loc[2] = c->loc[2];
            temp->h[0] = c->h[0]/2;
            temp->h[1] = c->h[1]/2;
            temp->h[2] = c->h[2]/2;
            if ( k & 4 )
                temp->loc[0] += temp->h[0];
            if ( k & 2 )
                temp->loc[1] += temp->h[1];
            if ( k & 1 )
                temp->loc[2] += temp->h[2];
            temp->depth = c->depth + 1;
            temp->split = 0;
            temp->h_max = 0.0;
            temp->parent = c;
            c->progeny[k] = temp;
            }
            
        /* Split the cell data. */
        cell_split( c );
            
        /* Recurse? */
        for ( k = 0 ; k < 8 ; k++ )
            space_split( s , c->progeny[k] );
            
        /* Remove any progeny with zero parts. */
        for ( k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k]->count == 0 ) {
                space_recycle( s , c->progeny[k] );
                c->progeny[k] = NULL;
                }
                
        }
        
    /* Otherwise, set the progeny to null. */
    else {
        bzero( c->progeny , sizeof(struct cell *) * 8 );
        c->split = 0;
        }

    }


/**
 * @brief Return a used cell to the cell buffer.
 *
 * @param s The #space.
 * @param c The #cell.
 */
 
void space_recycle ( struct space *s , struct cell *c ) {

    /* Clear the cell. */
    if ( lock_destroy( &c->lock ) != 0 )
        error( "Failed to destroy spinlock." );
    
    /* Hook this cell into the buffer. */
    c->next = s->cells_new;
    s->cells_new = c;
    s->tot_cells -= 1;
    
    }


/**
 * @brief Get a new empty cell.
 *
 * @param s The #space.
 */
 
struct cell *space_getcell ( struct space *s ) {

    struct cell *c;
    int k;
    
    /* Is the buffer empty? */
    if ( s->cells_new == NULL ) {
        if ( posix_memalign( (void *)&s->cells_new , 64 , space_cellallocchunk * sizeof(struct cell) ) != 0 )
            error( "Failed to allocate more cells." );
        bzero( s->cells_new , space_cellallocchunk * sizeof(struct cell) );
        for ( k = 0 ; k < space_cellallocchunk-1 ; k++ )
            s->cells_new[k].next = &s->cells_new[k+1];
        s->cells_new[ space_cellallocchunk-1 ].next = NULL;
        }

    /* Pick off the next cell. */
    c = s->cells_new;
    s->cells_new = c->next;
    s->tot_cells += 1;
    
    /* Init some things in the cell. */
    bzero( c , sizeof(struct cell) );
    if ( lock_init( &c->lock ) != 0 )
        error( "Failed to initialize cell spinlock." );
        
    return c;

    }


/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param The #space to initialize.
 * @param dim Spatial dimensions of the domain.
 * @param parts Pointer to an array of #part.
 * @param N The number of parts in the space.
 * @param periodic flag whether the domain is periodic or not.
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_maxppc
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */


void space_init ( struct space *s , double dim[3] , struct part *parts , int N , int periodic , double h_cells ) {

    int i, j, k;
    int nr_cells, cdim[3];
    double h_min, h_max, h[3], ih[3];
    struct cell *c, *cells;
    struct part *finger;
    int *ind;
    
    
    /* Get the minimum and maximum cutoff radii. */
    h_min = parts[0].h; h_max = h_min;
    for ( k = 1 ; k < N ; k++ )
        if ( parts[k].h < h_min )
            h_min = parts[k].h;
        else if ( parts[k].h > h_max )
            h_max = parts[k].h;
            
    /* Stretch the maximum smoothing length. */
    h_max *= space_stretch;
            
    /* Get the cell width. */
    if ( h_cells < h_max )
        h_cells = h_max;
    for ( k = 0 ; k < 3 ; k++ ) {
        cdim[k] = floor( dim[k] / h_cells );
        h[k] = dim[k] / cdim[k];
        ih[k] = 1.0 / h[k];
        }
        
    /* Allocate the highest level of cells. */
    nr_cells = cdim[0] * cdim[1] * cdim[2];
    if ( posix_memalign( (void *)&cells , 64 , nr_cells * sizeof(struct cell) ) != 0 )
        error( "Failed to allocate cells." );
    bzero( cells , nr_cells * sizeof(struct cell) );
    for ( k = 0 ; k < nr_cells ; k++ )
        if ( lock_init( &cells[k].lock ) != 0 )
            error( "Failed to init spinlock." );
        
    /* Set the cell locations. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                c = &cells[ cell_getid( cdim , i , j , k ) ];
                c->loc[0] = i*h[0]; c->loc[1] = j*h[1]; c->loc[2] = k*h[2];
                c->h[0] = h[0]; c->h[1] = h[1]; c->h[2] = h[2];
                }
        
    /* Run through the particles and get their cell index. */
    ind = (int *)alloca( sizeof(int) * N );
    for ( k = 0 ; k < N ; k++ )  {
        ind[k] = cell_getid( cdim , parts[k].x[0]*ih[0] , parts[k].x[1]*ih[1] , parts[k].x[2]*ih[2] );
        cells[ ind[k] ].count += 1;
        }
        
    /* Sort the parts according to their cells. */
    parts_sort( parts , ind , N , 0 , nr_cells );
        
    /* Hook the cells up to the parts. */
    for ( finger = parts , k = 0 ; k < nr_cells ; k++ ) {
        c = &cells[ k ];
        c->parts = finger;
        finger = &finger[ c->count ];
        }
        
    /* Store eveything in the space. */
    s->h_min = h_min; s->h_max = h_max;
    s->dim[0] = dim[0]; s->dim[1] = dim[1]; s->dim[2] = dim[2];
    s->periodic = periodic;
    s->nr_parts = N;
    s->parts = parts;
    s->h[0] = h[0]; s->h[1] = h[1]; s->h[2] = h[2];
    s->ih[0] = ih[0]; s->ih[1] = ih[1]; s->ih[2] = ih[2];
    s->cdim[0] = cdim[0]; s->cdim[1] = cdim[1]; s->cdim[2] = cdim[2];
    s->cells = cells;
    s->nr_cells = nr_cells;
    s->tot_cells = nr_cells;
    if ( lock_init( &s->task_lock ) != 0 )
        error( "Failed to create task spin-lock." );
    
    /* Loop over the cells and split them. */
    for ( k = 0 ; k < nr_cells ; k++ )
        space_split( s , &cells[k] );
        
    }

