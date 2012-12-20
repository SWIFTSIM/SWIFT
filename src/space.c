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
#define error(s) { fprintf( stderr , "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize = space_subsize_default;

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
 * @brief Get the shift-id of the given pair of cells, swapping them
 *      if need be.
 *
 * @param s The space
 * @param ci Pointer to first #cell.
 * @param cj Pointer second #cell.
 * @param shift Vector from ci to cj.
 *
 * @return The shift ID and set shift, may or may not swap ci and cj.
 */
 
int space_getsid ( struct space *s , struct cell **ci , struct cell **cj , double *shift ) {

    int k, sid = 0;
    struct cell *temp;
    double dx[3];

    /* Get the relative distance between the pairs, wrapping. */
    for ( k = 0 ; k < 3 ; k++ ) {
        dx[k] = (*cj)->loc[k] - (*ci)->loc[k];
        if ( dx[k] < -s->dim[k]/2 )
            shift[k] = s->dim[k];
        else if ( dx[k] > s->dim[k]/2 )
            shift[k] = -s->dim[k];
        else
            shift[k] = 0.0;
        dx[k] += shift[k];
        }
        
    /* Get the sorting index. */
    for ( k = 0 ; k < 3 ; k++ )
        sid = 3*sid + ( (dx[k] < 0.0) ? 0 : ( (dx[k] > 0.0) ? 2 : 1 ) );

    /* Switch the cells around? */
    if ( runner_flip[sid] ) {
        temp = *ci; *ci = *cj; *cj = temp;
        for ( k = 0 ; k < 3 ; k++ )
            shift[k] = -shift[k];
        }
    sid = sortlistID[sid];
    
    /* Return the sort ID. */
    return sid;

    }


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
    float h, h_limit, h_max = 0.0f, dt_min = c->parts[0].dt, dt_max = dt_min;
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
            h = c->cparts[k].h;
            if ( h <= h_limit )
                count += 1;
            if ( h > h_max )
                h_max = h;
            if ( c->cparts[k].dt < dt_min )
                dt_min = c->cparts[k].dt;
            if ( c->cparts[k].dt > dt_max )
                dt_max = c->cparts[k].dt;
            }
        c->h_max = h_max;
        c->dt_min = dt_min;
        c->dt_max = dt_max;
            
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
 
int space_rebuild ( struct space *s , int force , double cell_max ) {

    float h_max = s->cell_min, h_min = s->parts[0].h;
    int i, j, k, cdim[3];
    struct cell *c;
    struct part *finger;
    struct cpart *cfinger;
    int *ind, changes = 0;
    // ticks tic;
    
    /* Run through the parts and get the current h_max. */
    // tic = getticks();
    for ( k = 0 ; k < s->nr_parts ; k++ ) {
        if ( s->parts[k].h > h_max )
            h_max = s->parts[k].h;
        else if ( s->parts[k].h < h_min )
            h_min = s->parts[k].h;
        }
    s->h_min = h_min;
    s->h_max = h_max;
    // printf( "space_rebuild: getting h_min and h_max took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    /* Get the new putative cell dimensions. */
    for ( k = 0 ; k < 3 ; k++ )
        cdim[k] = floor( s->dim[k] / fmax( h_max*space_stretch , cell_max ) );
        
    /* Do we need to re-build the upper-level cells? */
    // tic = getticks();
    if ( force || cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] || cdim[2] < s->cdim[2] ) {
    
        /* Free the old cells, if they were allocated. */
        if ( s->cells != NULL ) {
            for ( k = 0 ; k < s->nr_cells ; k++ ) {
                space_rebuild_recycle( s , &s->cells[k] );
                if ( s->cells[k].sort != NULL )
                    free( s->cells[k].sort );
                }
            free( s->cells );
            s->maxdepth = 0;
            }
            
        /* Set the new cell dimensions only if smaller. */
        for ( k = 0 ; k < 3 ; k++ ) {
            s->cdim[k] = cdim[k];
            s->h[k] = s->dim[k] / cdim[k];
            s->ih[k] = 1.0 / s->h[k];
            }

        /* Allocate the highest level of cells. */
        s->tot_cells = s->nr_cells = cdim[0] * cdim[1] * cdim[2];
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
    // printf( "space_rebuild: rebuilding upper-level cells took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
        
        
    /* Run through the particles and get their cell index. */
    // tic = getticks();
    if ( ( ind = (int *)malloc( sizeof(int) * s->nr_parts ) ) == NULL )
        error( "Failed to allocate temporary particle indices." );
    for ( k = 0 ; k < s->nr_cells ; k++ )
        s->cells[ k ].count = 0;
    for ( k = 0 ; k < s->nr_parts ; k++ )  {
        ind[k] = cell_getid( s->cdim , s->parts[k].x[0]*s->ih[0] , s->parts[k].x[1]*s->ih[1] , s->parts[k].x[2]*s->ih[2] );
        s->cells[ ind[k] ].count += 1;
        }
    // printf( "space_rebuild: getting particle indices took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );

    /* Sort the parts according to their cells. */
    // tic = getticks();
    parts_sort( s->parts , ind , s->nr_parts , 0 , s->nr_cells );    
    // printf( "space_rebuild: parts_sort took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    /* We no longer need the indices as of here. */
    free( ind );    

    /* Update the condensed particle data. */         
    // tic = getticks();
    #pragma omp parallel for schedule(static)
    for ( k = 0 ; k < s->nr_parts ; k++ ) {
        s->cparts[k].x[0] = s->parts[k].x[0];
        s->cparts[k].x[1] = s->parts[k].x[1];
        s->cparts[k].x[2] = s->parts[k].x[2];
        s->cparts[k].h = s->parts[k].h;
        s->cparts[k].dt = s->parts[k].dt;
        }
    // printf( "space_rebuild: creating condensed parts took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );

    /* Hook the cells up to the parts. */
    // tic = getticks();
    finger = s->parts;
    cfinger = s->cparts;
    for ( k = 0 ; k < s->nr_cells ; k++ ) {
        c = &s->cells[ k ];
        c->parts = finger;
        c->cparts = cfinger;
        finger = &finger[ c->count ];
        cfinger = &cfinger[ c->count ];
        }
    // printf( "space_rebuild: hooking up cells took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
        
        
    /* At this point, we have the upper-level cells, old or new. Now make
       sure that the parts in each cell are ok. */
    // tic = getticks();
    #pragma omp parallel for schedule(dynamic) shared(s) reduction(+:changes)
    for ( k = 0 ; k < s->nr_cells ; k++ )
        changes += space_rebuild_recurse( s , &s->cells[k] );
    // printf( "space_rebuild: space_rebuild_recurse took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
        
    /* Now that we have the cell structre, re-build the tasks. */
    // tic = getticks();
    if ( changes )
        space_maketasks( s , 1 );
    // printf( "space_rebuild: maketasks took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    /* Return the number of changes. */
    return changes;

    }


/**
 * @brief Sort the particles and condensed particles according to the given indices.
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
    
    /* If N is small enough, just do insert sort. */
    if ( N < 16 ) {
    
        for ( i = 1 ; i < N ; i++ )
            if ( ind[i] < ind[i-1] ) {
                temp_i = ind[i];
                temp_p = parts[j];
                for ( j = i ; j > 0 && ind[j-1] > temp_i ; j-- ) {
                    ind[j] = ind[j-1];
                    parts[j] = parts[j-1];
                    }
                ind[j] = temp_i;
                parts[j] = temp_p;
                }
    
        }
        
    /* Otherwise, recurse with Quicksort. */
    else {
    
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

        /* Try to recurse in parallel. */
        if ( N < 100 ) {

            /* Recurse on the left? */
            if ( j > 0 && pivot > min )
                parts_sort( parts , ind , j+1 , min , pivot );

            /* Recurse on the right? */
            if ( i < N && pivot+1 < max )
                parts_sort( &parts[i], &ind[i], N-i , pivot+1 , max );

            }

        else
        // #pragma omp parallel sections
        {

            /* Recurse on the left? */
            // #pragma omp section
            if ( j > 0 && pivot > min )
                parts_sort( parts , ind , j+1 , min , pivot );

            /* Recurse on the right? */
            // #pragma omp section
            if ( i < N && pivot+1 < max )
                parts_sort( &parts[i], &ind[i], N-i , pivot+1 , max );

            }
            
        }
    
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

    struct task *t;
    
    /* Lock the space. */
    lock_lock( &s->lock );
    
    /* Get the next free task. */
    t = &s->tasks[ s->nr_tasks ];
    
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
    
    /* Unock the space. */
    lock_unlock_blind( &s->lock );
    
    /* Return a pointer to the new task. */
    return t;

    }



/**
 * @brief Split tasks that may be too large.
 *
 * @param s The #space we are working in.
 */
 
void space_splittasks ( struct space *s ) {

    int j, k, sid, tid;
    struct cell *ci, *cj;
    double hi, hj, shift[3];
    struct task *t;
    float dt_max = s->dt_max;
    int pts[7][8] = { { -1 , 12 , 10 , 9 , 4 , 3 , 1 , 0 } ,
                      { -1 , -1 , 11 , 10 , 5 , 4 , 2 , 1 } ,
                      { -1 , -1 , -1 , 12 , 7 , 6 , 4 , 3 } , 
                      { -1 , -1 , -1 , -1 , 8 , 7 , 5 , 4 } ,
                      { -1 , -1 , -1 , -1 , -1 , 12 , 10 , 9 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , 11 , 10 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , -1 , 12 } };

    /* Loop through the tasks... */
    for ( tid = 0 ; tid < s->nr_tasks ; tid++ ) {
    
        /* Get a pointer on the task. */
        t = &s->tasks[tid];
        
        /* Self-interaction? */
        if ( t->type == task_type_self ) {
        
            /* Get a handle on the cell involved. */
            ci = t->ci;
            
            /* Ingore this task? */
            if ( ci->dt_min > dt_max ) {
                t->type = task_type_none;
                continue;
                }
            
            /* Is this cell even split? */
            if ( !ci->split )
                continue;
            
            /* Make a sub? */
            if ( space_dosub && ci->count < space_subsize ) {
            
                /* convert to a self-subtask. */
                t->type = task_type_sub;
                
                /* Wait for this tasks sorts, as we will now have pairwise
                   components in this sub. */
                space_addsorts( s , t , ci , NULL , -1 );
            
                }
                
            /* Otherwise, make tasks explicitly. */
            else {
            
                /* Take a step back (we're going to recycle the current task)... */
                tid -= 1;

                /* Add the self taks. */
                for ( k = 0 ; ci->progeny[k] == NULL ; k++ );
                t->ci = ci->progeny[k];
                for ( k += 1 ; k < 8 ; k++ )
                    if ( ci->progeny[k] != NULL )
                        space_addtask( s , task_type_self , task_subtype_density , 0 , 0 , ci->progeny[k] , NULL , NULL , 0 , NULL , 0 );
            
                /* Make a task for each pair of progeny. */
                for ( j = 0 ; j < 8 ; j++ )
                    if ( ci->progeny[j] != NULL )
                        for ( k = j + 1 ; k < 8 ; k++ )
                            if ( ci->progeny[k] != NULL )
                                space_addtask( s , task_type_pair , task_subtype_density , pts[j][k] , 0 , ci->progeny[j] , ci->progeny[k] , NULL , 0 , NULL , 0 );
                }
        
            }
    
        /* Pair interaction? */
        else if ( t->type == task_type_pair ) {
            
            /* Get a handle on the cells involved. */
            ci = t->ci;
            cj = t->cj;
            hi = fmax( ci->h[0] , fmax( ci->h[1] , ci->h[2] ) );
            hj = fmax( cj->h[0] , fmax( cj->h[1] , cj->h[2] ) );

            /* Ingore this task? */
            if ( ci->dt_min > dt_max && cj->dt_min > dt_max ) {
                t->type = task_type_none;
                continue;
                }
            
            /* Get the sort ID, use space_getsid and not t->flags
               to make sure we get ci and cj swapped if needed. */
            sid = space_getsid( s , &ci , &cj , shift );
                
            /* Should this task be split-up? */
            if ( ci->split && cj->split &&
                 ci->h_max*space_stretch < hi/2 && cj->h_max*space_stretch < hj/2 ) {
                 
                /* Replace by a single sub-task? */
                if ( space_dosub &&
                     ci->count < space_subsize && cj->count < space_subsize &&
                     sid != 0 && sid != 2 && sid != 6 && sid != 8 ) {
                
                    /* Make this task a sub task. */
                    t->type = task_type_sub;
                    t->flags = sid;
                    t->ci = ci; t->cj = cj;
                    
                    /* Create the sorts recursively. */
                    space_addsorts( s , t , ci , cj , sid );
                    
                    /* Don't go any further. */
                    continue;
                
                    }

                /* Take a step back (we're going to recycle the current task)... */
                tid -= 1;

                /* For each different sorting type... */
                switch ( sid ) {

                    case 0: /* (  1 ,  1 ,  1 ) */
                        t->ci = ci->progeny[7]; t->cj = cj->progeny[0]; t->flags = 0;
                        break;

                    case 1: /* (  1 ,  1 ,  0 ) */
                        t->ci = ci->progeny[6]; t->cj = cj->progeny[0]; t->flags = 1;
                        t = space_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        break;

                    case 2: /* (  1 ,  1 , -1 ) */
                        t->ci = ci->progeny[6]; t->cj = cj->progeny[1]; t->flags = 2;
                        break;

                    case 3: /* (  1 ,  0 ,  1 ) */
                        t->ci = ci->progeny[5]; t->cj = cj->progeny[0]; t->flags = 3;
                        t = space_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        break;

                    case 4: /* (  1 ,  0 ,  0 ) */
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[0]; t->flags = 4;
                        t = space_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[4] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[5] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[6] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[7] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        break;

                    case 5: /* (  1 ,  0 , -1 ) */
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[1]; t->flags = 5;
                        t = space_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        break;

                    case 6: /* (  1 , -1 ,  1 ) */
                        t->ci = ci->progeny[5]; t->cj = cj->progeny[2]; t->flags = 6;
                        break;

                    case 7: /* (  1 , -1 ,  0 ) */
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[3]; t->flags = 6;
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        break;

                    case 8: /* (  1 , -1 , -1 ) */
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[3]; t->flags = 8;
                        break;

                    case 9: /* (  0 ,  1 ,  1 ) */
                        t->ci = ci->progeny[3]; t->cj = cj->progeny[0]; t->flags = 9;
                        t = space_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        break;

                    case 10: /* (  0 ,  1 ,  0 ) */
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[0]; t->flags = 10;
                        t = space_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[2] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[3] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[2] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[6] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[3] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[7] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        break;

                    case 11: /* (  0 ,  1 , -1 ) */
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[1]; t->flags = 11;
                        t = space_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        break;

                    case 12: /* (  0 ,  0 ,  1 ) */
                        t->ci = ci->progeny[1]; t->cj = cj->progeny[0]; t->flags = 12;
                        t = space_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[1] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[3] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[1] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[5] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[1] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[3] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[5] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        t = space_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[7] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        break;

                    }

                } /* split this task? */
                
            /* Otherwise, if not spilt, stitch-up the sorting. */
            else {
            
                /* Create the sort for ci. */
                if ( ci->sorts[0] == NULL )
                    ci->sorts[0] = space_addtask( s , task_type_sort , 0 , 1 << sid , 0 , ci , NULL , NULL , 0 , NULL , 0 );
                ci->sorts[0]->flags |= (1 << sid);
                task_addunlock( ci->sorts[0] , t );
                
                /* Create the sort for cj. */
                if ( cj->sorts[0] == NULL )
                    cj->sorts[0] = space_addtask( s , task_type_sort , 0 , 1 << sid , 0 , cj , NULL , NULL , 0 , NULL , 0 );
                cj->sorts[0]->flags |= (1 << sid);
                task_addunlock( cj->sorts[0] , t );
                
                }
                
            } /* pair interaction? */
    
        } /* loop over all tasks. */
        
    }
    
    
/**
 * @brief Generate the sorts for a sub recursively.
 *
 * @param s The #space we are working in.
 */
 
void space_addsorts ( struct space *s , struct task *t , struct cell *ci , struct cell *cj , int sid ) {

    float h;
    double shift[3];
    int j, k;

    /* Get the cell dimensions. */
    h = fmin( ci->h[0] , fmin( ci->h[1] , ci->h[2] ) );
    
    /* Single-cell sub? */
    if ( cj == NULL ) {
    
        /* If there is further splitting, add the pairs recursively. */
        if ( ci->split ) {
        
            /* Recurse for each progeny. */
            for ( j = 0 ; j < 8 ; j++ )
                if ( ci->progeny[j] != NULL )
                    space_addsorts( s , t , ci->progeny[j] , NULL , -1 );

            /* Recurse for each pair of progeny. */
            for ( j = 0 ; j < 8 ; j++ )
                if ( ci->progeny[j] != NULL )
                    for ( k = j + 1 ; k < 8 ; k++ )
                        if ( ci->progeny[k] != NULL )
                            space_addsorts( s , t , ci->progeny[j] , ci->progeny[k] , -1 );

            }

        }
        
    /* Otherwise, it's a pair. */
    else {
        
        /* Get the sort ID if not specified. */
        // if ( sid < 0 )
            sid = space_getsid( s , &ci , &cj , shift );
        
        /* If there is no further splitting, add the sorts. */
        if ( !ci->split || !cj->split ||
             ci->h_max*2 >= h || cj->h_max*2 >= h ) {
            
            /* Create and add the sort for ci. */
            if ( ci->sorts[0] == NULL )
                ci->sorts[0] = space_addtask( s , task_type_sort , 0 , 1 << sid , 0 , ci , NULL , NULL , 0 , NULL , 0 );
            ci->sorts[0]->flags |= (1 << sid);
            task_addunlock( ci->sorts[0] , t );
            
            /* Create and add the sort for cj. */
            if ( cj->sorts[0] == NULL )
                cj->sorts[0] = space_addtask( s , task_type_sort , 0 , 1 << sid , 0 , cj , NULL , NULL , 0 , NULL , 0 );
            cj->sorts[0]->flags |= (1 << sid);
            task_addunlock( cj->sorts[0] , t );

            }

        /* Otherwise, recurse. */
        else {
                
            /* For each different sorting type... */
            switch ( sid ) {

                case 0: /* (  1 ,  1 ,  1 ) */
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 0 );
                    break;

                case 1: /* (  1 ,  1 ,  0 ) */
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[0] , 1 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[1] , 1 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 0 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 2 );
                    break;

                case 2: /* (  1 ,  1 , -1 ) */
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 2 );
                    break;

                case 3: /* (  1 ,  0 ,  1 ) */
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[0] , 3 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[2] , 3 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[2] , 0 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 6 );
                    break;

                case 4: /* (  1 ,  0 ,  0 ) */
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[0] , 4 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[0] , 5 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[0] , 7 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 8 );
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[1] , 3 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[1] , 4 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 6 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[1] , 7 );
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[2] , 1 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[2] , 2 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[2] , 4 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[2] , 5 );
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[3] , 0 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[3] , 1 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[3] , 3 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[3] , 4 );
                    break;

                case 5: /* (  1 ,  0 , -1 ) */
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[1] , 5 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[3] , 5 );
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[3] , 2 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 8 );
                    break;

                case 6: /* (  1 , -1 ,  1 ) */
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[2] , 6 );
                    break;

                case 7: /* (  1 , -1 ,  0 ) */
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[3] , 6 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[2] , 8 );
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[2] , 7 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[3] , 7 );
                    break;

                case 8: /* (  1 , -1 , -1 ) */
                    space_addsorts( s , t , ci->progeny[4] , cj->progeny[3] , 8 );
                    break;

                case 9: /* (  0 ,  1 ,  1 ) */
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[0] , 9 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[4] , 9 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[4] , 0 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 8 );
                    break;

                case 10: /* (  0 ,  1 ,  0 ) */
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[0] , 10 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[0] , 11 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[0] , 7 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 6 );
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[1] , 9 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[1] , 10 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 8 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[1] , 7 );
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[4] , 1 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[4] , 2 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[4] , 10 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[4] , 11 );
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[5] , 0 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[5] , 1 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[5] , 9 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[5] , 10 );
                    break;

                case 11: /* (  0 ,  1 , -1 ) */
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[1] , 11 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[5] , 11 );
                    space_addsorts( s , t , ci->progeny[2] , cj->progeny[5] , 2 );
                    space_addsorts( s , t , ci->progeny[6] , cj->progeny[1] , 6 );
                    break;

                case 12: /* (  0 ,  0 ,  1 ) */
                    space_addsorts( s , t , ci->progeny[1] , cj->progeny[0] , 12 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[0] , 11 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[0] , 5 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[0] , 2 );
                    space_addsorts( s , t , ci->progeny[1] , cj->progeny[2] , 9 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[2] , 12 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[2] , 8 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[2] , 5 );
                    space_addsorts( s , t , ci->progeny[1] , cj->progeny[4] , 3 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[4] , 6 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[4] , 12 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[4] , 11 );
                    space_addsorts( s , t , ci->progeny[1] , cj->progeny[6] , 0 );
                    space_addsorts( s , t , ci->progeny[3] , cj->progeny[6] , 3 );
                    space_addsorts( s , t , ci->progeny[5] , cj->progeny[6] , 9 );
                    space_addsorts( s , t , ci->progeny[7] , cj->progeny[6] , 12 );
                    break;

                } /* switch. */

            } /* recurse. */
            
        } /* it's a pair. */

    }
    
    
/**
 * @brief Fill the #space's task list.
 *
 * @param s The #space we are working in.
 * @param do_sort Flag to add sorting tasks to the list.
 */
 
void space_maketasks ( struct space *s , int do_sort ) {

    int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd, sid;
    int *cdim = s->cdim;
    struct task *t, *t2;
    struct cell *ci, *cj;
    // float dt_max = s->dt_max;
    int counts[task_type_count];

    /* Allocate the task-list, if needed. */
    if ( s->tasks == NULL || s->tasks_size < s->tot_cells * space_maxtaskspercell ) {
        if ( s->tasks != NULL )
            free( s->tasks );
        if ( s->tasks_ind != NULL )
            free( s->tasks_ind );
        s->tasks_size = s->tot_cells * space_maxtaskspercell;
        if ( posix_memalign( (void *)&s->tasks , 64 , sizeof(struct task) * s->tasks_size ) != 0 )
            error( "Failed to allocate task list." );
        if ( ( s->tasks_ind = (int *)malloc( sizeof(int) * s->tasks_size ) ) == NULL )
            error( "Failed to allocate task indices." );
        }
    s->nr_tasks = 0;
    
    /* Run through the highest level of cells and add pairs. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                cid = cell_getid( cdim , i , j , k );
                if ( s->cells[cid].count == 0 )
                    continue;
                ci = &s->cells[cid];
                if ( ci->count == 0 )
                    continue;
                // if ( ci->dt_min <= dt_max )
                    space_addtask( s , task_type_self , task_subtype_density , 0 , 0 , ci , NULL , NULL , 0 , NULL , 0 );
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
                            cj = &s->cells[cjd];
                            if ( cid >= cjd || cj->count == 0 /* ||
                                 ( ci->dt_min > dt_max && cj->dt_min > dt_max ) */ )
                                continue;
                            sid = sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ];
                            t = space_addtask( s , task_type_pair , task_subtype_density , sid , 0 , ci , cj , NULL , 0 , NULL , 0 );
                            }
                        }
                    }
                }

    /* Split the tasks. */
    space_splittasks( s );
    
    /* Make each sort depend on the sorts of its progeny. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( t->type == task_type_sort && t->ci->split )
            for ( j = 0 ; j < 8 ; j++ )
                if ( t->ci->progeny[j] != NULL && t->ci->progeny[j]->sorts[0] != NULL )
                    task_addunlock( t->ci->progeny[j]->sorts[0] , t );
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
    
    }
    
    

/**
 * @brief Split cells that contain too many particles.
 *
 * @param s The #space we are working in.
 * @param c The #cell under consideration.
 */
 
void space_split ( struct space *s , struct cell *c ) {

    int k, count;
    double h, h_limit, h_max = 0.0, dt_min = c->parts[0].dt, dt_max = dt_min;
    struct cell *temp;
    
    /* Check the depth. */
    if ( c->depth > s->maxdepth )
        s->maxdepth = c->depth;
    
    /* Set the minimum cutoff. */
    h_limit = fmin( c->h[0] , fmin( c->h[1] , c->h[2] ) ) / 2;
    
    /* Count the particles below that. */
    for ( count = 0 , k = 0 ; k < c->count ; k++ ) {
        h = c->cparts[k].h;
        if ( h <= h_limit )
            count += 1;
        if ( h > h_max )
            h_max = h;
        if ( c->cparts[k].dt < dt_min )
            dt_min = c->cparts[k].dt;
        if ( c->cparts[k].dt > dt_max )
            dt_max = c->cparts[k].dt;
        }
    c->h_max = h_max;
    c->dt_min = dt_min;
    c->dt_max = dt_max;
            
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

    /* Lock the space. */
    lock_lock( &s->lock );
    
    /* Clear the cell. */
    if ( lock_destroy( &c->lock ) != 0 )
        error( "Failed to destroy spinlock." );
        
    /* Clear this cell's sort arrays. */
    if ( c->sort != NULL )
        free( c->sort );
        
    /* Clear the cell data. */
    bzero( c , sizeof(struct cell) );
    
    /* Hook this cell into the buffer. */
    c->next = s->cells_new;
    s->cells_new = c;
    s->tot_cells -= 1;
    
    /* Unlock the space. */
    lock_unlock_blind( &s->lock );
    
    }


/**
 * @brief Get a new empty cell.
 *
 * @param s The #space.
 */
 
struct cell *space_getcell ( struct space *s ) {

    struct cell *c;
    int k;
    
    /* Lock the space. */
    lock_lock( &s->lock );
    
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
        
    /* Unlock the space. */
    lock_unlock_blind( &s->lock );
    
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


void space_init ( struct space *s , double dim[3] , struct part *parts , int N , int periodic , double h_max ) {

    /* Store eveything in the space. */
    s->dim[0] = dim[0]; s->dim[1] = dim[1]; s->dim[2] = dim[2];
    s->periodic = periodic;
    s->nr_parts = N;
    s->parts = parts;
    s->cell_min = h_max;
    
    /* Allocate the cparts array. */
    if ( posix_memalign( (void *)&s->cparts , 32 ,  N * sizeof(struct cpart) ) != 0 )
        error( "Failed to allocate cparts." );
        
    /* Init the space lock. */
    if ( lock_init( &s->lock ) != 0 )
        error( "Failed to create space spin-lock." );
    
    /* Build the cells and the tasks. */
    space_rebuild( s , 1 , h_max );
        
    }

