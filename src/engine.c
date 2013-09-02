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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <sched.h>

/* MPI headers. */
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* METIS headers. */
#ifdef HAVE_METIS
    #include <metis.h>
#endif

/* Local headers. */
#include "const.h"
#include "cycle.h"
#include "atomic.h"
#include "timers.h"
#include "const.h"
#include "vector.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "debug.h"
#include "space.h"
#include "cell.h"
#include "queue.h"
#include "scheduler.h"
#include "engine.h"
#include "runner.h"
#include "proxy.h"
#include "error.h"

#ifdef LEGACY_GADGET2_SPH
#include "runner_iact_legacy.h"
#else
#include "runner_iact.h"
#endif


/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )


/** The rank of the engine as a global variable (for messages). */
int engine_rank;


/**
 * @brief Redistribute the particles amongst the nodes accorind
 *      to their cell's node IDs.
 *
 * @param e The #engine.
 */
 
void engine_redistribute ( struct engine *e ) {

#ifdef WITH_MPI

    int i, j, k, cid;
    int nr_nodes = e->nr_nodes, nodeID = e->nodeID;
    struct space *s = e->s;
    int my_cells = 0;
    int *cdim = s->cdim;
    struct cell *cells = s->cells;
    int nr_cells = s->nr_cells;

    /* Start by sorting the particles according to their nodes and
       getting the counts. */
    int *counts, *dest;
    struct part *parts = s->parts;
    double ih[3];
    ih[0] = s->ih[0]; ih[1] = s->ih[1]; ih[2] = s->ih[2];
    if ( ( counts = (int *)malloc( sizeof(int) * nr_nodes * nr_nodes ) ) == NULL ||
         ( dest = (int *)malloc( sizeof(int) * s->nr_parts ) ) == NULL )
        error( "Failed to allocate count and dest buffers." );
    bzero( counts , sizeof(int) * nr_nodes * nr_nodes );
    for ( k = 0 ; k < s->nr_parts ; k++ ) {
        cid = cell_getid( cdim , parts[k].x[0]*ih[0] , parts[k].x[1]*ih[1] , parts[k].x[2]*ih[2] );
        dest[k] = cells[ cid ].nodeID;
        counts[ nodeID*nr_nodes + dest[k] ] += 1;
        }
    parts_sort( s->parts , s->xparts , dest , s->nr_parts , 0 , nr_nodes-1 );
    
    /* Get all the counts from all the nodes. */
    if ( MPI_Allreduce( MPI_IN_PLACE , counts , nr_nodes * nr_nodes , MPI_INT , MPI_SUM , MPI_COMM_WORLD ) != MPI_SUCCESS )
        error( "Failed to allreduce particle transfer counts." );
        
    /* Get the new number of parts for this node, be generous in allocating. */
    int nr_parts = 0;
    for ( k = 0 ; k < nr_nodes ; k++ )
        nr_parts += counts[ k*nr_nodes + nodeID ];
    struct part *parts_new;
    struct xpart *xparts_new, *xparts = s->xparts;
    if ( posix_memalign( (void **)&parts_new , part_align , sizeof(struct part) * nr_parts * 2 ) != 0 ||
         posix_memalign( (void **)&xparts_new , part_align , sizeof(struct xpart) * nr_parts * 2 ) != 0 )
        error( "Failed to allocate new part data." );
        
    /* Emit the sends and recvs for the particle data. */
    MPI_Request *reqs;
    if ( ( reqs = (MPI_Request *)malloc( sizeof(MPI_Request) * 4 * nr_nodes ) ) == NULL )
        error( "Failed to allocate MPI request list." );
    for ( k = 0 ; k < 4*nr_nodes ; k++ )
        reqs[k] = MPI_REQUEST_NULL;
    for ( i = 0 , j = 0 , k = 0 ; k < nr_nodes ; k++ ) {
        if ( k == nodeID && counts[ nodeID*nr_nodes + k ] > 0 ) {
            memcpy( &parts_new[j] , &parts[i] , sizeof(struct part) * counts[ k*nr_nodes + nodeID ] );
            memcpy( &xparts_new[j] , &xparts[i] , sizeof(struct xpart) * counts[ k*nr_nodes + nodeID ] );
            i += counts[ nodeID*nr_nodes + k ];
            j += counts[ k*nr_nodes + nodeID ];
            }
        if ( k != nodeID && counts[ nodeID*nr_nodes + k ] > 0 ) {
            if ( MPI_Isend( &parts[i] , sizeof(struct part) * counts[ nodeID*nr_nodes + k ] , MPI_BYTE , k , 2*(nodeID*nr_nodes + k) + 0 , MPI_COMM_WORLD , &reqs[4*k] ) != MPI_SUCCESS )
                error( "Failed to isend parts to node %i." , k );
            if ( MPI_Isend( &xparts[i] , sizeof(struct xpart) * counts[ nodeID*nr_nodes + k ] , MPI_BYTE , k , 2*(nodeID*nr_nodes + k) + 1 , MPI_COMM_WORLD , &reqs[4*k+1] ) != MPI_SUCCESS )
                error( "Failed to isend xparts to node %i." , k );
            i += counts[ nodeID*nr_nodes + k ];
            }
        if ( k != nodeID && counts[ k*nr_nodes + nodeID ] > 0 ) {
            if ( MPI_Irecv( &parts_new[j] , sizeof(struct part) * counts[ k*nr_nodes + nodeID ] , MPI_BYTE , k , 2*(k*nr_nodes + nodeID) + 0 , MPI_COMM_WORLD , &reqs[4*k+2] ) != MPI_SUCCESS )
                error( "Failed to emit irecv of parts from node %i." , k );
            if ( MPI_Irecv( &xparts_new[j] , sizeof(struct xpart) * counts[ k*nr_nodes + nodeID ] , MPI_BYTE , k , 2*(k*nr_nodes + nodeID) + 1 , MPI_COMM_WORLD , &reqs[4*k+3] ) != MPI_SUCCESS )
                error( "Failed to emit irecv of parts from node %i." , k );
            j += counts[ k*nr_nodes + nodeID ];
            }
        }
        
    /* Wait for all the sends and recvs to tumble in. */
    if ( MPI_Waitall( 4*nr_nodes , reqs , MPI_STATUSES_IGNORE ) != MPI_SUCCESS )
        error( "Failed during waitall for part data." );
        
    /* Verify that all parts are in the right place. */
    /* for ( k = 0 ; k < nr_parts ; k++ ) {
        cid = cell_getid( cdim , parts_new[k].x[0]*ih[0] , parts_new[k].x[1]*ih[1] , parts_new[k].x[2]*ih[2] );
        if ( cells[ cid ].nodeID != nodeID )
            error( "Received particle (%i) that does not belong here (nodeID=%i)." , k , cells[ cid ].nodeID );
        } */
        
    /* Set the new part data, free the old. */
    free( parts );
    free( xparts );
    s->parts = parts_new;
    s->xparts = xparts_new;
    s->nr_parts = nr_parts;
    s->size_parts = 2*nr_parts;
    
    /* Be verbose about what just happened. */
    for ( k = 0 ; k < nr_cells ; k++ )
        if ( cells[k].nodeID == nodeID )
            my_cells += 1;
    message( "node %i now has %i parts in %i cells." , nodeID , nr_parts , my_cells );
    
    /* Clean up other stuff. */
    free( reqs );
    free( counts );
    free( dest );
        
#else
    error( "SWIFT was not compiled with MPI and METIS support." );
#endif

    }


/**
 * @breif Repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */
 
void engine_repartition ( struct engine *e ) {

#if defined(WITH_MPI) && defined(HAVE_METIS)

    int i, j, k, l, cid, cjd, ii, jj, kk, res, w;
    idx_t *inds;
    idx_t *weights_v, *weights_e;
    struct space *s = e->s;
    int nr_cells = s->nr_cells, my_cells = 0;
    struct cell *cells = s->cells;
    int ind[3], *cdim = s->cdim;
    struct task *t, *tasks = e->sched.tasks;
    struct cell *ci, *cj;
    int nr_nodes = e->nr_nodes, nodeID = e->nodeID, *nodeIDs;
    float wscale = 0.01;
    
    /* Clear the repartition flag. */
    e->forcerepart = 0;
    
    /* Allocate the inds and weights. */
    if ( ( inds = (idx_t *)malloc( sizeof(idx_t) * 26*nr_cells ) ) == NULL ||
         ( weights_v = (idx_t *)malloc( sizeof(idx_t) * nr_cells ) ) == NULL ||
         ( weights_e = (idx_t *)malloc( sizeof(idx_t) * 26*nr_cells ) ) == NULL ||
         ( nodeIDs = (idx_t *)malloc( sizeof(idx_t) * nr_cells ) ) == NULL )
        error( "Failed to allocate inds and weights arrays." );
        
    /* Fill the inds array. */
    for ( cid = 0 ; cid < nr_cells ; cid++ ) {
        ind[0] = cells[cid].loc[0] / s->cells[cid].h[0] + 0.5;
        ind[1] = cells[cid].loc[1] / s->cells[cid].h[1] + 0.5;
        ind[2] = cells[cid].loc[2] / s->cells[cid].h[2] + 0.5;
        l = 0;
        for ( i = -1 ; i <= 1 ; i++ ) {
            ii = ind[0] + i;
            if ( ii < 0 ) ii += cdim[0];
            else if ( ii >= cdim[0] ) ii -= cdim[0];
            for ( j = -1 ; j <= 1 ; j++ ) {
                jj = ind[1] + j;
                if ( jj < 0 ) jj += cdim[1];
                else if ( jj >= cdim[1] ) jj -= cdim[1];
                for ( k = -1 ; k <= 1 ; k++ ) {
                    kk = ind[2] + k;
                    if ( kk < 0 ) kk += cdim[2];
                    else if ( kk >= cdim[2] ) kk -= cdim[2];
                    if ( i || j || k ) {
                        inds[ cid*26 + l ] = cell_getid( cdim , ii , jj , kk );
                        l += 1;
                        }
                    }
                }
            }
        }
        
    /* Init the weights arrays. */
    /* bzero( weights_e , sizeof(idx_t) * 26*nr_cells );
    bzero( weights_v , sizeof(idx_t) * nr_cells ); */
    for ( k = 0 ; k < 26*nr_cells ; k++ )
        weights_e[k] = 1;
    for ( k = 0 ; k < nr_cells ; k++ )
        weights_v[k] = 1;
    
    /* Loop over the tasks... */
    for ( j = 0 ; j < e->sched.nr_tasks ; j++ ) {
    
        /* Get a pointer to the kth task. */
        t = &tasks[j];
        
        /* Skip un-interesting tasks. */
        if ( t->type != task_type_self &&
             t->type != task_type_pair &&
             t->type != task_type_sub &&
             t->type != task_type_ghost &&
             t->type != task_type_kick1 &&
             t->type != task_type_kick2 )
            continue;
            
        /* Get the task weight. */
        w = ( t->toc - t->tic ) * wscale;
        if ( w < 0 )
            error( "Bad task weight (%i)." , w );
        
        /* Get the top-level cells involved. */
        for ( ci = t->ci ; ci->parent != NULL ; ci = ci->parent );
        if ( t->cj != NULL )
            for ( cj = t->cj ; cj->parent != NULL ; cj = cj->parent );
        else
            cj = NULL;
            
        /* Get the cell IDs. */
        cid = ci - cells;
            
        /* Different weights for different tasks. */
        if ( t->type == task_type_ghost ||
             t->type == task_type_kick1 || 
             t->type == task_type_kick2 ) {
             
            /* Particle updates add only to vertex weight. */
            weights_v[cid] += w;
            
            }
        
        /* Self interaction? */     
        else if ( ( t->type == task_type_self && ci->nodeID == nodeID ) ||
             ( t->type == task_type_sub && cj == NULL && ci->nodeID == nodeID ) ) {
        
            /* Self interactions add only to vertex weight. */
            weights_v[cid] += w;
            
            }
            
        /* Pair? */
        else if ( t->type == task_type_pair ||
                  ( t->type == task_type_sub && cj != NULL ) ) {
                  
            /* In-cell pair? */
            if ( ci == cj ) {
            
                /* Add weight to vertex for ci. */
                weights_v[cid] += w;
            
                }
                
            /* Distinct cells with local ci? */
            else if ( ci->nodeID == nodeID ) {
            
                /* Index of the jth cell. */
                cjd = cj - cells;
                
                /* Add half of weight to each cell. */
                if ( ci->nodeID == nodeID )
                    weights_v[cid] += 0.5 * w;
                if ( cj->nodeID == nodeID )
                    weights_v[cjd] += 0.5 * w;
                    
                /* Add Weight to edge. */
                for ( k = 26*cid ; inds[k] != cjd ; k++ );
                weights_e[ k ] += w;
                for ( k = 26*cjd ; inds[k] != cid ; k++ );
                weights_e[ k ] += w;
            
                }
                  
            }
    
        }
        
    /* Merge the weights arrays accross all nodes. */
    if ( ( res = MPI_Reduce( ( nodeID == 0 ) ? MPI_IN_PLACE : weights_v , weights_v , nr_cells , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ) ) != MPI_SUCCESS ) {
        char buff[ MPI_MAX_ERROR_STRING ];
        MPI_Error_string( res , buff , &i );
        error( "Failed to allreduce vertex weights (%s)." , buff );
        }
    if ( MPI_Reduce( ( nodeID == 0 ) ? MPI_IN_PLACE : weights_e , weights_e , 26*nr_cells , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ) != MPI_SUCCESS )
        error( "Failed to allreduce edge weights." );
        
    /* As of here, only one node needs to compute the partition. */
    if ( nodeID == 0 ) {
    
        /* Check that the edge weights are fully symmetric. */
        /* for ( cid = 0 ; cid < nr_cells ; cid++ )
            for ( k = 0 ; k < 26 ; k++ ) {
                cjd = inds[ cid*26 + k ];
                for ( j = 26*cjd ; inds[j] != cid ; j++ );
                if ( weights_e[ cid*26+k ] != weights_e[ j ] )
                    error( "Unsymmetric edge weights detected (%i vs %i)." , weights_e[ cid*26+k ] , weights_e[ j ] );
                } */
        /* int w_min = weights_e[0], w_max = weights_e[0];
        for ( k = 1 ; k < 26*nr_cells ; k++ )
            if ( weights_e[k] < w_min )
                w_min = weights_e[k];
            else if ( weights_e[k] > w_max )
                w_max = weights_e[k];
        message( "edge weights in [ %i , %i ]." , w_min , w_max );
        w_min = weights_e[0], w_max = weights_e[0];
        for ( k = 1 ; k < nr_cells ; k++ )
            if ( weights_v[k] < w_min )
                w_min = weights_v[k];
            else if ( weights_v[k] > w_max )
                w_max = weights_v[k];
        message( "vertex weights in [ %i , %i ]." , w_min , w_max ); */
                
        /* Make sure there are no zero weights. */
        for ( k = 0 ; k < 26*nr_cells ; k++ )
            if ( weights_e[k] == 0 )
                weights_e[k] = 1;
    
        /* Allocate and fill the connection array. */
        idx_t *offsets;
        if ( ( offsets = (idx_t *)malloc( sizeof(idx_t) * (nr_cells + 1) ) ) == NULL )
            error( "Failed to allocate offsets buffer." );
        offsets[0] = 0;
        for ( k = 0 ; k < nr_cells ; k++ )
            offsets[k+1] = offsets[k] + 26;
            
        /* Set the METIS options. */
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions( options );
        options[ METIS_OPTION_OBJTYPE ] = METIS_OBJTYPE_CUT;
        options[ METIS_OPTION_NUMBERING ] = 0;
        options[ METIS_OPTION_CONTIG ] = 1;
        
        /* Set the initial partition, although this is probably ignored. */
        for ( k = 0 ; k < nr_cells ; k++ )
            nodeIDs[k] = cells[k].nodeID;
            
        /* Call METIS. */
        int one = 1;
        idx_t objval;
        if ( METIS_PartGraphKway( &nr_cells , &one , offsets , inds , weights_v , NULL , weights_e , &nr_nodes , NULL , NULL , options , &objval , nodeIDs ) != METIS_OK )
            error( "Call to METIS_PartGraphKway failed." );
    
        }
        
    /* Broadcast the result of the partition. */
    if ( MPI_Bcast( nodeIDs , nr_cells , MPI_INT , 0 , MPI_COMM_WORLD ) != MPI_SUCCESS )
        error( "Failed to bcast the node IDs." );
        
    /* Set the cell nodeIDs and clear any non-local parts. */
    for ( k = 0 ; k < nr_cells ; k++ ) {
        cells[k].nodeID = nodeIDs[k];
        if ( nodeIDs[k] == nodeID )
            my_cells += 1;
        }
        
    /* Clean up. */
    free( inds );
    free( weights_v );
    free( weights_e );
    free( nodeIDs );
        
    /* Now comes the tricky part: Exchange particles between all nodes.
       This is done in two steps, first allreducing a matrix of 
       how many particles go from where to where, then re-allocating
       the parts array, and emiting the sends and receives.
       Finally, the space, tasks, and proxies need to be rebuilt. */
       
    /* Redistribute the particles between the nodes. */
    engine_redistribute( e );
        
    /* Make the proxies. */
    engine_makeproxies( e );
        
    /* Tell the engine it should re-build whenever possible */
    e->forcerebuild = 1;
    
#else
    error( "SWIFT was not compiled with MPI and METIS support." );
#endif

    }


/**
 * @brief Add send tasks to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param t_xv The send_xv #task, if it has already been created.
 * @param t_rho The send_rho #task, if it has already been created.
 */

void engine_addtasks_send ( struct engine *e , struct cell *ci , struct cell *cj ) {

    int k;

    /* Check if any of the density tasks are for the target node. */
    for ( k = 0 ; k < ci->nr_density ; k++ )
        if ( ci->density[k]->ci->nodeID == cj->nodeID ||
             ( ci->density[k]->cj != NULL && ci->density[k]->cj->nodeID == cj->nodeID ) )
            break;

    /* If so, attach send tasks. */
    if ( k < ci->nr_density ) {

        /* Create the tasks. */
        struct task *t_xv = scheduler_addtask( &e->sched , task_type_send_xv , task_subtype_none , 2*ci->tag , 0 , ci , cj , 0 );
        struct task *t_rho = scheduler_addtask( &e->sched , task_type_send_rho , task_subtype_none , 2*ci->tag + 1 , 0 , ci , cj , 0 );

        /* The send_rho task depends on the cell's ghost task. */
        task_addunlock( ci->ghost , t_rho );

        /* The send_rho task should unlock the super-cell's kick2 task. */
        task_addunlock( t_rho , ci->super->kick2 );

        /* The send_xv task should unlock the super-cell's ghost task. */
        task_addunlock( t_xv , ci->super->ghost );

        }
        
    /* Recurse? */
    else if ( ci->split )
        for ( k = 0 ; k < 8 ; k++ )
            if ( ci->progeny[k] != NULL )
                engine_addtasks_send( e , ci->progeny[k] , cj );

    }


/**
 * @brief Add recv tasks to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param t_xv The recv_xv #task, if it has already been created.
 * @param t_rho The recv_rho #task, if it has already been created.
 */

void engine_addtasks_recv ( struct engine *e , struct cell *c , struct task *t_xv , struct task *t_rho ) {

    int k;

    /* Do we need to construct a recv task? */
    if ( t_xv != NULL || c->nr_density > 0 ) {
    
        /* Create the tasks. */
        c->recv_xv = scheduler_addtask( &e->sched , task_type_recv_xv , task_subtype_none , 2*c->tag , 0 , c , NULL , 0 );
        c->recv_rho = scheduler_addtask( &e->sched , task_type_recv_rho , task_subtype_none , 2*c->tag + 1 , 0 , c , NULL , 0 );
        
        /* If there has been a higher-up recv task, then these tasks
           are implicit and depend on the higher-up task. */
        if ( t_xv != NULL ) {
            task_addunlock( c->parent->recv_xv , c->recv_xv );
            task_addunlock( c->parent->recv_rho , c->recv_rho );
            c->recv_xv->implicit = 1;
            c->recv_rho->implicit = 1;
            }
        else {
            t_xv = c->recv_xv;
            t_rho = c->recv_rho;
            }
        
        /* Add dependencies if there are density/force tasks. */
        for ( k = 0 ; k < c->nr_density ; k++ ) {
            task_addunlock( c->recv_xv , c->density[k] );
            task_addunlock( c->density[k] , t_rho );
            }
        for ( k = 0 ; k < c->nr_force ; k++ )
            task_addunlock( c->recv_rho , c->force[k] );
        if ( c->sorts != NULL )
            task_addunlock( c->recv_xv , c->sorts );
            
        }
        
    /* Recurse? */
    if ( c->split )
        for ( k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL )
                engine_addtasks_recv( e , c->progeny[k] , t_xv , t_rho );

    }


/**
 * @brief Exchange cell structures with other nodes.
 *
 * @param e The #engine.
 */
 
void engine_exchange_cells ( struct engine *e ) {

#ifdef WITH_MPI

    int j, k, pid, count = 0;
    struct pcell *pcells;
    struct space *s = e->s;
    struct cell *cells = s->cells;
    int nr_cells = s->nr_cells;
    int nr_proxies = e->nr_proxies;
    int offset[ nr_cells ];
    MPI_Request reqs_in[ engine_maxproxies ];
    MPI_Request reqs_out[ engine_maxproxies ];
    MPI_Status status;
    struct part *parts = &s->parts[ s->nr_parts ];
    
    /* Run through the cells and get the size of the ones that will be sent off. */
    for ( k = 0 ; k < nr_cells ; k++ ) {
        offset[k] = count;
        if ( cells[k].sendto )
            count += ( cells[k].pcell_size = cell_getsize( &cells[k] ) );
        }
        
    /* Allocate the pcells. */
    if ( ( pcells = (struct pcell *)malloc( sizeof(struct pcell) * count ) ) == NULL )
        error( "Failed to allocate pcell buffer." );
        
    /* Pack the cells. */
    for ( k = 0 ; k < nr_cells ; k++ )
        if ( cells[k].sendto ) {
            cell_pack( &cells[k] , &pcells[ offset[k] ] );
            cells[k].pcell = &pcells[ offset[k] ];
            }

    /* Launch the proxies. */
    for ( k = 0 ; k < nr_proxies ; k++ ) {
        proxy_cells_exch1( &e->proxies[k] );
        reqs_in[k] = e->proxies[k].req_cells_count_in;
        reqs_out[k] = e->proxies[k].req_cells_count_out;
        }
        
    /* Wait for each count to come in and start the recv. */
    for ( k = 0 ; k < nr_proxies ; k++ ) {
        if ( MPI_Waitany( nr_proxies , reqs_in , &pid , &status ) != MPI_SUCCESS ||
             pid == MPI_UNDEFINED )
            error( "MPI_Waitany failed." );
        // message( "request from proxy %i has arrived." , pid );
        proxy_cells_exch2( &e->proxies[pid] );
        }
        
    /* Wait for all the sends to have finnished too. */
    if ( MPI_Waitall( nr_proxies , reqs_out , MPI_STATUSES_IGNORE ) != MPI_SUCCESS )
        error( "MPI_Waitall on sends failed." );
        
    /* Set the requests for the cells. */
    for ( k = 0 ; k < nr_proxies ; k++ ) {
        reqs_in[k] = e->proxies[k].req_cells_in;
        reqs_out[k] = e->proxies[k].req_cells_out;
        }
    
    /* Wait for each pcell array to come in from the proxies. */
    for ( k = 0 ; k < nr_proxies ; k++ ) {
        if ( MPI_Waitany( nr_proxies , reqs_in , &pid , &status ) != MPI_SUCCESS ||
             pid == MPI_UNDEFINED )
            error( "MPI_Waitany failed." );
        // message( "cell data from proxy %i has arrived." , pid );
        for ( count = 0 , j = 0 ; j < e->proxies[pid].nr_cells_in ; j++ )
            count += cell_unpack( &e->proxies[pid].pcells_in[count] , e->proxies[pid].cells_in[j] , e->s );
        }
        
    /* Wait for all the sends to have finnished too. */
    if ( MPI_Waitall( nr_proxies , reqs_out , MPI_STATUSES_IGNORE ) != MPI_SUCCESS )
        error( "MPI_Waitall on sends failed." );
        
    /* Count the number of particles we need to import and re-allocate
       the buffer if needed. */
    for ( count = 0 , k = 0 ; k < nr_proxies ; k++ )
        for ( j = 0 ; j < e->proxies[k].nr_cells_in ; j++ )
            count += e->proxies[k].cells_in[j]->count;
    if ( count > s->size_parts_foreign ) {
        if ( s->parts_foreign != NULL )
            free( s->parts_foreign );
        s->size_parts_foreign = 1.1 * count;
        if ( posix_memalign( (void **)&s->parts_foreign , part_align , sizeof(struct part) * s->size_parts_foreign ) != 0 )
            error( "Failed to allocate foreign part data." );
        }
        
    /* Unpack the cells and link to the particle data. */
    parts = s->parts_foreign;
    for ( k = 0 ; k < nr_proxies ; k++ ) {
        for ( count = 0 , j = 0 ; j < e->proxies[k].nr_cells_in ; j++ ) {
            count += cell_link( e->proxies[k].cells_in[j] , parts );
            parts = &parts[ e->proxies[k].cells_in[j]->count ];
            }
        }
    s->nr_parts_foreign = parts - s->parts_foreign;
        
    /* Is the parts buffer large enough? */
    if ( s->nr_parts_foreign > s->size_parts_foreign )
        error( "Foreign parts buffer too small." );
        
    /* Free the pcell buffer. */
    free( pcells );
    
#else
    error( "SWIFT was not compiled with MPI support." );
#endif

    }


/**
 * @brief Exchange straying parts with other nodes.
 *
 * @param e The #engine.
 * @param parts An array of straying parts.
 * @param xparts The corresponding xparts.
 * @param ind The ID of the foreign #cell.
 * @param N The number of stray parts.
 *
 * @return The number of arrived parts copied to parts and xparts.
 */
 
int engine_exchange_strays ( struct engine *e , struct part *parts , struct xpart *xparts , int *ind , int N ) {

#ifdef WITH_MPI

    int k, pid, count = 0, nr_in = 0, nr_out = 0;
    MPI_Request reqs_in[ 2*engine_maxproxies ];
    MPI_Request reqs_out[ 2*engine_maxproxies ];
    MPI_Status status;
    struct proxy *p;

    /* Re-set the proxies. */
    for ( k = 0 ; k < e->nr_proxies ; k++ )
        e->proxies[k].nr_parts_out = 0;
    
    /* Put the parts into the corresponding proxies. */
    for ( k = 0 ; k < N ; k++ ) {
        pid = e->proxy_ind[ e->s->cells[ ind[k] ].nodeID ];
        if ( pid < 0 )
            error( "Do not have a proxy for the requested nodeID." );
        proxy_parts_load( &e->proxies[pid] , &parts[k] , &xparts[k] , 1 );
        }
    
    /* Launch the proxies. */
    for ( k = 0 ; k < e->nr_proxies ; k++ ) {
        proxy_parts_exch1( &e->proxies[k] );
        reqs_in[k] = e->proxies[k].req_parts_count_in;
        reqs_out[k] = e->proxies[k].req_parts_count_out;
        }
        
    /* Wait for each count to come in and start the recv. */
    for ( k = 0 ; k < e->nr_proxies ; k++ ) {
        if ( MPI_Waitany( e->nr_proxies , reqs_in , &pid , &status ) != MPI_SUCCESS ||
             pid == MPI_UNDEFINED )
            error( "MPI_Waitany failed." );
        // message( "request from proxy %i has arrived." , pid );
        proxy_parts_exch2( &e->proxies[pid] );
        }
        
    /* Wait for all the sends to have finnished too. */
    if ( MPI_Waitall( e->nr_proxies , reqs_out , MPI_STATUSES_IGNORE ) != MPI_SUCCESS )
        error( "MPI_Waitall on sends failed." );
        
    /* Return the number of harvested parts. */
    /* Set the requests for the particle data. */
    for ( k = 0 ; k < e->nr_proxies ; k++ ) {
        if ( e->proxies[k].nr_parts_in > 0 ) {
            reqs_in[2*k] = e->proxies[k].req_parts_in;
            reqs_in[2*k+1] = e->proxies[k].req_xparts_in;
            nr_in += 1;
            }
        else
            reqs_in[2*k] = reqs_in[2*k+1] = MPI_REQUEST_NULL;
        if ( e->proxies[k].nr_parts_out > 0 ) {
            reqs_out[2*k] = e->proxies[k].req_parts_out;
            reqs_out[2*k+1] = e->proxies[k].req_xparts_out;
            nr_out += 1;
            }
        else
            reqs_out[2*k] = reqs_out[2*k+1] = MPI_REQUEST_NULL;
        }
    
    /* Wait for each part array to come in and collect the new
       parts from the proxies. */
    for ( k = 0 ; k < 2*(nr_in + nr_out) ; k++ ) {
        if ( MPI_Waitany( 2*e->nr_proxies , reqs_in , &pid , &status ) != MPI_SUCCESS )
            error( "MPI_Waitany failed." );
        if ( pid == MPI_UNDEFINED )
            break;
        // message( "request from proxy %i has arrived." , pid );
        if ( reqs_in[pid & ~1] == MPI_REQUEST_NULL &&
             reqs_in[pid | 1 ] == MPI_REQUEST_NULL ) {
            p = &e->proxies[pid/2];
            memcpy( &parts[count] , p->parts_in , sizeof(struct part) * p->nr_parts_in );
            memcpy( &xparts[count] , p->xparts_in , sizeof(struct xpart) * p->nr_parts_in );
            count += p->nr_parts_in;
            /* for ( int k = 0 ; k < p->nr_parts_in ; k++ )
                message( "received particle %lli, x=[%.3e %.3e %.3e], h=%.3e, from node %i." ,
                    p->parts_in[k].id , p->parts_in[k].x[0] , p->parts_in[k].x[1] , p->parts_in[k].x[2] ,
                    p->parts_in[k].h , p->nodeID ); */
            }
        }
    
    /* Wait for all the sends to have finnished too. */
    if ( nr_out > 0 )
        if ( MPI_Waitall( 2*e->nr_proxies , reqs_out , MPI_STATUSES_IGNORE ) != MPI_SUCCESS )
            error( "MPI_Waitall on sends failed." );
        
    /* Return the number of harvested parts. */
    return count;
    
#else
    error( "SWIFT was not compiled with MPI support." );
    return 0;
#endif

    }


/**
 * @brief Fill the #space's task list.
 *
 * @param s The #space we are working in.
 */
 
void engine_maketasks ( struct engine *e ) {

    struct space *s = e->s;
    struct scheduler *sched = &e->sched;
    int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd, sid;
    int *cdim = s->cdim;
    struct task *t, *t2;
    struct cell *ci, *cj;

    /* Re-set the scheduler. */
    scheduler_reset( sched , s->tot_cells * engine_maxtaskspercell );
    
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
                if ( ci->nodeID == e->nodeID )
                    scheduler_addtask( sched , task_type_self , task_subtype_density , 0 , 0 , ci , NULL , 0 );
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
                            if ( cid >= cjd || cj->count == 0 || 
                                 ( ci->nodeID != e->nodeID && cj->nodeID != e->nodeID ) )
                                continue;
                            sid = sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ];
                            t = scheduler_addtask( sched , task_type_pair , task_subtype_density , sid , 0 , ci , cj , 1 );
                            }
                        }
                    }
                }

    /* Split the tasks. */
    scheduler_splittasks( sched );
    
    /* Count the number of tasks associated with each cell and
       store the density tasks in each cell, and make each sort
       depend on the sorts of its progeny. */
    // #pragma omp parallel for private(t,j)
    for ( k = 0 ; k < sched->nr_tasks ; k++ ) {
        t = &sched->tasks[k];
        if ( t->skip )
            continue;
        if ( t->type == task_type_sort && t->ci->split )
            for ( j = 0 ; j < 8 ; j++ )
                if ( t->ci->progeny[j] != NULL && t->ci->progeny[j]->sorts != NULL ) {
                    t->ci->progeny[j]->sorts->skip = 0;
                    task_addunlock( t->ci->progeny[j]->sorts , t );
                    }
        if ( t->type == task_type_self ) {
            atomic_inc( &t->ci->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                }
            }
        else if ( t->type == task_type_pair ) {
            atomic_inc( &t->ci->nr_tasks );
            atomic_inc( &t->cj->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                t->cj->density[ atomic_inc( &t->cj->nr_density ) ] = t;
                }
            }
        else if ( t->type == task_type_sub ) {
            atomic_inc( &t->ci->nr_tasks );
            if ( t->cj != NULL )
                atomic_inc( &t->cj->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                if ( t->cj != NULL )
                    t->cj->density[ atomic_inc( &t->cj->nr_density ) ] = t;
                }
            }
        }
        
    /* Append a ghost task to each cell. */
    if ( e->policy & engine_policy_fixdt )
        space_map_cells_pre( s , 1 , &scheduler_map_mkghosts_nokick1 , sched );
    else
        space_map_cells_pre( s , 1 , &scheduler_map_mkghosts , sched );
    
    /* Run through the tasks and make force tasks for each density task.
       Each force task depends on the cell ghosts and unlocks the kick2 task
       of its super-cell. */
    kk = sched->nr_tasks;
    // #pragma omp parallel for private(t,t2)
    for ( k = 0 ; k < kk ; k++ ) {
    
        /* Get a pointer to the task. */
        t = &sched->tasks[k];
        
        /* Skip? */
        if ( t->skip )
            continue;
        
        /* Self-interaction? */
        if ( t->type == task_type_self && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            t2 = scheduler_addtask( sched , task_type_self , task_subtype_force , 0 , 0 , t->ci , NULL , 0 );
            task_addunlock( t->ci->ghost , t2 );
            task_addunlock( t2 , t->ci->super->kick2 );
            t->ci->force[ atomic_inc( &t->ci->nr_force ) ] = t2;
            }
            
        /* Otherwise, pair interaction? */
        else if ( t->type == task_type_pair && t->subtype == task_subtype_density ) {
            t2 = scheduler_addtask( sched , task_type_pair , task_subtype_force , 0 , 0 , t->ci , t->cj , 0 );
            if ( t->ci->nodeID == e->nodeID ) {
                task_addunlock( t->ci->ghost , t2 );
                task_addunlock( t , t->ci->super->ghost );
                task_addunlock( t2 , t->ci->super->kick2 );
                }
            if ( t->cj->nodeID == e->nodeID ) {
                task_addunlock( t->cj->ghost , t2 );
                if ( t->ci->super != t->cj->super ) {
                    task_addunlock( t , t->cj->super->ghost );
                    task_addunlock( t2 , t->cj->super->kick2 );
                    }
                }
            t->ci->force[ atomic_inc( &t->ci->nr_force ) ] = t2;
            t->cj->force[ atomic_inc( &t->cj->nr_force ) ] = t2;
            }
    
        /* Otherwise, sub interaction? */
        else if ( t->type == task_type_sub && t->subtype == task_subtype_density ) {
            t2 = scheduler_addtask( sched , task_type_sub , task_subtype_force , t->flags , 0 , t->ci , t->cj , 0 );
            if ( t->ci->nodeID == e->nodeID ) {
                task_addunlock( t , t->ci->super->ghost );
                task_addunlock( t->ci->ghost , t2 );
                task_addunlock( t2 , t->ci->super->kick2 );
                }
            if ( t->cj != NULL && t->cj->nodeID == e->nodeID ) {
                task_addunlock( t->cj->ghost , t2 );
                if ( t->ci->super != t->cj->super ) {
                    task_addunlock( t , t->cj->super->ghost );
                    task_addunlock( t2 , t->cj->super->kick2 );
                    }
                }
            t->ci->force[ atomic_inc( &t->ci->nr_force ) ] = t2;
            if ( t->cj != NULL )
                t->cj->force[ atomic_inc( &t->cj->nr_force ) ] = t2;
            }
            
        }
        
    /* Add the communication tasks if MPI is being used. */
    #ifdef WITH_MPI
        
        /* Loop over the proxies. */
        for ( int pid = 0 ; pid < e->nr_proxies ; pid++ ) {
        
            /* Get a handle on the proxy. */
            struct proxy *p = &e->proxies[pid];
            
            /* Loop through the proxy's incomming cells and add the
               recv tasks. */
            for ( k = 0 ; k < p->nr_cells_in ; k++ )
                engine_addtasks_recv( e , p->cells_in[k] , NULL , NULL );
            
            /* Loop through the proxy's outgoing cells and add the
               send tasks. */
            for ( k = 0 ; k < p->nr_cells_out ; k++ )
                engine_addtasks_send( e , p->cells_out[k] , p->cells_in[0] );
            
            }
        
    #endif
        
    /* Rank the tasks. */
    scheduler_ranktasks( sched );
    
    /* Weight the tasks. */
    scheduler_reweight( sched );
            
    /* Set the tasks age. */
    e->tasks_age = 0;
            
    }
    
    

/**
 * @brief Mark tasks to be skipped and set the sort flags accordingly.
 * 
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */
 
int engine_marktasks ( struct engine *e ) {

    struct scheduler *s = &e->sched;
    int k, nr_tasks = s->nr_tasks, *ind = s->tasks_ind;
    struct task *t, *tasks = s->tasks;
    float dt_step = e->dt_step;
    struct cell *ci, *cj;
    // ticks tic = getticks();
    
    /* Muc less to do here if we're on a fixed time-step. */
    if ( !( e->policy & engine_policy_multistep ) ) {
    
        /* Run through the tasks and mark as skip or not. */
        for ( k = 0 ; k < nr_tasks ; k++ ) {

            /* Get a handle on the kth task. */
            t = &tasks[ ind[k] ];

            /* Pair? */
            if ( t->type == task_type_pair || ( t->type == task_type_sub && t->cj != NULL ) ) {

                /* Local pointers. */
                ci = t->ci;
                cj = t->cj;

                /* Too much particle movement? */
                if ( t->tight &&
                     ( fmaxf( ci->h_max , cj->h_max ) + ci->dx_max + cj->dx_max > cj->dmin || 
                       ci->dx_max > space_maxreldx*ci->h_max || cj->dx_max > space_maxreldx*cj->h_max ) )
                    return 1;

                }
                
            /* Sort? */
            else if ( t->type == task_type_sort ) {
            
                /* If all the sorts have been done, make this task implicit. */
                if ( !( t->flags & (t->flags ^ t->ci->sorted ) ) )
                    t->implicit = 1;
            
                }

            }
            
        }
    
    else {
    
        /* Run through the tasks and mark as skip or not. */
        for ( k = 0 ; k < nr_tasks ; k++ ) {

            /* Get a handle on the kth task. */
            t = &tasks[ ind[k] ];

            /* Sort-task? Note that due to the task ranking, the sorts
               will all come before the pairs. */
            if ( t->type == task_type_sort ) {

                /* Re-set the flags. */
                t->flags = 0;
                t->skip = 1;

                }

            /* Single-cell task? */
            else if ( t->type == task_type_self ||
                      t->type == task_type_ghost ||
                    ( t->type == task_type_sub && t->cj == NULL ) ) {

                /* Set this task's skip. */
                t->skip = ( t->ci->dt_min > dt_step );

                }

            /* Pair? */
            else if ( t->type == task_type_pair || ( t->type == task_type_sub && t->cj != NULL ) ) {

                /* Local pointers. */
                ci = t->ci;
                cj = t->cj;

                /* Set this task's skip. */
                t->skip = ( ci->dt_min > dt_step && cj->dt_min > dt_step );

                /* Too much particle movement? */
                if ( t->tight &&
                     ( fmaxf( ci->h_max , cj->h_max ) + ci->dx_max + cj->dx_max > cj->dmin || 
                       ci->dx_max > space_maxreldx*ci->h_max || cj->dx_max > space_maxreldx*cj->h_max ) )
                    return 1;

                /* Set the sort flags. */
                if ( !t->skip && t->type == task_type_pair ) {
                    if ( !( ci->sorted & ( 1 << t->flags ) ) ) {
                        ci->sorts->flags |= (1 << t->flags);
                        ci->sorts->skip = 0;
                        }
                    if ( !( cj->sorted & ( 1 << t->flags ) ) ) {
                        cj->sorts->flags |= (1 << t->flags);
                        cj->sorts->skip = 0;
                        }
                    }

                }

            /* Kick2? */
            else if ( t->type == task_type_kick2 )
                t->skip = 0;

            /* None? */
            else if ( t->type == task_type_none )
                t->skip = 1;

            }
            
        }
        
    // message( "took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );
    
    /* All is well... */
    return 0;
    
    }
    

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 */
 
void engine_rebuild ( struct engine *e ) {

    int k;
    struct scheduler *sched = &e->sched;
    
    /* Clear the forcerebuild flag, whatever it was. */
    e->forcerebuild = 0;

    /* Re-build the space. */
    // tic = getticks();
    space_rebuild( e->s , 0.0 );
    // message( "space_rebuild took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );

    /* If in parallel, exchange the cell structure. */
    #ifdef WITH_MPI
        // tic = getticks();
        engine_exchange_cells( e );
        // message( "engine_exchange_cells took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );
    #endif

    /* Re-build the tasks. */
    // tic = getticks();
    engine_maketasks( e );
    // message( "engine_maketasks took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );

    /* Run through the tasks and mark as skip or not. */
    // tic = getticks();
    if ( engine_marktasks( e ) )
        error( "engine_marktasks failed after space_rebuild." );
    // message( "engine_marktasks took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );

    /* Count and print the number of each task type. */
    int counts[ task_type_count+1 ];
    for ( k = 0 ; k <= task_type_count ; k++ )
        counts[k] = 0;
    for ( k = 0 ; k < sched->nr_tasks ; k++ )
        if ( !sched->tasks[k].skip )
            counts[ (int)sched->tasks[k].type ] += 1;
        else
            counts[ task_type_count ] += 1;
    #ifdef WITH_MPI
        printf( "[%03i] engine_prepare: task counts are [ %s=%i" , e->nodeID , taskID_names[0] , counts[0] );
    #else
        printf( "engine_prepare: task counts are [ %s=%i" , taskID_names[0] , counts[0] );
    #endif
    for ( k = 1 ; k < task_type_count ; k++ )
        printf( " %s=%i" , taskID_names[k] , counts[k] );
    printf( " skipped=%i ]\n" , counts[ task_type_count ] ); fflush(stdout);
    message( "nr_parts = %i." , e->s->nr_parts );
    
    }


/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 */
 
void engine_prepare ( struct engine *e ) {
    
    int rebuild;
    struct scheduler *sched = &e->sched;
    
    TIMER_TIC

    /* Run through the tasks and mark as skip or not. */
    // tic = getticks();
    rebuild = ( e->forcerebuild || engine_marktasks( e ) );
    // message( "space_marktasks took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );
        
    /* Collect the values of rebuild from all nodes. */
    #ifdef WITH_MPI
        int buff;
        if ( MPI_Allreduce( &rebuild , &buff , 1 , MPI_INT , MPI_MAX , MPI_COMM_WORLD ) != MPI_SUCCESS )
            error( "Failed to aggreggate the rebuild flag accross nodes." );
        rebuild = buff;
    #endif
    
    /* Did this not go through? */
    if ( rebuild )
        engine_rebuild( e );
        
    /* Re-rank the tasks every now and then. */
    if ( e->tasks_age % engine_tasksreweight == 1 )
        scheduler_reweight( &e->sched );
    e->tasks_age += 1;

    /* Start the scheduler. */
    // ticks tic2 = getticks();
    scheduler_start( sched , (1 << task_type_sort) | 
                             (1 << task_type_self) |
                             (1 << task_type_pair) | 
                             (1 << task_type_sub) |
                             (1 << task_type_ghost) | 
                             (1 << task_type_kick2) |
                             (1 << task_type_send_xv) |
                             (1 << task_type_recv_xv) |
                             (1 << task_type_send_rho) |
                             (1 << task_type_recv_rho) );
    // message( "scheduler_start took %.3f ms." , (double)(getticks() - tic2)/CPU_TPS*1000 );
    
    TIMER_TOC( timer_prepare );
    
    }


/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 */
 
void engine_barrier ( struct engine *e , int tid ) {

    /* First, get the barrier mutex. */
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to get barrier mutex." );
        
    /* This thread is no longer running. */
    e->barrier_running -= 1;
        
    /* If all threads are in, send a signal... */
    if ( e->barrier_running == 0 )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast barrier full condition." );
        
    /* Wait for the barrier to open. */
    while ( e->barrier_launch == 0 || tid >= e->barrier_launchcount )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Eror waiting for barrier to close." );
        
    /* This thread has been launched. */
    e->barrier_running += 1;
    e->barrier_launch -= 1;
    
    /* If I'm the last one out, signal the condition again. */
    if ( e->barrier_launch == 0 )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast empty barrier condition." );
            
    /* Last but not least, release the mutex. */
    if ( pthread_mutex_unlock( &e->barrier_mutex ) != 0 )
        error( "Failed to get unlock the barrier mutex." );

    }
    
    
/**
 * @brief Mapping function to collect the data from the second kick.
 */

void engine_collect_kick2 ( struct cell *c ) {

    int k, updated = 0;
    float dt_min = FLT_MAX, dt_max = 0.0f;
    double ekin = 0.0, epot = 0.0;
    float mom[3] = { 0.0f , 0.0f , 0.0f }, ang[3] = { 0.0f , 0.0f , 0.0f };
    struct cell *cp;
    
    /* If I am a super-cell, return immediately. */
    if ( c->kick2 != NULL || c->count == 0 )
        return;
        
    /* If this cell is not split, I'm in trouble. */
    if ( !c->split )
        error( "Cell has no super-cell." );
        
    /* Collect the values from the progeny. */
    for ( k = 0 ; k < 8 ; k++ )
        if ( ( cp = c->progeny[k] ) != NULL ) {
            engine_collect_kick2( cp );
            dt_min = fminf( dt_min , cp->dt_min );
            dt_max = fmaxf( dt_max , cp->dt_max );
            updated += cp->updated;
            ekin += cp->ekin;
            epot += cp->epot;
            mom[0] += cp->mom[0]; mom[1] += cp->mom[1]; mom[2] += cp->mom[2];
            ang[0] += cp->ang[0]; ang[1] += cp->ang[1]; ang[2] += cp->ang[2];
            }
    
    /* Store the collected values in the cell. */
    c->dt_min = dt_min;
    c->dt_max = dt_max;
    c->updated = updated;
    c->ekin = ekin;
    c->epot = epot;
    c->mom[0] = mom[0]; c->mom[1] = mom[1]; c->mom[2] = mom[2];
    c->ang[0] = ang[0]; c->ang[1] = ang[1]; c->ang[2] = ang[2];
        
    }


/**
 * @brief Compute the force on a single particle brute-force.
 */

// void engine_single_density ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {
// 
//     int i, k;
//     double r2, dx[3];
//     float fdx[3], ih;
//     struct part p;
//     
//     /* Find "our" part. */
//     for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
//     if ( k == N )
//         error( "Part not found." );
//     p = parts[k];
//     
//     /* Clear accumulators. */
//     ih = 1.0f / p.h;
//     p.rho = 0.0f; p.rho_dh = 0.0f;
//     p.density.wcount = 0.0f; p.density.wcount_dh = 0.0f;
// 	p.density.div_v = 0.0;
// 	for ( k=0 ; k < 3 ; k++)
// 		p.density.curl_v[k] = 0.0;
//             
//     /* Loop over all particle pairs (force). */
//     for ( k = 0 ; k < N ; k++ ) {
//         if ( parts[k].id == p.id )
//             continue;
//         for ( i = 0 ; i < 3 ; i++ ) {
//             dx[i] = p.x[i] - parts[k].x[i];
//             if ( periodic ) {
//                 if ( dx[i] < -dim[i]/2 )
//                     dx[i] += dim[i];
//                 else if ( dx[i] > dim[i]/2 )
//                     dx[i] -= dim[i];
//                 }
//             fdx[i] = dx[i];
//             }
//         r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
//         if ( r2 < p.h*p.h*kernel_gamma2 ) {
//             runner_iact_nonsym_density( r2 , fdx , p.h , parts[k].h , &p , &parts[k] );
//             }
//         }
//         
//     /* Dump the result. */
//     p.rho = ih * ih * ih * ( p.rho + p.mass*kernel_root );
//     p.rho_dh = p.rho_dh * ih * ih * ih * ih;
//     p.density.wcount = ( p.density.wcount + kernel_root ) * ( 4.0f / 3.0 * M_PI * kernel_gamma3 );
//     message( "part %lli (h=%e) has wcount=%e, rho=%e, rho_dh=%e." , p.id , p.h , p.density.wcount , p.rho , p.rho_dh );
//     fflush(stdout);
//     
//     }


// void engine_single_force ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {
// 
//     int i, k;
//     double r2, dx[3];
//     float fdx[3];
//     struct part p;
//     
//     /* Find "our" part. */
//     for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
//     if ( k == N )
//         error( "Part not found." );
//     p = parts[k];
//     
//     /* Clear accumulators. */
//     p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
//     p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
//             
//     /* Loop over all particle pairs (force). */
//     for ( k = 0 ; k < N ; k++ ) {
//     // for ( k = N-1 ; k >= 0 ; k-- ) {
//         if ( parts[k].id == p.id )
//             continue;
//         for ( i = 0 ; i < 3 ; i++ ) {
//             dx[i] = p.x[i] - parts[k].x[i];
//             if ( periodic ) {
//                 if ( dx[i] < -dim[i]/2 )
//                     dx[i] += dim[i];
//                 else if ( dx[i] > dim[i]/2 )
//                     dx[i] -= dim[i];
//                 }
//             fdx[i] = dx[i];
//             }
//         r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
//         if ( r2 < p.h*p.h*kernel_gamma2 || r2 < parts[k].h*parts[k].h*kernel_gamma2 ) {
//             p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
//             p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
//             runner_iact_nonsym_force( r2 , fdx , p.h , parts[k].h , &p , &parts[k] );
//             double dvdr = ( (p.v[0]-parts[k].v[0])*fdx[0] + (p.v[1]-parts[k].v[1])*fdx[1] + (p.v[2]-parts[k].v[2])*fdx[2] ) / sqrt(r2);
//             message( "part %lli and %lli interact (r=%.3e,dvdr=%.3e) with a=[%.3e,%.3e,%.3e], dudt=%.3e." ,
//                 p.id , parts[k].id , sqrt(r2) , dvdr , p.a[0] , p.a[1], p.a[2] , p.force.u_dt );
//             }
//         }
//         
//     /* Dump the result. */
//     // message( "part %lli (h=%e) has a=[%.3e,%.3e,%.3e], udt=%e." , p.id , p.h , p.a[0] , p.a[1] , p.a[2] , p.force.u_dt );
//     fflush(stdout);
//     
//     }
    
    
/**
 * @breif Launch the runners.
 *
 * @param e The #engine.
 * @param nr_runners The number of #runners to let loose.
 */
 
void engine_launch ( struct engine *e , int nr_runners ) {

    /* Cry havoc and let loose the dogs of war. */
    e->barrier_launch = nr_runners;
    e->barrier_launchcount = nr_runners;
    if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
        error( "Failed to broadcast barrier open condition." );
        
    /* Sit back and wait for the runners to come home. */
    while ( e->barrier_launch || e->barrier_running )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for barrier." );
            
    }
    
    
void hassorted ( struct cell *c ) {

    if ( c->sorted )
        error( "Suprious sorted flags." );
        
    if ( c->split )
        for ( int k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL )
                hassorted( c->progeny[k] );
                    
    }


/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 * @param sort_queues Flag to try to sort the queues topologically.
 */
 
void engine_step ( struct engine *e ) {

    int k;
    float dt = e->dt, dt_step, dt_max = 0.0f, dt_min = FLT_MAX;
    double epot = 0.0, ekin = 0.0;
    float mom[3] = { 0.0 , 0.0 , 0.0 };
    float ang[3] = { 0.0 , 0.0 , 0.0 };
    int count = 0;
    struct cell *c;
    struct space *s = e->s;
    
    TIMER_TIC2

    /* Get the maximum dt. */
    if ( e->policy & engine_policy_multistep ) {
        dt_step = 2.0f*dt;
        for ( k = 0 ; k < 32 && (e->step & (1 << k)) == 0 ; k++ )
            dt_step *= 2;
        }
    else
        dt_step = FLT_MAX;
        
    /* Set the maximum dt. */
    e->dt_step = dt_step;
    e->s->dt_step = dt_step;
    // message( "dt_step set to %.3e (dt=%.3e)." , dt_step , e->dt ); fflush(stdout);
    
    // printParticle( parts , 432626 );
    
    /* First kick. */
    if ( e->step == 0 || !( e->policy & engine_policy_fixdt ) ) {
        TIMER_TIC
        scheduler_start( &e->sched , (1 << task_type_kick1) );
        engine_launch( e , ( e->nr_threads > 8 ) ? 8 : e->nr_threads );
        TIMER_TOC( timer_kick1 );
        }
    
    /* Check if all the kick1 threads have executed. */
    /* for ( k = 0 ; k < e->sched.nr_tasks ; k++ )
        if ( e->sched.tasks[k].type == task_type_kick1 &&
             e->sched.tasks[k].toc == 0 )
            error( "Not all kick1 tasks completed." ); */
        
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( e->s->parts , 3392063069037 , e->s->nr_parts );
 
    /* Re-distribute the particles amongst the nodes? */
    if ( e->forcerepart )
        engine_repartition( e );
    
    /* Prepare the space. */
    engine_prepare( e );
    
    // engine_single_density( e->s->dim , 3392063069037 , e->s->parts , e->s->nr_parts , e->s->periodic );

    /* Send off the runners. */
    TIMER_TIC
    engine_launch( e , e->nr_threads );
    TIMER_TOC(timer_runners);
    
    // engine_single_force( e->s->dim , 8328423931905 , e->s->parts , e->s->nr_parts , e->s->periodic );
    
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( parts , 432626 );
    // printParticle( e->s->parts , 3392063069037 , e->s->nr_parts );
    // printParticle( e->s->parts , 8328423931905 , e->s->nr_parts );

    /* Collect the cell data from the second kick. */
    for ( k = 0 ; k < s->nr_cells ; k++ )
        if ( s->cells[k].nodeID == e->nodeID ) {
            c = &s->cells[k];
            engine_collect_kick2( c );
            dt_min = fminf( dt_min , c->dt_min );
            dt_max = fmaxf( dt_max , c->dt_max );
            ekin += c->ekin;
            epot += c->epot;
            count += c->updated;
            mom[0] += c->mom[0]; mom[1] += c->mom[1]; mom[2] += c->mom[2];
            ang[0] += c->ang[0]; ang[1] += c->ang[1]; ang[2] += c->ang[2];
            }
        
    /* Aggregate the data from the different nodes. */
    #ifdef WITH_MPI
        double in[3], out[3];
        out[0] = dt_min;
        if ( MPI_Allreduce( out , in , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD ) != MPI_SUCCESS )
            error( "Failed to aggregate dt_min." );
        dt_min = in[0];
        out[0] = dt_max;
        if ( MPI_Allreduce( out , in , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD ) != MPI_SUCCESS )
            error( "Failed to aggregate dt_max." );
        dt_max = in[0];
        out[0] = count; out[1] = ekin; out[2] = epot;
        if ( MPI_Allreduce( out , in , 3 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ) != MPI_SUCCESS )
            error( "Failed to aggregate energies." );
        count = in[0]; ekin = in[1]; epot = in[2];
        /* int nr_parts;
        if ( MPI_Allreduce( &s->nr_parts , &nr_parts , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ) != MPI_SUCCESS )
            error( "Failed to aggregate particle count." );
        if ( e->nodeID == 0 )
            message( "nr_parts=%i." , nr_parts ); */
    #endif
    
    e->dt_min = dt_min;
    e->dt_max = dt_max;
    e->count_step = count;
    e->ekin = ekin;
    e->epot = epot;
    // printParticle( e->s->parts , 382557 , e->s->nr_parts );
    // message( "dt_min/dt_max is %e/%e." , dt_min , dt_max ); fflush(stdout);
    // message( "etot is %e (ekin=%e, epot=%e)." , ekin+epot , ekin , epot ); fflush(stdout);
    // message( "total momentum is [ %e , %e , %e ]." , mom[0] , mom[1] , mom[2] ); fflush(stdout);
    // message( "total angular momentum is [ %e , %e , %e ]." , ang[0] , ang[1] , ang[2] ); fflush(stdout);
    // message( "updated %i parts (dt_step=%.3e)." , count , dt_step ); fflush(stdout);
        
    /* Increase the step. */
    e->step += 1;

    /* Does the time step need adjusting? */
    if ( e->policy & engine_policy_fixdt ) {
        dt = e->dt_orig;
        }
    else {
        if ( dt == 0 ) {
            e->nullstep += 1;
            if ( e->dt_orig > 0.0 ) {
                dt = e->dt_orig;
                while ( dt_min < dt )
                    dt *= 0.5;
                while ( dt_min > 2*dt )
                    dt *= 2.0;
                }
            else
                dt = dt_min;
            for ( k = 0 ; k < s->nr_parts ; k++ ) {
                /* struct part *p = &s->parts[k];
                struct xpart *xp = &s->xparts[k];
                float dt_curr = dt;
                for ( int j = (int)( p->dt / dt ) ; j > 1 ; j >>= 1 )
                    dt_curr *= 2.0f; 
                xp->dt_curr = dt_curr; */
                s->parts[k].dt = dt;
                s->xparts[k].dt_curr = dt;
                }
            // message( "dt_min=%.3e, adjusting time step to dt=%e." , dt_min , e->dt );
            }
        else {
            while ( dt_min < dt ) {
                dt *= 0.5;
                e->step *= 2;
                e->nullstep *= 2;
                // message( "dt_min dropped below time step, adjusting to dt=%e." , e->dt );
                }
            while ( dt_min > 2*dt && (e->step & 1) == 0 ) {
                dt *= 2.0;
                e->step /= 2;
                e->nullstep /= 2;
                // message( "dt_min is larger than twice the time step, adjusting to dt=%e." , e->dt );
                }
            }
        } 
    e->dt = dt;
    
    /* Set the system time. */
    e->time = dt * (e->step - e->nullstep);
        
    TIMER_TOC2(timer_step);
    
    }
    
    
/**
 * @brief Create and fill the proxies.
 *
 * @param e The #engine.
 */
 
void engine_makeproxies ( struct engine *e ) {

    int i, j, k, ii, jj, kk;
    int cid, cjd, pid, ind[3], *cdim = e->s->cdim;
    struct space *s = e->s;
    struct cell *cells = s->cells;
    struct proxy *proxies = e->proxies;
    
    /* Prepare the proxies and the proxy index. */
    if ( e->proxy_ind == NULL )
        if ( ( e->proxy_ind = (int *)malloc( sizeof(int) * e->nr_nodes ) ) == NULL )
            error( "Failed to allocate proxy index." );
    for ( k = 0 ; k < e->nr_nodes ; k++ )
        e->proxy_ind[k] = -1;
    e->nr_proxies = 0;
    
    /* The following loop is super-clunky, but it's necessary
       to ensure that the order of the send and recv cells in
       the proxies is identical for all nodes! */
    
    /* Loop over each cell in the space. */
    for ( ind[0] = 0 ; ind[0] < cdim[0] ; ind[0]++ )
        for ( ind[1] = 0 ; ind[1] < cdim[1] ; ind[1]++ )
            for ( ind[2] = 0 ; ind[2] < cdim[2] ; ind[2]++ ) {
            
                /* Get the cell ID. */
                cid = cell_getid( cdim , ind[0] , ind[1] , ind[2] );
                
                /* Loop over all its neighbours (periodic). */
                for ( i = -1 ; i <= 1 ; i++ ) {
                    ii = ind[0] + i;
                    if ( ii >= cdim[0] )
                        ii -= cdim[0];
                    else if ( ii < 0 )
                        ii += cdim[0];
                    for ( j = -1 ; j <= 1 ; j++ ) {
                        jj = ind[1] + j;
                        if ( jj >= cdim[1] )
                            jj -= cdim[1];
                        else if ( jj < 0 )
                            jj += cdim[1];
                        for ( k = -1 ; k <= 1 ; k++ ) {
                            kk = ind[2] + k;
                            if ( kk >= cdim[2] )
                                kk -= cdim[2];
                            else if ( kk < 0 )
                                kk += cdim[2];
                            
                            /* Get the cell ID. */
                            cjd = cell_getid( cdim , ii , jj , kk );
                            
                            /* Add to proxies? */
                            if ( cells[cid].nodeID == e->nodeID && cells[cjd].nodeID != e->nodeID ) {
                                pid = e->proxy_ind[ cells[cjd].nodeID ];
                                if ( pid < 0 ) {
                                    if ( e->nr_proxies == engine_maxproxies )
                                        error( "Maximum number of proxies exceeded." );
                                    proxy_init( &proxies[ e->nr_proxies ] , e->nodeID , cells[cjd].nodeID );
                                    e->proxy_ind[ cells[cjd].nodeID ] = e->nr_proxies;
                                    pid = e->nr_proxies;
                                    e->nr_proxies += 1;
                                    }
                                proxy_addcell_in( &proxies[pid] , &cells[cjd] );
                                proxy_addcell_out( &proxies[pid] , &cells[cid] );
                                cells[cid].sendto |= ( 1ULL << pid );
                                }
                                
                            if ( cells[cjd].nodeID == e->nodeID && cells[cid].nodeID != e->nodeID ) {
                                pid = e->proxy_ind[ cells[cid].nodeID ];
                                if ( pid < 0 ) {
                                    if ( e->nr_proxies == engine_maxproxies )
                                        error( "Maximum number of proxies exceeded." );
                                    proxy_init( &proxies[ e->nr_proxies ] , e->nodeID , cells[cid].nodeID );
                                    e->proxy_ind[ cells[cid].nodeID ] = e->nr_proxies;
                                    pid = e->nr_proxies;
                                    e->nr_proxies += 1;
                                    }
                                proxy_addcell_in( &proxies[pid] , &cells[cid] );
                                proxy_addcell_out( &proxies[pid] , &cells[cjd] );
                                cells[cjd].sendto |= ( 1ULL << pid );
                                }
                            }
                        }
                    }
                }
        
    }
    
    
/** 
 * @brief Split the underlying space according to the given grid.
 *
 * @param e The #engine.
 * @param grid The grid.
 */
 
void engine_split ( struct engine *e , int *grid ) {

    int j, k;
    float scale[3];
    int ind[3];
    struct space *s = e->s;
    struct cell *c;
    struct part *p;
    
    /* If we've got the wrong number of nodes, fail. */
    if ( e->nr_nodes != grid[0]*grid[1]*grid[2] )
        error( "Grid size does not match number of nodes." );
        
    /* Get the scale. */
    for ( j = 0 ; j < 3 ; j++ )
        scale[j] = ((float)grid[j]) / s->cdim[j];
    
    /* Run through the cells and set their nodeID. */
    for ( k = 0 ; k < s->nr_cells ; k++ ) {
        c = &s->cells[k];
        for ( j = 0 ; j < 3 ; j++ )
            ind[j] = c->loc[j] * s->ih[j] * scale[j];
        c->nodeID = ind[0] + grid[0]*( ind[1] + grid[1]*ind[2] );
        }
        
    /* Make the proxies. */
    engine_makeproxies( e );
        
    /* For now, just kill any particle outside of our grid. */
    for ( k = 0 ; k < s->nr_parts ; k++ ) {
        p = &s->parts[k];
        if ( s->cells[ cell_getid( s->cdim , p->x[0]*s->ih[0] , p->x[1]*s->ih[1] , p->x[2]*s->ih[2] ) ].nodeID != e->nodeID ) {
            s->nr_parts -= 1;
            s->parts[k] = s->parts[ s->nr_parts ];
            s->xparts[k] = s->xparts[ s->nr_parts ];
            k -= 1;
            }
        }
    
    /* Re-allocate the local parts. */
    s->size_parts = s->nr_parts * 2;
    struct part *parts_new;
    struct xpart *xparts_new;
    if ( posix_memalign( (void **)&parts_new , part_align , sizeof(struct part) * s->size_parts ) != 0 ||
         posix_memalign( (void **)&xparts_new , part_align , sizeof(struct xpart) * s->size_parts ) != 0 )
        error( "Failed to allocate new part data." );
    memcpy( parts_new , s->parts , sizeof(struct part) * s->nr_parts );
    memcpy( xparts_new , s->xparts , sizeof(struct xpart) * s->nr_parts );
    free( s->parts );
    free( s->xparts );
    s->parts = parts_new;
    s->xparts = xparts_new;

    }
    
    
/**
 * @brief init an engine with the given number of threads, queues, and
 *      the given policy.
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param dt The initial time step to use.
 * @param nr_threads The number of threads to spawn.
 * @param nr_queues The number of task queues to create.
 * @param policy The queueing policy to use.
 */
 
void engine_init ( struct engine *e , struct space *s , float dt , int nr_threads , int nr_queues , int nr_nodes , int nodeID , int policy ) {

    int k;
    float dt_min = dt;
    #if defined(HAVE_SETAFFINITY)
        int nr_cores = sysconf( _SC_NPROCESSORS_ONLN );
        int i, j, cpuid[ nr_cores ];
        cpu_set_t cpuset;
        if ( policy & engine_policy_cputight ) {
            for ( k = 0 ; k < nr_cores ; k++ )
                cpuid[k] = k;
            }
        else {
            cpuid[0] = 0;
            k = 1;
            for ( i = 1 ; i < nr_cores ; i *= 2 )
                for ( j = nr_cores / i / 2 ; j < nr_cores ; j += nr_cores / i )
                    cpuid[k++] = j;
            #ifdef WITHMPI
                printf( "engine_init: cpu map is [ " );
            #else
                printf( "[%03i] engine_init: cpu map is [ " , nodeID );
            #endif
            for ( i = 0 ; i < nr_cores ; i++ )
                printf( "%i " , cpuid[i] );
            printf( "].\n" );
            }
    #endif
    
    /* Store the values. */
    e->s = s;
    e->nr_threads = nr_threads;
    e->policy = policy;
    e->step = 0;
    e->nullstep = 0;
    e->time = 0.0;
    e->nr_nodes = nr_nodes;
    e->nodeID = nodeID;
    e->proxy_ind = NULL;
    e->nr_proxies = 0;
    e->forcerebuild = 1;
    e->forcerepart = 0;
    engine_rank = nodeID;
    
    /* Make the space link back to the engine. */
    s->e = e;
    
    /* Are we doing stuff in parallel? */
    if ( nr_nodes > 1 ) {
        #ifndef HAVE_MPI
            error( "SWIFT was not compiled with MPI support." );
        #else
            e->policy |= engine_policy_mpi;
            if ( ( e->proxies = (struct proxy *)malloc( sizeof(struct proxy) * engine_maxproxies ) ) == NULL )
                error( "Failed to allocate memory for proxies." );
            bzero( e->proxies , sizeof(struct proxy) * engine_maxproxies );
            e->nr_proxies = 0;
        #endif
        }
    
    /* First of all, init the barrier and lock it. */
    if ( pthread_mutex_init( &e->barrier_mutex , NULL ) != 0 )
        error( "Failed to initialize barrier mutex." );
    if ( pthread_cond_init( &e->barrier_cond , NULL ) != 0 )
        error( "Failed to initialize barrier condition variable." );
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to lock barrier mutex." );
    e->barrier_running = 0;
    e->barrier_launch = 0;
    e->barrier_launchcount = 0;
    
    /* Run through the parts and get the minimum time step. */
    e->dt_orig = dt;
    for ( k = 0 ; k < s->nr_parts ; k++ )
        if ( s->parts[k].dt < dt_min )
            dt_min = s->parts[k].dt;
    if ( dt_min == 0.0f )
        dt = 0.0f;
    else
        while ( dt > dt_min )
            dt *= 0.5f;
    e->dt = dt;
    
    /* Init the scheduler. */
    scheduler_init( &e->sched , e->s , nr_queues , scheduler_flag_steal , e->nodeID );
    s->nr_queues = nr_queues;
        
    /* Append a kick1 task to each cell. */
    scheduler_reset( &e->sched , s->tot_cells );
    space_map_cells_pre( e->s , 1 , &scheduler_map_mkkick1 , &e->sched );
    scheduler_ranktasks( &e->sched );
    
    /* Allocate and init the threads. */
    if ( ( e->runners = (struct runner *)malloc( sizeof(struct runner) * nr_threads ) ) == NULL )
        error( "Failed to allocate threads array." );
    for ( k = 0 ; k < nr_threads ; k++ ) {
        e->runners[k].id = k;
        e->runners[k].e = e;
        e->barrier_running += 1;
        if ( pthread_create( &e->runners[k].thread , NULL , &runner_main , &e->runners[k] ) != 0 )
            error( "Failed to create runner thread." );
        if ( e->policy & engine_policy_setaffinity ) {
            #if defined(HAVE_SETAFFINITY)

                /* Set a reasonable queue ID. */
                e->runners[k].cpuid = cpuid[ k % nr_cores ];
                if ( nr_queues < nr_threads )
                    e->runners[k].qid = cpuid[ k % nr_cores ] * nr_queues / nr_cores;
                else
                    e->runners[k].qid = k;

                /* Set the cpu mask to zero | e->id. */
                CPU_ZERO( &cpuset );
                CPU_SET( cpuid[ k % nr_cores ] , &cpuset );

                /* Apply this mask to the runner's pthread. */
                if ( pthread_setaffinity_np( e->runners[k].thread , sizeof(cpu_set_t) , &cpuset ) != 0 )
                    error( "Failed to set thread affinity." );

            #else
                error( "SWIFT was not compiled with affinity enabled." );
            #endif
            }
        else {
            e->runners[k].cpuid = k;
            e->runners[k].qid = k * nr_queues / nr_threads;
            }
        // message( "runner %i on cpuid=%i with qid=%i." , e->runners[k].id , e->runners[k].cpuid , e->runners[k].qid );
        }
        
    /* Wait for the runner threads to be in place. */
    while ( e->barrier_running || e->barrier_launch )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for runner threads to get in place." );
    
    }
    
    
    
