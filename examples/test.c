/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include <fenv.h>
#include <omp.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
    #include <zlib.h>
#endif

/* MPI headers. */
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* Local headers. */
#include "swift.h"

/* Ticks per second on this machine. */
#ifndef CPU_TPS
    #define CPU_TPS 2.40e9
#endif

/* Engine policy flags. */
#ifndef ENGINE_POLICY
    #define ENGINE_POLICY engine_policy_none
#endif


/**
 * @brief Mapping function to draw a specific cell (gnuplot).
 */

void map_cells_plot ( struct cell *c , void *data ) {

    int depth = *(int *)data;
    double *l = c->loc, *h = c->h;

    if ( c->depth <= depth ) {
    
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0] , l[1]+h[1] , l[2] );
    
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0] , l[1]+h[1] , l[2]+h[2] );
    
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1]+h[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1]+h[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0] , l[1] , l[2]+h[2] );
    
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0]+h[0] , l[1] , l[2] +h[2]);
    
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0]+h[0] , l[1] , l[2] );
    
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1]+h[1] , l[2] );
        printf( "%.16e %.16e %.16e\n" , l[0] , l[1]+h[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2]+h[2] );
        printf( "%.16e %.16e %.16e\n\n\n" , l[0]+h[0] , l[1]+h[1] , l[2] );
        
        if ( !c->split ) {
            for ( int k = 0 ; k < c->count ; k++ )
                printf( "0 0 0 %.16e %.16e %.16e\n" ,
                    c->parts[k].x[0] , c->parts[k].x[1] , c->parts[k].x[2] );
            printf( "\n\n" );
            }
        /* else
            for ( int k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    map_cells_plot( c->progeny[k] , data ); */
    
        }

    }


/**
 * @brief Mapping function for checking if each part is in its box.
 */

/* void map_check ( struct part *p , struct cell *c , void *data ) {

    if ( p->x[0] < c->loc[0] || p->x[0] > c->loc[0]+c->h[0] ||
         p->x[0] < c->loc[0] || p->x[0] > c->loc[0]+c->h[0] ||
         p->x[0] < c->loc[0] || p->x[0] > c->loc[0]+c->h[0] )
        printf( "map_check: particle %i is outside of its box.\n" , p->id );

    } */


/**
 * @brief Mapping function for neighbour count.
 */

void map_cellcheck ( struct cell *c , void *data ) {

    int k, *count = (int *)data;
    struct part *p;
    
    __sync_fetch_and_add( count , c->count );
    
    /* Loop over all parts and check if they are in the cell. */
    for ( k = 0 ; k < c->count ; k++ ) {
        p = &c->parts[k];
        if ( p->x[0] < c->loc[0] || p->x[1] < c->loc[1] || p->x[2] < c->loc[2] ||
             p->x[0] > c->loc[0] + c->h[0] || p->x[1] > c->loc[1] + c->h[1] || p->x[2] > c->loc[2] + c->h[2] ) {
            printf( "map_cellcheck: particle at [ %.16e %.16e %.16e ] outside of cell [ %.16e %.16e %.16e ] - [ %.16e %.16e %.16e ].\n" ,
                p->x[0] , p->x[1] , p->x[2] ,
                c->loc[0] , c->loc[1] , c->loc[2] ,
                c->loc[0] + c->h[0] , c->loc[1] + c->h[1] , c->loc[2] + c->h[2] );
            error( "particle out of bounds!" );
            }
        }

    }


/**
 * @brief Mapping function for maxdepth cell count.
 */

void map_maxdepth ( struct cell *c , void *data ) {

    int maxdepth = ((int *)data)[0];
    int *count = &((int *)data)[1];
    
    // printf( "%e\n" , p->count );

    if ( c->depth == maxdepth )
        *count += 1;

    }


/**
 * @brief Mapping function for neighbour count.
 */

void map_count ( struct part *p , struct cell *c , void *data ) {

    double *wcount = (double *)data;
    
    // printf( "%i %e %e\n" , p->id , p->count , p->count_dh );

    *wcount += p->density.wcount;

    }

void map_wcount_min ( struct part *p , struct cell *c , void *data ) {

    struct part **p2 = (struct part **)data;
    
    if ( p->density.wcount < (*p2)->density.wcount )
        *p2 = p;

    }

void map_wcount_max ( struct part *p , struct cell *c , void *data ) {

    struct part **p2 = (struct part **)data;
    
    if ( p->density.wcount > (*p2)->density.wcount )
        *p2 = p;

    }

void map_h_min ( struct part *p , struct cell *c , void *data ) {

    struct part **p2 = (struct part **)data;
    
    if ( p->h < (*p2)->h )
        *p2 = p;

    }

void map_h_max ( struct part *p , struct cell *c , void *data ) {

    struct part **p2 = (struct part **)data;
    
    if ( p->h > (*p2)->h )
        *p2 = p;

    }


/**
 * @brief Mapping function for neighbour count.
 */

void map_icount ( struct part *p , struct cell *c , void *data ) {

    // int *count = (int *)data;
    
    // printf( "%i\n" , p->icount );

    // *count += p->icount;

    }


/**
 * @brief Mapping function to print the particle position.
 */

void map_dump ( struct part *p , struct cell *c , void *data ) {

    double *shift = (double *)data;

    printf( "%g\t%g\t%g\n" , p->x[0]-shift[0] , p->x[1]-shift[1] , p->x[2]-shift[2] );

    }


/**
 * @brief Compute the average number of pairs per particle using
 *      a brute-force O(N^2) computation.
 *
 * @param dim The space dimensions.
 * @param parts The #part array.
 * @param N The number of parts.
 * @param periodic Periodic boundary conditions flag.
 */

void pairs_n2 ( double *dim , struct part *__restrict__ parts , int N , int periodic ) {

    int i, j, k, count = 0;
    // int mj, mk;
    // double maxratio = 1.0;
    double r2, dx[3], rho = 0.0;
    double rho_max = 0.0, rho_min = 100;
    
    /* Loop over all particle pairs. */
    #pragma omp parallel for schedule(dynamic), default(none), private(k,i,dx,r2), shared(periodic,parts,dim,N,stdout)
    for ( j = 0 ; j < N ; j++ ) {
        if ( j % 1000 == 0 ) {
            printf( "pairs_n2: j=%i.\n" , j );
            fflush(stdout);
            }
        for ( k = j+1 ; k < N ; k++ ) {
            for ( i = 0 ; i < 3 ; i++ ) {
                dx[i] = parts[j].x[i] - parts[k].x[i];
                if ( periodic ) {
                    if ( dx[i] < -dim[i]/2 )
                        dx[i] += dim[i];
                    else if ( dx[i] > dim[i]/2 )
                        dx[i] -= dim[i];
                    }
                }
            r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            if ( r2 < parts[j].h*parts[j].h || r2 < parts[k].h*parts[k].h ) {
                runner_iact_density( r2 , NULL , parts[j].h , parts[k].h , &parts[j] , &parts[k] );
                /* if ( parts[j].h / parts[k].h > maxratio )
                    #pragma omp critical
                    {
                    maxratio = parts[j].h / parts[k].h;
                    mj = j; mk = k;
                    }
                else if ( parts[k].h / parts[j].h > maxratio )
                    #pragma omp critical
                    {
                    maxratio = parts[k].h / parts[j].h;
                    mj = j; mk = k;
                    } */
                }
            }
        }
        
    /* Aggregate the results. */
    for ( k = 0 ; k < N ; k++ ) {
        // count += parts[k].icount;
        rho += parts[k].density.wcount;
        rho_min = fmin( parts[k].density.wcount , rho_min );
        rho_min = fmax( parts[k].density.wcount , rho_max );
        }
            
    /* Dump the result. */
    printf( "pairs_n2: avg. density per part is %.3f (nr. pairs %.3f).\n" , rho/N + 32.0/3 , ((double)count)/N );
    printf( "pairs_n2: densities are in [ %e , %e ].\n" , rho_min/N + 32.0/3 , rho_max/N + 32.0/3 );
    /* printf( "pairs_n2: maximum ratio between parts %i [%e,%e,%e] and %i [%e,%e,%e] is %.3f/%.3f\n" ,
        mj , parts[mj].x[0] , parts[mj].x[1] , parts[mj].x[2] ,
        mk , parts[mk].x[0] , parts[mk].x[1] , parts[mk].x[2] ,
        parts[mj].h , parts[mk].h ); fflush(stdout); */
    fflush(stdout);
            
    }


void pairs_single_density ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {

    int i, k;
    // int mj, mk;
    // double maxratio = 1.0;
    double r2, dx[3];
    float fdx[3];
    struct part p;
    // double ih = 12.0/6.25;
    
    /* Find "our" part. */
    for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
    if ( k == N )
        error( "Part not found." );
    p = parts[k];
    printf( "pairs_single: part[%i].id == %lli.\n" , k , pid );
    
    p.rho = 0.0;
    p.density.wcount = 0.0;
    // p.icount = 0;
    p.rho_dh = 0.0;
            
    /* Loop over all particle pairs. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( parts[k].id == p.id )
            continue;
        for ( i = 0 ; i < 3 ; i++ ) {
            dx[i] = p.x[i] - parts[k].x[i];
            if ( periodic ) {
                if ( dx[i] < -dim[i]/2 )
                    dx[i] += dim[i];
                else if ( dx[i] > dim[i]/2 )
                    dx[i] -= dim[i];
                }
            fdx[i] = dx[i];
            }
        r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
        if ( r2 < p.h*p.h ) {
            runner_iact_nonsym_density( r2 , fdx , p.h , parts[k].h , &p , &parts[k] );
            /* printf( "pairs_simple: interacting particles %lli [%i,%i,%i] and %lli [%i,%i,%i], r=%e.\n" ,
                pid , (int)(p.x[0]*ih) , (int)(p.x[1]*ih) , (int)(p.x[2]*ih) ,
                parts[k].id , (int)(parts[k].x[0]*ih) , (int)(parts[k].x[1]*ih) , (int)(parts[k].x[2]*ih) ,
                sqrtf(r2) ); */
            }
        }
        
    /* Dump the result. */
    printf( "pairs_single: wcount of part %lli (h=%e) is %f.\n" , p.id , p.h , p.density.wcount + 32.0/3 );
    fflush(stdout);
    
    }


void pairs_single_grav ( double *dim , long long int pid , struct gpart *__restrict__ parts , int N , int periodic ) {

    int i, k;
    // int mj, mk;
    // double maxratio = 1.0;
    double r2, dx[3];
    float fdx[3], a[3] = { 0.0 , 0.0 , 0.0 }, aabs[3] = { 0.0 , 0.0 , 0.0 };
    struct gpart pi, pj;
    // double ih = 12.0/6.25;
    
    /* Find "our" part. */
    for ( k = 0 ; k < N ; k++ )
        if ( ( parts[k].id > 0 && parts[k].part->id == pid ) || parts[k].id == -pid )
            break;
    if ( k == N )
        error( "Part not found." );
    pi = parts[k];
    pi.a[0] = 0.0f; pi.a[1] = 0.0f; pi.a[2] = 0.0f;
    
    /* Loop over all particle pairs. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( parts[k].id == pi.id )
            continue;
        pj = parts[k];
        for ( i = 0 ; i < 3 ; i++ ) {
            dx[i] = pi.x[i] - pj.x[i];
            if ( periodic ) {
                if ( dx[i] < -dim[i]/2 )
                    dx[i] += dim[i];
                else if ( dx[i] > dim[i]/2 )
                    dx[i] -= dim[i];
                }
            fdx[i] = dx[i];
            }
        r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
        runner_iact_grav( r2 , fdx , &pi , &pj );
        a[0] += pi.a[0]; a[1] += pi.a[1]; a[2] += pi.a[2];
        aabs[0] += fabsf( pi.a[0] ); aabs[1] += fabsf( pi.a[1] ); aabs[2] += fabsf( pi.a[2] );
        pi.a[0] = 0.0f; pi.a[1] = 0.0f; pi.a[2] = 0.0f;
        }
        
    /* Dump the result. */
    message( "acceleration on gpart %lli is a=[ %e %e %e ], |a|=[ %.2e %.2e %.2e ].\n" , pi.part->id , a[0] , a[1] , a[2] , aabs[0] , aabs[1] , aabs[2] );
    
    }


/**
 * @brief Test the kernel function by dumping it in the interval [0,1].
 *
 * @param N number of intervals in [0,1].
 */
 
void kernel_dump ( int N ) {

    int k;
    float x, w, dw_dx;
    float x4[4] = {0.0f,0.0f,0.0f,0.0f};
    float w4[4] = {0.0f,0.0f,0.0f,0.0f};
    // float dw_dx4[4] __attribute__ ((aligned (16)));

    for ( k = 0 ; k <= N ; k++ ) {
        x = ((float)k) / N;
        x4[3] = x4[2]; x4[2] = x4[1]; x4[1] = x4[0]; x4[0] = x;
        kernel_deval( x , &w , &dw_dx );
        // kernel_deval_vec( (vector *)x4 , (vector *)w4 , (vector *)dw_dx4 );
        printf( " %e %e %e %e %e %e %e\n" , x , w , dw_dx , w4[0] , w4[1] , w4[2] , w4[3] );
        }

    }


void gravity_dump ( float r_max , int N ) {

    int k;
    float x, w;
    float x4[4] = {0.0f,0.0f,0.0f,0.0f};
    float w4[4] = {0.0f,0.0f,0.0f,0.0f};
    // float dw_dx4[4] __attribute__ ((aligned (16)));
    
    float gadget ( float r ) {
        float fac, h_inv, u, r2 = r*r;
        if ( r >= const_epsilon )
            fac = 1.0f / (r2 * r);
        else {
            h_inv = 1. / const_epsilon;
            u = r * h_inv;
            if ( u < 0.5 )
                fac = const_iepsilon3 * (10.666666666667 + u * u * (32.0 * u - 38.4));
            else
                fac = const_iepsilon3 * (21.333333333333 - 48.0 * u +
                                     38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
            }
        return const_G * fac;
        }

    for ( k = 1 ; k <= N ; k++ ) {
        x = (r_max * k) / N;
        x4[3] = x4[2]; x4[2] = x4[1]; x4[1] = x4[0]; x4[0] = x;
        kernel_grav_eval( x , &w );
        w *= const_G / ( x*x*x );
        // blender_deval_vec( (vector *)x4 , (vector *)w4 , (vector *)dw_dx4 );
        printf( " %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n" , x , w*x , w4[0] , w4[1] , w4[2] , w4[3] , gadget(x)*x );
        }

    }


/**
 * @brief Test the density function by dumping it for two random parts.
 *
 * @param N number of intervals in [0,1].
 */
 
void density_dump ( int N ) {

    int k;
    float r2[4] = {0.0f,0.0f,0.0f,0.0f}, hi[4], hj[4];
    struct part *pi[4], *pj[4], Pi[4], Pj[4];
    
    /* Init the interaction parameters. */
    for ( k = 0 ; k < 4 ; k++ ) {
        Pi[k].mass = 1.0f; Pi[k].rho = 0.0f; Pi[k].density.wcount = 0.0f;
        Pj[k].mass = 1.0f; Pj[k].rho = 0.0f; Pj[k].density.wcount = 0.0f;
        hi[k] = 1.0;
        hj[k] = 1.0;
        pi[k] = &Pi[k];
        pj[k] = &Pj[k];
        }

    for ( k = 0 ; k <= N ; k++ ) {
        r2[3] = r2[2]; r2[2] = r2[1]; r2[1] = r2[0];
        r2[0] = ((float)k) / N;
        Pi[0].density.wcount = 0; Pj[0].density.wcount = 0;
        runner_iact_density( r2[0] , NULL , hi[0] , hj[0] , &Pi[0] , &Pj[0] );
        printf( " %e %e %e" , r2[0] , Pi[0].density.wcount , Pj[0].density.wcount );
        Pi[0].density.wcount = 0; Pj[0].density.wcount = 0;
        Pi[1].density.wcount = 0; Pj[1].density.wcount = 0;
        Pi[2].density.wcount = 0; Pj[2].density.wcount = 0;
        Pi[3].density.wcount = 0; Pj[3].density.wcount = 0;
        runner_iact_vec_density( r2 , NULL , hi , hj , pi , pj );
        printf( " %e %e %e %e\n" , Pi[0].density.wcount , Pi[1].density.wcount , Pi[2].density.wcount , Pi[3].density.wcount );
        }

    }


/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
 
int main ( int argc , char *argv[] ) {

    int c, icount, j, k, N = -1, periodic = 1;
    int nr_threads = 1, nr_queues = -1, runs = INT_MAX;
    int data[2];
    double dim[3] = { 1.0 , 1.0 , 1.0 }, shift[3] = { 0.0 , 0.0 , 0.0 };
    double h_max = -1.0 , scaling = 1.0;
    double clock = DBL_MAX;
    struct part *parts = NULL;
    struct space s;
    struct engine e;
    struct UnitSystem us;
    char ICfileName[200];
    float dt_max = 0.0f;
    ticks tic;
    int nr_nodes = 1, myrank = 0, grid[3] = { 1 , 1 , 1 };
    
    /* Choke on FP-exceptions. */
    // feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
    
    /* Just dump the gravity potential and leave. */
    /* gravity_dump( 0.005 , 1000 );
    return 0; */
    
#ifdef WITH_MPI
    /* Start by initializing MPI. */
    int res, prov;
    if ( ( res = MPI_Init_thread( &argc , &argv , MPI_THREAD_MULTIPLE , &prov ) ) != MPI_SUCCESS )
        error( "Call to MPI_Init failed with error %i." , res );
    if ( prov != MPI_THREAD_MULTIPLE )
        error( "MPI does not provide the level of threading required (MPI_THREAD_MULTIPLE)." );
    if ( ( res = MPI_Comm_size( MPI_COMM_WORLD , &nr_nodes ) != MPI_SUCCESS ) )
        error( "MPI_Comm_size failed with error %i." , res );
    if ( ( res = MPI_Comm_rank( MPI_COMM_WORLD , &myrank ) ) != MPI_SUCCESS )
        error( "Call to MPI_Comm_rank failed with error %i." , res );
    if ( ( res = MPI_Comm_set_errhandler( MPI_COMM_WORLD , MPI_ERRORS_RETURN ) ) != MPI_SUCCESS )
        error( "Call to MPI_Comm_set_errhandler failed with error %i." , res );
    if ( myrank == 0 )
        message( "MPI is up and running with %i nodes." , nr_nodes );
    fflush(stdout);
#endif

    /* Greeting message */
    message( "This is %s\n", package_description() );
    
    /* Init the space. */
    bzero( &s , sizeof(struct space) );

    /* Parse the options */
    while ( ( c = getopt( argc , argv  , "a:c:d:f:g:m:q:r:s:t:w:z:" ) ) != -1 )
      switch( c )
	{
	case 'a':
	  if ( sscanf( optarg , "%lf" , &scaling ) != 1 )
	    error( "Error parsing cutoff scaling." );
	  if ( myrank == 0 )
        message( "scaling cutoff by %.3f." , scaling ); fflush(stdout);
	  break;
	case 'c':
	  if ( sscanf( optarg , "%lf" , &clock ) != 1 )
	    error( "Error parsing clock." );
	  if ( myrank == 0 )
        message( "clock set to %.3e." , clock ); fflush(stdout);
	  break;
	case 'd':
	  if ( sscanf( optarg , "%f" , &dt_max ) != 1 )
	    error( "Error parsing timestep." );
	  if ( myrank == 0 )
        message( "dt set to %e." , dt_max ); fflush(stdout);
	  break;
	case 'f':
	  if( !strcpy(ICfileName, optarg))
	    error("Error parsing IC file name.");
	  // if ( myrank == 0 )
        // message("IC to be read from file '%s'.", ICfileName);
	  break;
	case 'g':
	  if ( sscanf( optarg , "%i %i %i" , &grid[0] , &grid[1] , &grid[2] ) != 3 )
	    error( "Error parsing grid." );
	  if ( myrank == 0 )
        message( "grid set to [ %i %i %i ]." , grid[0] , grid[1] , grid[2] ); fflush(stdout);
	  break;
	case 'm':
	  if ( sscanf( optarg , "%lf" , &h_max ) != 1 )
	    error( "Error parsing h_max." );
	  if ( myrank == 0 )
        message( "maximum h set to %e." , h_max ); fflush(stdout);
	  break;
	case 'q':
	  if ( sscanf( optarg , "%d" , &nr_queues ) != 1 )
	    error( "Error parsing number of queues." );
	  break;
	case 'r':
	  if ( sscanf( optarg , "%d" , &runs ) != 1 )
	    error( "Error parsing number of runs." );
	  break;
	case 's':
	  if ( sscanf( optarg , "%lf %lf %lf" , &shift[0] , &shift[1] , &shift[2] ) != 3 )
	    error( "Error parsing shift." );
	  if ( myrank == 0 )
        message( "will shift parts by [ %.3f %.3f %.3f ]." , shift[0] , shift[1] , shift[2] );
	  break;
	case 't':
	  if ( sscanf( optarg , "%d" , &nr_threads ) != 1 )
	    error( "Error parsing number of threads." );
	  omp_set_num_threads( nr_threads );
	  break;
	case 'w':
	  if ( sscanf( optarg , "%d" , &space_subsize ) != 1 )
	    error( "Error parsing sub size." );
	  if ( myrank == 0 )
        message( "sub size set to %i." , space_subsize );
	  break;
	case 'z':
	  if ( sscanf( optarg , "%d" , &space_splitsize ) != 1 )
	    error( "Error parsing split size." );
	  if ( myrank == 0 )
        message( "split size set to %i." , space_splitsize );
	  break;
	case '?':
	  error( "Unknown option." );
	  break;

	}
    

    /* How large are the parts? */
    if ( myrank == 0 ) {
        message( "sizeof(struct part) is %li bytes." , (long int)sizeof( struct part ));
        message( "sizeof(struct gpart) is %li bytes." , (long int)sizeof( struct gpart ));
        }

    /* Initilaize unit system */
    initUnitSystem(&us);
    if ( myrank == 0 )
      {
	message( "Unit system: U_M = %e g.", us.UnitMass_in_cgs );
	message( "Unit system: U_L = %e cm.", us.UnitLength_in_cgs );
	message( "Unit system: U_t = %e s.", us.UnitTime_in_cgs );
	message( "Unit system: U_I = %e A.", us.UnitCurrent_in_cgs );
	message( "Unit system: U_T = %e K.", us.UnitTemperature_in_cgs );
	message( "Density units: %e a^%f h^%f.", conversionFactor(&us, UNIT_CONV_DENSITY), aFactor(&us, UNIT_CONV_DENSITY), hFactor(&us, UNIT_CONV_DENSITY) );
	message( "Entropy units: %e a^%f h^%f.", conversionFactor(&us, UNIT_CONV_ENTROPY), aFactor(&us, UNIT_CONV_ENTROPY), hFactor(&us, UNIT_CONV_ENTROPY) );
      }

    /* Read particles and space information from (GADGET) IC */
    tic = getticks();
#if defined( WITH_MPI ) 
#if defined( HAVE_PARALLEL_HDF5 )
    read_ic_parallel( ICfileName , dim , &parts , &N , &periodic, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL );
#else
    read_ic_serial( ICfileName , dim , &parts , &N , &periodic, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL );
#endif
#else
    read_ic_single( ICfileName , dim , &parts , &N , &periodic );
#endif

    if ( myrank == 0 )
        message( "reading particle properties took %.3f ms." , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);
    
    /* Apply h scaling */
    if( scaling != 1.0 )
      for ( k = 0 ; k < N ; k++ )
	    parts[k].h *= scaling;

    /* Apply shift */
    if(shift[0] !=0 || shift[1] !=0 || shift[2] !=0 )
      for ( k = 0 ; k < N ; k++ ) {
	    parts[k].x[0] += shift[0];
	    parts[k].x[1] += shift[1];
	    parts[k].x[2] += shift[2];
      }

    /* printParticle( parts , 10312237508790 , N );
    printParticle( parts , 10312286091950 , N ); */
    /* for ( k = 0 ; k < N ; k++ )
	if ( parts[k].id == 10312286091950 ||
	     parts[k].id == 10286889371446  ||
	     parts[k].id == 9536045071298 ||
	     parts[k].id == 12726778692106  ||
	     parts[k].id == 9479892852626  ||
	     parts[k].id == 9535843125514  ||
	     parts[k].id == 14151507889834  ||
	     parts[k].id == 14144038209438  ||
	     parts[k].id == 14121890205050  ||
	     parts[k].id == 5868762382714  ||
	     parts[k].id == 12840527117206 ||
	     parts[k].id == 10292087642778  ||
	     parts[k].id == 9465178320650  ||
	     parts[k].id == 2834846537770  ||
	     parts[k].id == 9483000048314  ||
	     parts[k].id == 10247332828902  ||
	     parts[k].id == 10223834653674  ||
	     parts[k].id == 16719632108962  ||
	     parts[k].id == 16759192850622  ||
	     parts[k].id == 9483599082554  ||
	     parts[k].id == 10247340329226 )
            parts[k] = parts[--N]; */

    /* Set default number of queues. */
    if ( nr_queues < 0 )
        nr_queues = nr_threads;
            
    /* Initialize the space with this data. */
    tic = getticks();
    space_init( &s , dim , parts , N , periodic , h_max );
    if ( myrank == 0 )
        message( "space_init took %.3f ms." , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);

    /* Set the default time step to 1.0f. */
    if ( myrank == 0 )
        message( "dt_max is %e." , dt_max );
    
    /* Say a few nice things about the space we just created. */
    if ( myrank == 0 ) {
        message( "space dimensions are [ %.3f %.3f %.3f ]." , s.dim[0] , s.dim[1] , s.dim[2] );
        message( "space %s periodic." , s.periodic ? "is" : "isn't" );
        message( "highest-level cell dimensions are [ %i %i %i ]." , s.cdim[0] , s.cdim[1] , s.cdim[2] );
        message( "%i parts in %i cells." , s.nr_parts , s.tot_cells );
        message( "maximum depth is %d." , s.maxdepth );
        // message( "cutoffs in [ %g %g ]." , s.h_min , s.h_max ); fflush(stdout);
        }
    
    /* Verify that each particle is in it's propper cell. */
    if ( myrank == 0 ) {
        icount = 0;
        space_map_cells_pre( &s , 0 , &map_cellcheck , &icount );
        message( "map_cellcheck picked up %i parts." , icount );
        }
    
    if ( myrank == 0 ) {
        data[0] = s.maxdepth; data[1] = 0;
        space_map_cells_pre( &s , 0 , &map_maxdepth , data );
        message( "nr of cells at depth %i is %i." , data[0] , data[1] );
        }
    
    /* Dump the particle positions. */
    // space_map_parts( &s , &map_dump , shift );
    
    
    /* Initialize the engine with this space. */
    tic = getticks();
    message( "nr_nodes is %i." , nr_nodes );
    engine_init( &e , &s , dt_max , nr_threads , nr_queues , nr_nodes , myrank , ENGINE_POLICY | engine_policy_steal );
    if ( myrank == 0 )
        message( "engine_init took %.3f ms." , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);

#ifdef WITH_MPI
    /* Split the space. */
    engine_split( &e , grid );
    engine_redistribute ( &e );
#endif

    message("Before write !");
    /* Write the state of the system as it is before starting time integration. */
    tic = getticks();
#if defined( WITH_MPI ) 
#if defined( HAVE_PARALLEL_HDF5 )
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
    message( "writing particle properties took %.3f ms." , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);
    
    /* Init the runner history. */
    #ifdef HIST
    for ( k = 0 ; k < runner_hist_N ; k++ )
        runner_hist_bins[k] = 0;
    #endif
    
    /* Inauguration speech. */
    if ( runs < INT_MAX )
        message( "starting for %i steps with %i threads and %i queues..." , runs , e.nr_threads , e.sched.nr_queues );
    else
        message( "starting for t=%.3e with %i threads and %i queues..." , clock , e.nr_threads , e.sched.nr_queues );
    fflush(stdout);
    
    /* Set a target particle. */
    /* long long int pid[5];
    unsigned int seed = 6178;
    for ( k = 0 ; k < 5 ; k++ )
        pid[k] = s.gparts[ rand_r( &seed ) % N ].part->id;
    for ( k = 0 ; k < 5 ; k++ )
        pairs_single_grav( dim , pid[k] , s.gparts , N , 0 ); */
    
    /* Legend. */
    if ( myrank == 0 )
        printf( "# step time e_tot e_kin e_temp dt dt_step count dt_min dt_max\n" );
    
    /* Let loose a runner on the space. */
    for ( j = 0 ; j < runs && e.time < clock ; j++ ) {
    
        /* Repartition the space amongst the nodes? */
        #if defined(WITH_MPI) && defined(HAVE_METIS)
            if ( j % 100 == 2 )
                e.forcerepart = 1;
        #endif
        
        /* Force a rebuild for testing. */
        /* if ( j % 4 == 3 )
            e.forcerebuild = 1; */
        
        // message( "starting run %i/%i (t=%.3e) with %i threads and %i queues..." , j+1 , runs , e.time , e.nr_threads , e.nr_queues ); fflush(stdout);
        timers_reset( timers_mask_all );
        #ifdef COUNTER
            for ( k = 0 ; k < runner_counter_count ; k++ )
                runner_counter[k] = 0;
        #endif
        
        /* Take a step. */
        engine_step( &e );
        
        if ( j % 100 == 0 )
	  {

#if defined( WITH_MPI ) 
#if defined( HAVE_PARALLEL_HDF5 )
	    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
	    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
	    write_output_single(&e, &us);
#endif

          }
                
        /* Dump a line of agregate output. */
        if ( myrank == 0 ) {
            printf( "%i %e %.16e %.16e %.16e %.3e %.3e %i %.3e %.3e" ,
                j , e.time ,
                e.ekin+e.epot , e.ekin , e.epot ,
                e.dt , e.dt_step , e.count_step ,
                e.dt_min , e.dt_max );
            for ( k = 0 ; k < timer_count ; k++ )
                printf( " %.3f" , ((double)timers[k])/CPU_TPS*1000 );
            printf( "\n" ); fflush(stdout);
            }
        /* for ( k = 0 ; k < 5 ; k++ )
            printgParticle( s.gparts , pid[k] , N ); */
        
    }
        
    /* Print the values of the runner histogram. */
#ifdef HIST
    printf( "main: runner histogram data:\n" );
    for ( k = 0 ; k < runner_hist_N ; k++ )
      printf( " %e %e %e\n" ,
	      runner_hist_a + k * (runner_hist_b - runner_hist_a) / runner_hist_N ,
	      runner_hist_a + (k + 1) * (runner_hist_b - runner_hist_a) / runner_hist_N ,
	      (double)runner_hist_bins[k] );
#endif

    /* Loop over the parts directly. */
    // for ( k = 0 ; k < N ; k++ )
    //     printf( " %i %e %e\n" , s.parts[k].id , s.parts[k].count , s.parts[k].count_dh );
    
    /* Dump the task data. */
    /* #ifdef WITH_MPI
    for ( j = 0 ; j < nr_nodes ; j++ ) {
    MPI_Barrier( MPI_COMM_WORLD );
    if ( j == myrank ) {
    printf( " %03i 0 0 0 0 %lli 0 0 0 0\n" , myrank , e.tic_step );
    for ( k = 0 ; k < e.sched.nr_tasks ; k++ )
        if ( !e.sched.tasks[k].skip && !e.sched.tasks[k].implicit )
            printf( " %03i %i %i %i %i %lli %lli %i %i %i\n" ,
                myrank ,
                e.sched.tasks[k].rid , e.sched.tasks[k].type , e.sched.tasks[k].subtype ,
                (e.sched.tasks[k].cj == NULL) , e.sched.tasks[k].tic , e.sched.tasks[k].toc ,
		e.sched.tasks[k].ci->count , (e.sched.tasks[k].cj!=NULL)?e.sched.tasks[k].cj->count:0 , e.sched.tasks[k].flags); 
    fflush(stdout);
    sleep(1);
    }
    }
    #else
    for ( k = 0 ; k < e.sched.nr_tasks ; k++ )
        if ( !e.sched.tasks[k].skip && !e.sched.tasks[k].implicit )
                printf( " %i %i %i %i %lli %lli %i %i\n" ,
                e.sched.tasks[k].rid , e.sched.tasks[k].type , e.sched.tasks[k].subtype , 
                (e.sched.tasks[k].cj == NULL) , e.sched.tasks[k].tic , e.sched.tasks[k].toc ,
                e.sched.tasks[k].ci->count , 
                (e.sched.tasks[k].cj==NULL)?0:e.sched.tasks[k].cj->count ); 
    #endif */
    
    /* Write final output. */
#if defined( WITH_MPI ) 
#if defined( HAVE_PARALLEL_HDF5 )
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
        
#ifdef WITH_MPI
    if ( MPI_Finalize() != MPI_SUCCESS )
      error( "call to MPI_Finalize failed with error %i." , res );
#endif
    
    /* Say goodbye. */
    message( "done." );
    
    /* All is calm, all is bright. */
    return 0;
    
}
