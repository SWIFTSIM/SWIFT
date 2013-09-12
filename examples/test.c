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
    #define CPU_TPS 2.67e9
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
 * @brief Read coordinates from a text file.
 *
 * @param fname The name of the coordinate file.
 * @param parts An array of #part in which to store the coordinates.
 * @param N The number of parts to read.
 */
 
void read_coords ( char *fname , struct part *parts , int N ) {

#ifdef HAVE_LIBZ
    gzFile fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open coordinate file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading coordinate file." );
        if ( sscanf( buff , "%lf %lf %lf" , &parts[k].x[0] , &parts[k].x[1] , &parts[k].x[2] ) != 3 ) {
            printf( "read_coords: failed to parse %ith entry.\n" , k );
            error( "Error parsing coordinate file." );
            }
        }
        
    /* Wrap it up. */
    gzclose( fd );
#else
    FILE *fd;
    int k;
    
    /* Open the given file. */
    if ( ( fd = fopen( fname , "r" ) ) == NULL )
        error( "Failed to open coordinate file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( fscanf( fd , "%lf %lf %lf" , &parts[k].x[0] , &parts[k].x[1] , &parts[k].x[2] ) != 3 ) {
            printf( "read_coords: failed to read %ith entry.\n" , k );
            error( "Error reading coordinate file." );
            }
        }
        
    /* Wrap it up. */
    fclose( fd );
#endif

    }


/**
 * @brief Read cutoffs from a text file.
 *
 * @param fname The name of the cutoffs file.
 * @param parts An array of #part in which to store the cutoffs.
 * @param N The number of parts to read.
 */
 
void read_cutoffs ( char *fname , struct part *parts , int N ) {

#ifdef HAVE_LIBZ
    gzFile fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open cutoff file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading cutoff file." );
        if ( sscanf( buff , "%ef" , &parts[k].h ) != 1 ) {
            printf( "read_cutoffs: failed to parse %ith entry.\n" , k );
            error( "Error parsing cutoff file." );
            }
        }
        
    /* Wrap it up. */
    gzclose( fd );
#else
    FILE *fd;
    int k;
    
    /* Open the given file. */
    if ( ( fd = fopen( fname , "r" ) ) == NULL )
        error( "Failed to open cutoff file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( fscanf( fd , "%ef" , &parts[k].h ) != 1 ) {
            printf( "read_cutoffs: failed to read %ith entry.\n" , k );
            error( "Error reading cutoff file." );
            }
        }
        
    /* Wrap it up. */
    fclose( fd );
#endif

    }
    
    
/**
 * @brief Read id from a text file.
 *
 * @param fname The name of the id file.
 * @param parts An array of #part in which to store the dt.
 * @param N The number of parts to read.
 */
 
void read_id ( char *fname , struct part *parts , int N ) {

#ifdef HAVE_LIBZ
    gzFile fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open id file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading id file." );
        if ( sscanf( buff , "%lli" , &parts[k].id ) != 1 ) {
            printf( "read_id: failed to parse %ith entry.\n" , k );
            error( "Error parsing id file." );
            }
        }
        
    /* Wrap it up. */
    gzclose( fd );
#else
    FILE *fd;
    int k;
    
    /* Open the given file. */
    if ( ( fd = fopen( fname , "r" ) ) == NULL )
        error( "Failed to open id file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( fscanf( fd , "%lli" , &parts[k].id ) != 1 ) {
            printf( "read_id: failed to read %ith entry.\n" , k );
            error( "Error reading id file." );
            }
        }

    /* Wrap it up. */
    fclose( fd );
#endif

    }
    
    
/**
 * @brief Read dt from a text file.
 *
 * @param fname The name of the dt file.
 * @param parts An array of #part in which to store the dt.
 * @param N The number of parts to read.
 */
 
void read_dt ( char *fname , struct part *parts , int N ) {

#ifdef HAVE_LIBZ
    gzFile fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open dt file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading id file." );
        if ( sscanf( buff , "%f" , &parts[k].dt ) != 1 )
            error( "Error parsing dt file." );
        }
        
    /* Wrap it up. */
    gzclose( fd );
#else
    FILE *fd;
    int k;
    
    /* Open the given file. */
    if ( ( fd = fopen( fname , "r" ) ) == NULL )
        error( "Failed to open dt file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( fscanf( fd , "%ef" , &parts[k].dt ) != 1 )
            error( "Error reading dt file." );
        }

    /* Wrap it up. */
    fclose( fd );
#endif

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


void pairs_single ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {

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


/**
 * @brief Find the pairs of a single particle
 *
 * @param dim The space dimensions.
 * @param parts The #part array.
 * @param N The number of parts.
 * @param periodic Periodic boundary conditions flag.
 * @param target the index of the target particle.
 */

void pairs_single_old ( double *dim , struct part *__restrict__ parts , int N , int periodic , int target ) {

    int i, k, tid;
    double r, tx[3], th, dx[3];
    
    /* Get the target position and radius. */
    for ( k = 0 ; k < 3 ; k++ )
        tx[k] = parts[target].x[k];
    th = parts[target].h;
    tid = parts[target].id;
    
    /* Loop over all particle pairs. */
    #pragma omp parallel for schedule(dynamic), default(none), private(k,i,dx,r), shared(target,tx,th,tid,periodic,parts,dim,N)
    for ( k = 0 ; k < N ; k++ ) {
        if ( k == target )
            continue;
        for ( i = 0 ; i < 3 ; i++ ) {
            dx[i] = tx[i] - parts[k].x[i];
            if ( periodic ) {
                if ( dx[i] < -dim[i]/2 )
                    dx[i] += dim[i];
                else if ( dx[i] > dim[i]/2 )
                    dx[i] -= dim[i];
                }
            }
        r = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
        if ( r < th )
            printf( "pairs_single: %i %lli [%e,%e,%e] %e\n" ,
                tid , parts[k].id , dx[0] , dx[1] , dx[2] , r );
        }
            
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

    int c, icount, j, k, N = 100, periodic = 1;
    int nr_threads = 1, nr_queues = -1, runs = INT_MAX;
    int data[2];
    double dim[3] = { 1.0 , 1.0 , 1.0 }, shift[3] = { 0.0 , 0.0 , 0.0 };
    double h_max = -1.0 , scaling = 1.0;
    double clock = DBL_MAX;
    struct part *parts = NULL;
    struct space s;
    struct engine e;
    char ICfileName[200];
    float dt_max = 0.0f;
    ticks tic;
    int nr_nodes = 1, myrank = 0, grid[3] = { 1 , 1 , 1 };
    
    /* Choke on FP-exceptions. */
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
    
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
	  if ( myrank == 0 )
        message("IC to be read from file '%s'.", ICfileName);
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
    if ( myrank == 0 )
        message( "sizeof(struct part) is %li bytes." , (long int)sizeof( struct part ));

    /* Read particles and space information from (GADGET) IC */
    tic = getticks();
#ifdef WITH_MPI
    read_ic_parallel( ICfileName , dim , &parts , &N , &periodic );
#else
    read_ic( ICfileName , dim , &parts , &N , &periodic );
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
    engine_init( &e , &s , dt_max , nr_threads , nr_queues , nr_nodes , myrank , ENGINE_POLICY | engine_policy_steal );
    if ( myrank == 0 )
        message( "engine_init took %.3f ms." , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);

#ifdef WITH_MPI
    /* Split the space. */
    if ( nr_nodes != grid[0]*grid[1]*grid[2] )
        error( "Grid size does not match number of nodes." );
    engine_split( &e , grid );
#endif

    /* Write the state of the system as it is before starting time integration. */
    tic = getticks();
#ifdef WITH_MPI
    write_output_parallel(&e, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    write_output(&e);
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
    
    /* Legend. */
    if ( myrank == 0 )
        printf( "# step time e_tot e_kin e_temp dt dt_step count dt_min dt_max\n" );
    
    /* Let loose a runner on the space. */
    for ( j = 0 ; j < runs && e.time < clock ; j++ ) {
    
        /* Repartition the space amongst the nodes? */
        #if defined(WITH_MPI) && defined(HAVE_METIS)
            if ( j == 1 )
                e.forcerepart = 1;
        #endif
        
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
#ifdef WITH_MPI
	    write_output_parallel(&e, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
             write_output(&e);
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
    /* for ( k = 0 ; k < e.sched.nr_tasks ; k++ )
        if ( !e.sched.tasks[k].skip && !e.sched.tasks[k].implicit )
            printf( " %i %i %i %i %lli %lli %i %i\n" ,
                e.sched.tasks[k].rid , e.sched.tasks[k].type , e.sched.tasks[k].subtype , 
                (e.sched.tasks[k].cj == NULL) , e.sched.tasks[k].tic , e.sched.tasks[k].toc ,
                e.sched.tasks[k].ci->count , 
                (e.sched.tasks[k].cj==NULL)?0:e.sched.tasks[k].cj->count ); */
    
    /* Write final output. */
#ifdef WITH_MPI
	write_output_parallel( &e, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL );
#else
	write_output( &e );
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
