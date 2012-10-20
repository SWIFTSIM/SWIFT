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
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <omp.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
    #include <zlib.h>
#endif

/* Local headers. */
#include "gadgetsmp.h"

/* Ticks per second on this machine. */
#ifndef CPU_TPS
    #define CPU_TPS 2.67e9
#endif

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/**
 * @brief Mapping function to draw a specific cell (gnuplot).
 */

void map_cells_plot ( struct cell *c , void *data ) {

    int k, depth = *(int *)data;
    double *l = c->loc, *h = c->h;

    if ( c->depth >= depth ) {
    
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
        
        for ( k = 0 ; k < c->count ; k++ )
            printf( "%.16e %.16e %.16e %.16e %.16e %.16e\n" , l[0]+h[0] , l[1]+h[1] , l[2] ,
                c->parts[k].x[0] , c->parts[k].x[1] , c->parts[k].x[2] );
        printf( "\n\n" );
    
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
    
    /* Loop over all parts and check if they are in the cell. */
    for ( k = 0 ; k < c->count ; k++ ) {
        *count += 1;
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

    double *count = (double *)data;
    
    // printf( "%i %e %e\n" , p->id , p->count , p->count_dh );

    *count += p->count;

    }


/**
 * @brief Mapping function for neighbour count.
 */

void map_icount ( struct part *p , struct cell *c , void *data ) {

    int *count = (int *)data;
    
    // printf( "%i\n" , p->icount );

    *count += p->icount;

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
    gzFile *fd;
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
    gzFile *fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open cutoff file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading cutoff file." );
        if ( sscanf( buff , "%ef" , &parts[k].r ) != 1 ) {
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
        if ( fscanf( fd , "%ef" , &parts[k].r ) != 1 ) {
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
    gzFile *fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open id file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading id file." );
        if ( sscanf( buff , "%i" , &parts[k].id ) != 1 ) {
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
        if ( fscanf( fd , "%i" , &parts[k].id ) != 1 ) {
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
    gzFile *fd;
    char buff[1024];
    int k;
    
    /* Open the given file. */
    if ( ( fd = gzopen( fname , "r" ) ) == NULL )
        error( "Failed to open dt file" );
        
    /* Read the coordinates into the part positions. */
    for ( k = 0 ; k < N ; k++ ) {
        if ( gzgets( fd , buff , 1024 ) == NULL )
            error( "Error reading id file." );
        if ( sscanf( buff , "%ef" , &parts[k].dt ) != 1 )
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

    int i, j, k, count = 0, mj, mk;
    double r2, dx[3], dcount = 0.0, maxratio = 1.0;
    float f1, f2;
    
    /* Loop over all particle pairs. */
    #pragma omp parallel for schedule(dynamic), default(none), private(k,i,dx,r2,f1,f2), shared(maxratio,mj,mk,periodic,parts,dim,N,stdout), reduction(+:count,dcount)
    for ( j = 0 ; j < N ; j++ ) {
        if ( j % 1000 == 0 ) {
            printf( "pairs_n2: j=%i, count=%i.\n" , j , count );
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
            if ( r2 < parts[j].r*parts[j].r || r2 < parts[k].r*parts[k].r ) {
                f1 = 0.0; f2 = 0.0;
                iact_nopart( r2 , parts[j].r , parts[k].r , &f1 , &f2 , &count , &count );
                dcount += f1 + f2;
                if ( parts[j].r / parts[k].r > maxratio )
                    #pragma omp critical
                    {
                    maxratio = parts[j].r / parts[k].r;
                    mj = j; mk = k;
                    }
                else if ( parts[k].r / parts[j].r > maxratio )
                    #pragma omp critical
                    {
                    maxratio = parts[k].r / parts[j].r;
                    mj = j; mk = k;
                    }
                }
            }
        }
            
    /* Dump the result. */
    printf( "pairs_n2: avg. nr. of pairs per part is %.3f (%.3f).\n" , ((double)count)/N , dcount/N + 32.0/3 );
    printf( "pairs_n2: maximum ratio between parts %i [%e,%e,%e] and %i [%e,%e,%e] is %.3f/%.3f\n" ,
        mj , parts[mj].x[0] , parts[mj].x[1] , parts[mj].x[2] ,
        mk , parts[mk].x[0] , parts[mk].x[1] , parts[mk].x[2] ,
        parts[mj].r , parts[mk].r ); fflush(stdout);
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

void pairs_single ( double *dim , struct part *__restrict__ parts , int N , int periodic , int target ) {

    int i, k, tid;
    double r, tx[3], tr, dx[3];
    
    /* Get the target position and radius. */
    for ( k = 0 ; k < 3 ; k++ )
        tx[k] = parts[target].x[k];
    tr = parts[target].r;
    tid = parts[target].id;
    
    /* Loop over all particle pairs. */
    #pragma omp parallel for schedule(dynamic), default(none), private(k,i,dx,r), shared(target,tx,tr,tid,periodic,parts,dim,N)
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
        if ( r < tr )
            printf( "pairs_single: %i %i [%e,%e,%e] %e\n" ,
                tid , parts[k].id , dx[0] , dx[1] , dx[2] , r );
        }
            
    }


/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
 
int main ( int argc , char *argv[] ) {

    int c, icount, j, k, N = 100, periodic = 1;
    int nr_threads = 1, nr_queues = -1, runs = 1;
    int data[2];
    double dim[3] = { 1.0 , 1.0 , 1.0 }, shift[3] = { 0.0 , 0.0 , 0.0 };
    double r_min = 0.01, r_max = 0.1, h_max = -1.0 , scaling = 1.0, count = 0.0;
    struct part *parts = NULL;
    struct space s;
    struct runner r;
    ticks tic;
    
    /* Init the space. */
    bzero( &s , sizeof(struct space) );
    
    /* Parse the options. */
    while ( ( c = getopt( argc , argv  , "a:b:p:d:N:c:h:v:m:s:t:q:r:i:m:z:" ) ) != -1 )
        switch ( c ) {
            case 'N':
                if ( sscanf( optarg , "%d" , &N ) != 1 )
                    error( "Error parsing number of particles." );
                if ( posix_memalign( (void *)&parts , 16 , N * sizeof(struct part) ) != 0 )
                    error( "Call to posix_memalign failed." );
                for ( k = 0 ; k < N ; k++ ) {
                    parts[k].x[0] = ((double)rand()) / RAND_MAX * dim[0];
                    parts[k].x[1] = ((double)rand()) / RAND_MAX * dim[1];
                    parts[k].x[2] = ((double)rand()) / RAND_MAX * dim[2];
                    parts[k].id = k;
                    parts[k].r = r_min + ((r_max - r_min)*rand())/RAND_MAX;
                    }
                printf( "main: allocated memory for %i parts.\n" , N ); fflush(stdout);
                break;
            case 'a':
                if ( sscanf( optarg , "%lf" , &scaling ) != 1 )
                    error( "Error parsing cutoff scaling." );
                printf( "main: scaling cutoff by %.3f.\n" , scaling ); fflush(stdout);
                for ( k = 0 ; k < N ; k++ )
                    parts[k].r *= scaling;
                break;
            case 'b':
                if ( sscanf( optarg , "%lf %lf %lf" , &dim[0] , &dim[1] , &dim[2] ) != 3 )
                    error( "Error parsing box dimensions." );
                break;
            case 'c':
                printf( "main: reading parts from %s...\n" , optarg ); fflush(stdout);
                if ( parts == NULL && posix_memalign( (void *)&parts , 16 , N * sizeof(struct part) ) != 0 )
                    error( "Call to calloc failed." );
                read_coords( optarg , parts , N );
                break;
            case 'd':
                printf( "main: reading dt from %s...\n" , optarg ); fflush(stdout);
                read_dt( optarg , parts , N );
                break;
            case 'h':
                printf( "main: reading cutoffs from %s...\n" , optarg ); fflush(stdout);
                read_cutoffs( optarg , parts , N );
                break;
            case 'i':
                printf( "main: reading ids from %s...\n" , optarg ); fflush(stdout);
                read_id( optarg , parts , N );
                break;
            case 'm':
                if ( sscanf( optarg , "%lf" , &h_max ) != 1 )
                    error( "Error parsing h_max." );
                printf( "main: maximum h set to %e.\n" , h_max );
                break;
            case 'p':
                if ( sscanf( optarg , "%d" , &periodic ) != 1 )
                    error( "Error parsing periodicity." );
                printf( "main: periodicity switched %s.\n" , periodic ? "on" : "off" );
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
                for ( k = 0 ; k < N ; k++ ) {
                    parts[k].x[0] += shift[0];
                    parts[k].x[1] += shift[1];
                    parts[k].x[2] += shift[2];
                    }
                printf( "main: shifted parts by [ %.3f %.3f %.3f ].\n" , shift[0] , shift[1] , shift[2] );
                break;
            case 't':
                if ( sscanf( optarg , "%d" , &nr_threads ) != 1 )
                    error( "Error parsing number of threads." );
                omp_set_num_threads( nr_threads );
                break;
            case 'z':
                if ( sscanf( optarg , "%d" , &space_splitsize ) != 1 )
                    error( "Error parsing split size." );
                printf( "main: split size set to %i.\n" , space_splitsize );
                break;
            case '?':
                error( "Unknown option." );
                break;
            }
    
    /* Get the brute-force number of pairs. */
    // pairs_n2( dim , parts , N , periodic );
    // pairs_single( dim , parts , N , periodic , 63628 );
    fflush( stdout );
    
    /* Set default number of queues. */
    if ( nr_queues < 0 )
        nr_queues = nr_threads;
            
    /* Initialize the space with this data. */
    tic = getticks();
    space_init( &s , dim , parts , N , periodic , h_max );
    printf( "main: space_init took %.3f ms.\n" , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);
    
    /* Say a few nice things about the space we just created. */
    printf( "main: space dimensions are [ %.3f %.3f %.3f ].\n" , s.dim[0] , s.dim[1] , s.dim[2] );
    printf( "main: space %s periodic.\n" , s.periodic ? "is" : "isn't" );
    printf( "main: highest-level cell dimensions are [ %i %i %i ].\n" , s.cdim[0] , s.cdim[1] , s.cdim[2] );
    printf( "main: %i parts in %i cells.\n" , s.nr_parts , s.tot_cells );
    printf( "main: maximum depth is %d.\n" , s.maxdepth );
    printf( "main: cutoffs in [ %.3f %.3f ].\n" , s.r_min , s.r_max ); fflush(stdout);
    
    /* Verify that each particle is in it's propper cell. */
    // icount = 0;
    // space_map_cells( &s , &map_cellcheck , &icount );
    // printf( "main: map_cellcheck picked up %i parts.\n" , icount );
    
    data[0] = s.maxdepth; data[1] = 0;
    space_map_cells( &s , &map_maxdepth , data );
    printf( "main: nr of cells at depth %i is %i.\n" , data[0] , data[1] );
    
    /* Dump the particle positions. */
    // space_map_parts( &s , &map_dump , shift );
    
    /* Generate the tasks. */
    tic = getticks();
    space_maketasks( &s , 1 );
    printf( "main: generated %i tasks.\n" , s.nr_tasks );
    printf( "main: space_maketasks took %.3f ms.\n" , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);
    
    /* Initialize the runner with this space. */
    tic = getticks();
    runner_init( &r , &s , nr_threads , nr_queues , runner_policy_steal | runner_policy_keep );
    printf( "main: runner_init took %.3f ms.\n" , ((double)(getticks() - tic)) / CPU_TPS * 1000 ); fflush(stdout);
    
    /* Init the runner history. */
    #ifdef HIST
    for ( k = 0 ; k < runner_hist_N ; k++ )
        runner_hist_bins[k] = 0;
    #endif
    
    /* Let loose a runner on the space. */
    for ( j = 0 ; j < runs ; j++ ) {
        printf( "main: starting run %i/%i with %i threads and %i queues...\n" , j+1 , runs , r.nr_threads , r.nr_queues ); fflush(stdout);
        tic = getticks();
        #ifdef TIMER
            for ( k = 0 ; k < runner_timer_count ; k++ )
                runner_timer[k] = 0;
        #endif
        #ifdef COUNTER
            for ( k = 0 ; k < runner_counter_count ; k++ )
                runner_counter[k] = 0;
        #endif
        runner_run( &r , 0 );
        #ifdef TIMER
            printf( "main: runner timers are [ %.3f" , runner_timer[0]/CPU_TPS*1000 );
            for ( k = 1 ; k < runner_timer_count ; k++ )
                printf( " %.3f" , ((double)runner_timer[k])/CPU_TPS*1000 );
            printf( " %.3f ] ms.\n" , ((double)(getticks() - tic)) / CPU_TPS * 1000 );
        #else
            printf( "main: runner_run with %i threads took %.3f ms.\n" , nr_threads , ((double)(getticks() - tic)) / CPU_TPS * 1000 );
        #endif
        #ifdef COUNTER
            printf( "main: runner counters are [ %d" , runner_counter[0] );
            for ( k = 1 ; k < runner_counter_count ; k++ )
                printf( " %d" , runner_counter[k] );
            printf( " ].\n" );
        #endif
        printf( "main: runner queue lengths are [ %i" , r.queues[0].count );
        for ( k = 1 ; k < r.nr_queues ; k++ )
            printf( " %i" , r.queues[k].count );
        printf( " ].\n" );
        fflush(stdout);
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
    
    /* Get the average interactions per particle. */
    count = 0;
    space_map_parts( &s , &map_count , &count );
    printf( "main: average interactions per particle is %.3f.\n" , count / s.nr_parts / runs + 32.0/3 );
    
    /* Get the average interactions per particle. */
    icount = 0;
    space_map_parts( &s , &map_icount , &icount );
    printf( "main: average neighbours per particle is %.3f.\n" , (double)icount / s.nr_parts / runs );
    
    /* Get all the cells of a certain depth. */
    /* count = 11;
    space_map_cells( &s , &map_cells_plot , &count ); */
    
    /* Check for outliers. */
    // space_map_parts( &s , &map_check , NULL );
    
    /* All is calm, all is bright. */
    return 0;
    
    }
