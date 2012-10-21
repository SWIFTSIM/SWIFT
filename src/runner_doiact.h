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


/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#define PASTE(x,y) x ## _ ## y

#define DOPAIR2(f) PASTE(runner_dopair,f)
#define DOPAIR DOPAIR2(FUNCTION)

#define DOPAIR_NAIVE2(f) PASTE(runner_dopair_naive,f)
#define DOPAIR_NAIVE DOPAIR_NAIVE2(FUNCTION)

#define DOSELF2(f) PASTE(runner_doself,f)
#define DOSELF DOSELF2(FUNCTION)

#define DOSUB2(f) PASTE(runner_dosub,f)
#define DOSUB DOSUB2(FUNCTION)

#define IACT2(f) PASTE(runner_iact,f)
#define IACT IACT2(FUNCTION)


/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
 
void DOPAIR_NAIVE ( struct runner_thread *rt , struct cell *ci , struct cell *cj ) {

    struct runner *r = rt->r;
    int pid, pjd, k, count_i = ci->count, count_j = cj->count;
    double shift[3] = { 0.0 , 0.0 , 0.0 };
    struct part *pi, *pj, *parts_i = ci->parts, *parts_j = cj->parts;
    double dx[3], pix[3], hi, hi2, r2;
    TIMER_TIC
    
    /* Get the relative distance between the pairs, wrapping. */
    for ( k = 0 ; k < 3 ; k++ ) {
        if ( cj->loc[k] - ci->loc[k] < -r->s->dim[k]/2 )
            shift[k] = r->s->dim[k];
        else if ( cj->loc[k] - ci->loc[k] > r->s->dim[k]/2 )
            shift[k] = -r->s->dim[k];
        }
        
    /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with %i/%i parts and shift = [ %g %g %g ].\n" ,
        ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] , cj->loc[2] ,
        ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
    tic = getticks(); */
    
    /* Loop over the parts in ci. */
    for ( pid = 0 ; pid < count_i ; pid++ ) {
    
        /* Get a hold of the ith part in ci. */
        pi = &parts_i[ pid ];
        for ( k = 0 ; k < 3 ; k++ )
            pix[k] = pi->x[k] - shift[k];
        hi = pi->h;
        hi2 = hi * hi;
        
        /* Loop over the parts in cj. */
        for ( pjd = 0 ; pjd < count_j ; pjd++ ) {
        
            /* Get a pointer to the jth particle. */
            pj = &parts_j[ pjd ];
        
            /* Compute the pairwise distance. */
            r2 = 0.0;
            for ( k = 0 ; k < 3 ; k++ ) {
                dx[k] = pix[k] - pj->x[k];
                r2 += dx[k]*dx[k];
                }
                
            /* Hit or miss? */
            if ( r2 < hi2 || r2 < pj->h*pj->h ) {
            
                IACT( r2 , hi , pj->h , pi , pj );
            
                }
        
            } /* loop over the parts in cj. */
    
        } /* loop over the parts in ci. */
        
    #ifdef TIMER_VERBOSE
        printf( "runner_dopair_naive[%02i]: %i/%i parts at depth %i (r_max=%.3f/%.3f) took %.3f ms.\n" , rt->id , count_i , count_j , ci->depth , ci->r_max , cj->r_max , ((double)TIMER_TOC(runner_timer_dopair)) / CPU_TPS * 1000 );
    #else
        TIMER_TOC(runner_timer_dopair);
    #endif


    }


/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
 
void DOPAIR ( struct runner_thread *rt , struct cell *ci , struct cell *cj ) {

    struct runner *r = rt->r;
    int pid, pjd, k, sid;
    double rshift, shift[3] = { 0.0 , 0.0 , 0.0 };
    struct cell *temp;
    struct entry *sort_i, *sort_j;
    struct part *pi, *pj, *parts_i, *parts_j;
    double dx[3], pix[3], pjx[3], hi, hi2, hj, hj2, r2, di, dj;
    double hi_max, hj_max, di_max, dj_min;
    int count_i, count_j;
    TIMER_TIC
    
    /* Get the relative distance between the pairs, wrapping. */
    for ( k = 0 ; k < 3 ; k++ ) {
        if ( cj->loc[k] - ci->loc[k] < -r->s->dim[k]/2 )
            shift[k] = r->s->dim[k];
        else if ( cj->loc[k] - ci->loc[k] > r->s->dim[k]/2 )
            shift[k] = -r->s->dim[k];
        }
        
    /* Get the sorting index. */
    for ( sid = 0 , k = 0 ; k < 3 ; k++ )
        sid = 3*sid + ( (cj->loc[k] - ci->loc[k] + shift[k] < 0) ? 0 : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1 );

    /* Switch the cells around? */
    if ( runner_flip[sid] ) {
        temp = ci; ci = cj; cj = temp;
        for ( k = 0 ; k < 3 ; k++ )
            shift[k] = -shift[k];
        }
    sid = sortlistID[sid];
    
    /* Get the cutoff shift. */
    for ( rshift = 0.0 , k = 0 ; k < 3 ; k++ )
        rshift += shift[k]*runner_shift[ 3*sid + k ];
        
    /* printf( "runner_dopair: doing pair [ %g %g %g ]/[ %g %g %g ] with %i/%i parts and shift = [ %g %g %g ].\n" ,
        ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] , cj->loc[2] ,
        ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout); */
    /* for ( hi = 0 , k = 0 ; k < ci->count ; k++ )
        hi += ci->parts[k].r;
    for ( hj = 0 , k = 0 ; k < cj->count ; k++ )
        hj += cj->parts[k].r;
    printf( "runner_dopair: avg. radii %g/%g for h=%g at depth=%i.\n" , hi/ci->count , hj/cj->count , ci->h[0] , ci->depth ); fflush(stdout); */
    
    /* Pick-out the sorted lists. */
    sort_i = &ci->sort[ sid*(ci->count + 1) ];
    sort_j = &cj->sort[ sid*(cj->count + 1) ];
    
    /* Get some other useful values. */
    hi_max = ci->r_max - rshift; hj_max = cj->r_max - rshift;
    count_i = ci->count; count_j = cj->count;
    parts_i = ci->parts; parts_j = cj->parts;
    di_max = sort_i[count_i-1].d - rshift;
    dj_min = sort_j[0].d;
    
    /* if ( ci->split && cj->split && sid == 4 )
        printf( "boing!\n" ); */
    
    /* Loop over the parts in ci. */
    for ( pid = count_i-1 ; pid >= 0 && sort_i[pid].d + hi_max > dj_min ; pid-- ) {
    
        /* Get a hold of the ith part in ci. */
        pi = &parts_i[ sort_i[ pid ].i ];
        hi = pi->h;
        di = sort_i[pid].d + hi - rshift;
        if ( di < dj_min )
            continue;
            
        hi2 = pi->h * pi->h;
        for ( k = 0 ; k < 3 ; k++ )
            pix[k] = pi->x[k] - shift[k];
        
        /* Loop over the parts in cj. */
        for ( pjd = 0 ; pjd < count_j && sort_j[pjd].d < di ; pjd++ ) {
        
            /* Get a pointer to the jth particle. */
            pj = &parts_j[ sort_j[pjd].i ];
        
            /* Compute the pairwise distance. */
            r2 = 0.0;
            for ( k = 0 ; k < 3 ; k++ ) {
                dx[k] = pix[k] - pj->x[k];
                r2 += dx[k]*dx[k];
                }
                
            /* Hit or miss? */
            if ( r2 < hi2 ) {
            
                IACT( r2 , hi , pj->h , pi , pj );
            
                }
        
            } /* loop over the parts in cj. */
    
        } /* loop over the parts in ci. */
        
    /* printf( "runner_dopair: first half took %.3f ms...\n" , ((double)(getticks() - tic)) / CPU_TPS * 1000 );
    tic = getticks(); */

    /* Loop over the parts in cj. */
    for ( pjd = 0 ; pjd < count_j && sort_j[pjd].d - hj_max < di_max ; pjd++ ) {
    
        /* Get a hold of the jth part in cj. */
        pj = &parts_j[ sort_j[ pjd ].i ];
        hj = pj->h;
        dj = sort_j[pjd].d - hj - rshift;
        if ( dj > di_max )
            continue;
            
        for ( k = 0 ; k < 3 ; k++ )
            pjx[k] = pj->x[k] + shift[k];
        hj2 = pj->h * pj->h;
        
        /* Loop over the parts in ci. */
        for ( pid = count_i-1 ; pid >= 0 && sort_i[pid].d > dj ; pid-- ) {
        
            /* Get a pointer to the jth particle. */
            pi = &parts_i[ sort_i[pid].i ];
            
            /* Compute the pairwise distance. */
            r2 = 0.0;
            for ( k = 0 ; k < 3 ; k++ ) {
                dx[k] = pi->x[k] - pjx[k];
                r2 += dx[k]*dx[k];
                }
                
            /* Hit or miss? */
            if ( r2 < hj2 && r2 > pi->h*pi->h ) {
            
                IACT( r2 , pi->h , hj , pi , pj );
            
                }
        
            } /* loop over the parts in cj. */
    
        } /* loop over the parts in ci. */

    #ifdef TIMER_VERBOSE
        printf( "runner_dopair[%02i]: %i/%i parts at depth %i (r_max=%.3f/%.3f, h=%.3f) took %.3f ms.\n" , rt->id , count_i , count_j , ci->depth , ci->r_max , cj->r_max , fmax(ci->h[0],fmax(ci->h[1],ci->h[2])) , ((double)(TIMER_TOC(runner_timer_dopair))) / CPU_TPS * 1000 );
    #else
        TIMER_TOC(runner_timer_dopair);
    #endif

    }


/**
 * @brief Compute the cell self-interaction.
 *
 * @param r The #runner.
 * @param c The #cell.
 */

void DOSELF ( struct runner_thread *rt , struct cell *c ) {

    int k, pid, pjd, count = c->count;
    double pix[3], dx[3], hi, hi2, r2;
    struct part *pi, *pj, *parts = c->parts;
    TIMER_TIC
    
    if ( c->split )
        error( "Split cell should not have self-interactions." );
    
    /* Loop over the particles in the cell. */
    for ( pid = 0 ; pid < count ; pid++ ) {
    
        /* Get a pointer to the ith particle. */
        pi = &parts[pid];
    
        /* Get the particle position and radius. */
        for ( k = 0 ; k < 3 ; k++ )
            pix[k] = pi->x[k];
        hi = pi->h;
        hi2 = hi * hi;
            
        /* Loop over the other particles .*/
        for ( pjd = pid+1 ; pjd < count ; pjd++ ) {
        
            /* Get a pointer to the jth particle. */
            pj = &parts[pjd];
        
            /* Compute the pairwise distance. */
            r2 = 0.0;
            for ( k = 0 ; k < 3 ; k++ ) {
                dx[k] = pix[k] - pj->x[k];
                r2 += dx[k]*dx[k];
                }
                
            /* Hit or miss? */
            if ( r2 < hi2 || r2 < pj->h*pj->h ) {
            
                IACT( r2 , hi , pj->h , pi , pj );
            
                }
        
            } /* loop over all other particles. */
    
        } /* loop over all particles. */

    #ifdef TIMER_VERBOSE
        printf( "runner_doself[%02i]: %i parts at depth %i took %.3f ms.\n" , rt->id , count , c->depth , ((double)TIMER_TOC(runner_timer_doself)) / CPU_TPS * 1000 );
    #else
        TIMER_TOC(runner_timer_doself);
    #endif

    }


/**
 * @brief Compute grouped sub-cell interactions
 *
 * @param r The #runner.
 * @param c The #cell.
 */

void DOSUB ( struct runner_thread *rt , struct cell *ci , struct cell *cj , int flags ) {

    int j, k;

    TIMER_TIC
    
    /* Different types of flags. */
    switch ( flags ) {
    
        /* Regular sub-cell interactions of a single cell. */
        case 0:
            for ( j = 0 ; j < 7 ; j++ )
                for ( k = j + 1 ; k < 8 ; k++ )
                    if ( ci->progeny[j] != NULL && ci->progeny[k] != NULL )
                        DOPAIR( rt , ci->progeny[j] , ci->progeny[k] );
            break;
            
        case 1: /* (  1 ,  1 ,  0 ) */
            if ( ci->progeny[6] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[0] );
            if ( ci->progeny[6] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[1] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[1] );
            break;
    
        case 3: /* (  1 ,  0 ,  1 ) */
            if ( ci->progeny[5] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[0] );
            if ( ci->progeny[5] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[2] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[2] );
            break;
                    
        case 4: /* (  1 ,  0 ,  0 ) */
            if ( ci->progeny[4] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[0] );
            if ( ci->progeny[4] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[1] );
            if ( ci->progeny[4] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[2] );
            if ( ci->progeny[4] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[3] );
            if ( ci->progeny[5] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[0] );
            if ( ci->progeny[5] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[1] );
            if ( ci->progeny[5] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[2] );
            if ( ci->progeny[5] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[3] );
            if ( ci->progeny[6] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[0] );
            if ( ci->progeny[6] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[1] );
            if ( ci->progeny[6] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[2] );
            if ( ci->progeny[6] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[3] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[1] );
            if ( ci->progeny[7] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[2] );
            if ( ci->progeny[7] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[3] );
            break;
            
        case 5: /* (  1 ,  0 , -1 ) */
            if ( ci->progeny[4] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[1] );
            if ( ci->progeny[4] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[3] );
            if ( ci->progeny[6] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[1] );
            if ( ci->progeny[6] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[3] );
            break;
                    
        case 7: /* (  1 , -1 ,  0 ) */
            if ( ci->progeny[4] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[2] );
            if ( ci->progeny[4] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[4] , cj->progeny[3] );
            if ( ci->progeny[5] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[2] );
            if ( ci->progeny[5] != NULL && cj->progeny[3] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[3] );
            break;
                    
        case 9: /* (  0 ,  1 ,  1 ) */
            if ( ci->progeny[3] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[0] );
            if ( ci->progeny[3] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[4] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[4] );
            break;
                    
        case 10: /* (  0 ,  1 ,  0 ) */
            if ( ci->progeny[2] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[0] );
            if ( ci->progeny[2] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[1] );
            if ( ci->progeny[2] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[4] );
            if ( ci->progeny[2] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[5] );
            if ( ci->progeny[3] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[0] );
            if ( ci->progeny[3] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[1] );
            if ( ci->progeny[3] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[4] );
            if ( ci->progeny[3] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[5] );
            if ( ci->progeny[6] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[0] );
            if ( ci->progeny[6] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[1] );
            if ( ci->progeny[6] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[4] );
            if ( ci->progeny[6] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[5] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[1] );
            if ( ci->progeny[7] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[4] );
            if ( ci->progeny[7] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[5] );
            break;
                    
        case 11: /* (  0 ,  1 , -1 ) */
            if ( ci->progeny[2] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[1] );
            if ( ci->progeny[2] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[2] , cj->progeny[5] );
            if ( ci->progeny[6] != NULL && cj->progeny[1] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[1] );
            if ( ci->progeny[6] != NULL && cj->progeny[5] != NULL )
                DOPAIR( rt , ci->progeny[6] , cj->progeny[5] );
            break;
                    
        case 12: /* (  0 ,  0 ,  1 ) */
            if ( ci->progeny[1] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[1] , cj->progeny[0] );
            if ( ci->progeny[1] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[1] , cj->progeny[2] );
            if ( ci->progeny[1] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[1] , cj->progeny[4] );
            if ( ci->progeny[1] != NULL && cj->progeny[6] != NULL )
                DOPAIR( rt , ci->progeny[1] , cj->progeny[6] );
            if ( ci->progeny[3] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[0] );
            if ( ci->progeny[3] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[2] );
            if ( ci->progeny[3] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[4] );
            if ( ci->progeny[3] != NULL && cj->progeny[6] != NULL )
                DOPAIR( rt , ci->progeny[3] , cj->progeny[6] );
            if ( ci->progeny[5] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[0] );
            if ( ci->progeny[5] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[2] );
            if ( ci->progeny[5] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[4] );
            if ( ci->progeny[5] != NULL && cj->progeny[6] != NULL )
                DOPAIR( rt , ci->progeny[5] , cj->progeny[6] );
            if ( ci->progeny[7] != NULL && cj->progeny[0] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[0] );
            if ( ci->progeny[7] != NULL && cj->progeny[2] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[2] );
            if ( ci->progeny[7] != NULL && cj->progeny[4] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[4] );
            if ( ci->progeny[7] != NULL && cj->progeny[6] != NULL )
                DOPAIR( rt , ci->progeny[7] , cj->progeny[6] );
            break;
                
        }
    

    #ifdef TIMER_VERBOSE
        printf( "runner_dosub[%02i]: flags=%i at depth %i took %.3f ms.\n" , rt->id , flags , ci->depth , ((double)TIMER_TOC(runner_timer_dosub)) / CPU_TPS * 1000 );
    #else
        TIMER_TOC(runner_timer_dosub);
    #endif

    }


