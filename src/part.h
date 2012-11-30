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


/* Some constants. */
#define part_maxwait                    3
#define part_maxunlock                  39
#define part_dtmax                      10


/* Data of a single particle. */
struct part {

    /* Particle cutoff radius. */
    float h;
    
    /* Particle time-step. */
    int dt;
    
    /* Particle mass. */
    float mass;
    
    /* Particle ID. */
    long id;
    
    /* Particle position. */
    double x[3];
    
    /* Particle velocity. */
    float v[3];
    
    /* Particle acceleration. */
    float a[3];
    
    /* Particle density. */
    float rho;
    
    /* Particle pressure. */
    // float P;
    
    /* Aggregate quantities. */
    float POrho2;
    
    /* Particle internal energy. */
    float u;
    
    /* Change in particle energy over time. */
    float u_dt;
    
    /* Change in smoothing length over time. */
    float h_dt;
    
    /* Derivative of the density with respect to this particle's smoothing length. */
    float rho_dh;
    
    /* Particle number density. */
    int icount;
    float wcount;
    float wcount_dh;
    
    } __attribute__((aligned (32)));
    

