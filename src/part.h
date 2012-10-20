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


/* Data of a single particle. */
struct part {

    /* Particle position. */
    double x[3];
    
    /* Particle cutoff radius. */
    float r;
    
    /* Particle time-step. */
    float dt;
    
    /* Particle ID. */
    int id;
    
    /* Number of pairwise interactions. */
    double count, count_dh;
    int icount;
    
    } __attribute__((aligned (32)));
    

