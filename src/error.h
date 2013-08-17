/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include <stdio.h>


/**
 * @brief Error macro. Prints the message given in argument and aborts.
 *
 */
#ifdef WITH_MPI
    extern int engine_rank;
    #define error(s, ...) { fprintf( stderr , "[%03i] %s:%s():%i: " s "\n" , engine_rank , __FILE__ , __FUNCTION__ , __LINE__ , ##__VA_ARGS__ ); abort(); }
#else
    #define error(s, ...) { fprintf( stderr , "%s:%s():%i: " s "\n" , __FILE__ , __FUNCTION__ , __LINE__ , ##__VA_ARGS__ ); abort(); }
#endif


/**
 * @brief Macro to print a localized message with variable arguments.
 *
 */
#ifdef WITH_MPI
    extern int engine_rank;
    #define message(s, ...) printf( "%s[%03i]: " s "\n" , __FUNCTION__ , engine_rank , ##__VA_ARGS__ )
#else
    #define message(s, ...) printf( "%s: " s "\n" , __FUNCTION__ , ##__VA_ARGS__ )
#endif
