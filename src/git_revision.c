/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "git_revision.h"

const char* git_revision( void )
{
    static const char *revision = GIT_REVISION;
    return revision;
}

const char* package_version( void )
{
    static char buf[256];
    static int initialised = 0;
    if ( ! initialised ) {
        sprintf( buf, "SWIFT version: %s, at revision: %s", 
                 PACKAGE_VERSION, GIT_REVISION );
        initialised = 1;
    }
    return buf;
}
