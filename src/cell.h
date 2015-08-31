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
#ifndef SWIFT_CELL_H
#define SWIFT_CELL_H

/* Includes. */
#include "lock.h"
#include "part.h"
#include "multipole.h"

/* Forward declaration of space, needed for cell_unpack. */
struct space;

/* Some constants. */
#define cell_sid_dt 13
#define cell_max_tag (1 << 16)

/* Global variables. */
extern int cell_next_tag;

/* Packed cell. */
struct pcell {

  /* Stats on this cell's particles. */
  double h_max, dt_min, dt_max;

  /* Number of particles in this cell. */
  int count;

  /* tag used for MPI communication. */
  int tag;

  /* Relative indices of the cell's progeny. */
  int progeny[8];
};

/* Structure to store the data of a single cell. */
struct cell {

  /* The cell location on the grid. */
  double loc[3];

  /* The cell dimensions. */
  double h[3];

  /* Max radii in this cell. */
  double h_max;

  /* Minimum and maximum dt in this cell. */
  double dt_min, dt_max;

  /* Minimum dimension, i.e. smallest edge of this cell. */
  float dmin;

  /* Maximum slack allowed for particle movement. */
  float slack;

  /* Maximum particle movement in this cell. */
  float dx_max;

  /* The depth of this cell in the tree. */
  int depth, split, maxdepth;

  /* Nr of parts. */
  int count, gcount;

  /* Pointers to the particle data. */
  struct part *parts;

  /* Pointers to the extra particle data. */
  struct xpart *xparts;

  /* Pointers to the gravity particle data. */
  struct gpart *gparts;

  /* Pointers for the sorted indices. */
  struct entry *sort, *gsort;
  unsigned int sorted, gsorted;

  /* Pointers to the next level of cells. */
  struct cell *progeny[8];

  /* Parent cell. */
  struct cell *parent;

  /* Super cell, i.e. the highest-level supercell that has interactions. */
  struct cell *super;

  /* The task computing this cell's sorts. */
  struct task *sorts, *gsorts;
  int sortsize, gsortsize;

  /* The tasks computing this cell's density. */
  struct link *density, *force, *grav;
  int nr_density, nr_force, nr_grav;

  /* The ghost task to link density to interactions. */
  struct task *ghost, *kick1, *kick2;

  /* Task receiving data. */
  struct task *recv_xv, *recv_rho;

  /* Tasks for gravity tree. */
  struct task *grav_up, *grav_down;

  /* Number of tasks that are associated with this cell. */
  int nr_tasks;

  /* Is the data of this cell being used in a sub-cell? */
  int hold, ghold;

  /* Spin lock for various uses. */
  lock_type lock, glock;

  /* ID of the previous owner, e.g. runner. */
  int owner;

  /* Momentum of particles in cell. */
  float mom[3], ang[3];

  /* Potential and kinetic energy of particles in this cell. */
  double epot, ekin;

  /* Number of particles updated in this cell. */
  int updated;

  /* Linking pointer for "memory management". */
  struct cell *next;

  /* ID of the node this cell lives on. */
  int nodeID;

  /* Bit mask of the proxies this cell is registered with. */
  unsigned long long int sendto;

  /* Pointer to this cell's packed representation. */
  struct pcell *pcell;
  int pcell_size;
  int tag;

  /* This cell's multipole. */
  struct multipole multipole;

} __attribute__((aligned(64)));

/* Function prototypes. */
void cell_split(struct cell *c);
int cell_locktree(struct cell *c);
void cell_unlocktree(struct cell *c);
int cell_glocktree(struct cell *c);
void cell_gunlocktree(struct cell *c);
int cell_pack(struct cell *c, struct pcell *pc);
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s);
int cell_getsize(struct cell *c);
int cell_link(struct cell *c, struct part *parts);

#endif /* SWIFT_CELL_H */
