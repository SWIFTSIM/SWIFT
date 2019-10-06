/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_PROXY_H
#define SWIFT_PROXY_H

/* Includes. */
#include "cell.h"
#include "part.h"
#include "space.h"

/* Some constants. */
#define proxy_buffgrow 1.5
#define proxy_buffinit 100

/* Proxy tag arithmetic. */
#define proxy_tag_shift 8
#define proxy_tag_count 0
#define proxy_tag_parts 1
#define proxy_tag_xparts 2
#define proxy_tag_gparts 3
#define proxy_tag_sparts 4
#define proxy_tag_bparts 5
#define proxy_tag_cells 6

/**
 * @brief The different reasons a cell can be in a proxy
 */
enum proxy_cell_type {
  proxy_cell_type_none = 0,
  proxy_cell_type_hydro = (1 << 0),
  proxy_cell_type_gravity = (1 << 1),
};

/* Data structure for the proxy. */
struct proxy {

  /* ID of the node this proxy represents. */
  int mynodeID, nodeID;

  /* Incoming cells. */
  struct cell **cells_in;
  int *cells_in_type;
  struct pcell *pcells_in;
  int nr_cells_in, size_cells_in, size_pcells_in;

  /* Outgoing cells. */
  struct cell **cells_out;
  int *cells_out_type;
  struct pcell *pcells_out;
  int nr_cells_out, size_cells_out, size_pcells_out;

  /* The parts and xparts buffers for input and output. */
  struct part *parts_in, *parts_out;
  struct xpart *xparts_in, *xparts_out;
  struct gpart *gparts_in, *gparts_out;
  struct spart *sparts_in, *sparts_out;
  struct bpart *bparts_in, *bparts_out;
  int size_parts_in, size_parts_out;
  int nr_parts_in, nr_parts_out;
  int size_gparts_in, size_gparts_out;
  int nr_gparts_in, nr_gparts_out;
  int size_sparts_in, size_sparts_out;
  int nr_sparts_in, nr_sparts_out;
  int size_bparts_in, size_bparts_out;
  int nr_bparts_in, nr_bparts_out;

  /* Buffer to hold the incomming/outgoing particle counts. */
  int buff_out[4], buff_in[4];

/* MPI request handles. */
#ifdef WITH_MPI
  MPI_Request req_parts_count_out, req_parts_count_in;
  MPI_Request req_parts_out, req_parts_in;
  MPI_Request req_xparts_out, req_xparts_in;
  MPI_Request req_gparts_out, req_gparts_in;
  MPI_Request req_sparts_out, req_sparts_in;
  MPI_Request req_bparts_out, req_bparts_in;
  MPI_Request req_cells_count_out, req_cells_count_in;
  MPI_Request req_cells_out, req_cells_in;
#endif
};

/* Function prototypes. */
void proxy_init(struct proxy *p, int mynodeID, int nodeID);
void proxy_clean(struct proxy *p);
void proxy_parts_load(struct proxy *p, const struct part *parts,
                      const struct xpart *xparts, int N);
void proxy_gparts_load(struct proxy *p, const struct gpart *gparts, int N);
void proxy_sparts_load(struct proxy *p, const struct spart *sparts, int N);
void proxy_bparts_load(struct proxy *p, const struct bpart *bparts, int N);
void proxy_parts_exchange_first(struct proxy *p);
void proxy_parts_exchange_second(struct proxy *p);
void proxy_addcell_in(struct proxy *p, struct cell *c, int type);
void proxy_addcell_out(struct proxy *p, struct cell *c, int type);
void proxy_cells_exchange(struct proxy *proxies, int num_proxies,
                          struct space *s, int with_gravity);
void proxy_tags_exchange(struct proxy *proxies, int num_proxies,
                         struct space *s);
void proxy_create_mpi_type(void);
void proxy_free_mpi_type(void);

#endif /* SWIFT_PROXY_H */
