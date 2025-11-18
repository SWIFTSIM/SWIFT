/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "proxy.h"

#ifdef WITH_MPI
/* Number of particle types to wait for after launching the proxies. We have
   parts, xparts, gparts, sparts, bparts and sinks to exchange, hence 6 types.
 */
#define MPI_REQUEST_NUMBER_PARTICLE_TYPES 6
#endif

/**
 * @brief Exchange straying particles with other nodes.
 *
 * @param e The #engine.
 * @param offset_parts The index in the parts array as of which the foreign
 *        parts reside (i.e. the current number of local #part).
 * @param ind_part The foreign #cell ID of each part.
 * @param Npart The number of stray parts, contains the number of parts received
 *        on return.
 * @param offset_gparts The index in the gparts array as of which the foreign
 *        parts reside (i.e. the current number of local #gpart).
 * @param ind_gpart The foreign #cell ID of each gpart.
 * @param Ngpart The number of stray gparts, contains the number of gparts
 *        received on return.
 * @param offset_sparts The index in the sparts array as of which the foreign
 *        parts reside (i.e. the current number of local #spart).
 * @param ind_spart The foreign #cell ID of each spart.
 * @param Nspart The number of stray sparts, contains the number of sparts
 *        received on return.
 * @param offset_bparts The index in the bparts array as of which the foreign
 *        parts reside (i.e. the current number of local #bpart).
 * @param ind_bpart The foreign #cell ID of each bpart.
 * @param Nbpart The number of stray bparts, contains the number of bparts
 *        received on return.
 * @param offset_sinks The index in the sinks array as of which the foreign
 *        parts reside (i.e. the current number of local #sink).
 * @param ind_sink The foreign #cell ID of each sink.
 * @param Nsink The number of stray sinks, contains the number of sinks
 *        received on return.
 *
 * Note that this function does not mess-up the linkage between parts and
 * gparts, i.e. the received particles have correct linkeage.
 */
void engine_exchange_strays(
    struct engine* e, const size_t offset_parts, const int* restrict ind_part,
    size_t* Npart, const size_t offset_gparts, const int* restrict ind_gpart,
    size_t* Ngpart, const size_t offset_sparts, const int* restrict ind_spart,
    size_t* Nspart, const size_t offset_bparts, const int* restrict ind_bpart,
    size_t* Nbpart, const size_t offset_sinks, const int* restrict ind_sink,
    size_t* Nsink) {

#ifdef WITH_MPI
  struct space* s = e->s;
  ticks tic = getticks();

  /* Re-set the proxies. */
  for (int k = 0; k < e->nr_proxies; k++) {
    e->proxies[k].nr_parts_out = 0;
    e->proxies[k].nr_gparts_out = 0;
    e->proxies[k].nr_sparts_out = 0;
    e->proxies[k].nr_bparts_out = 0;
    e->proxies[k].nr_sinks_out = 0;
  }

  /* Put the parts into the corresponding proxies. */
  for (size_t k = 0; k < *Npart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_part[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_part[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->parts[offset_parts + k].id,
          s->parts[offset_parts + k].x[0], s->parts[offset_parts + k].x[1],
          s->parts[offset_parts + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the part. */
    if (s->parts[offset_parts + k].gpart != NULL) {
      s->parts[offset_parts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_parts_out;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->parts[offset_parts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the part and xpart into the proxy. */
    proxy_parts_load(&e->proxies[pid], &s->parts[offset_parts + k],
                     &s->xparts[offset_parts + k], 1);

#ifdef WITH_CSDS
    if (e->policy & engine_policy_csds) {
      /* Log the particle when leaving a rank. */
      csds_log_part(e->csds, &s->parts[offset_parts + k],
                    &s->xparts[offset_parts + k], e, /* log_all_fields */ 1,
                    csds_flag_mpi_exit, node_id);
    }
#endif
  }

  /* Put the sparts into the corresponding proxies. */
  for (size_t k = 0; k < *Nspart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_spart[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_spart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->sparts[offset_sparts + k].id,
          s->sparts[offset_sparts + k].x[0], s->sparts[offset_sparts + k].x[1],
          s->sparts[offset_sparts + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the spart. */
    if (s->sparts[offset_sparts + k].gpart != NULL) {
      s->sparts[offset_sparts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_sparts_out;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->sparts[offset_sparts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the spart into the proxy */
    proxy_sparts_load(&e->proxies[pid], &s->sparts[offset_sparts + k], 1);

#ifdef WITH_CSDS
    if (e->policy & engine_policy_csds) {
      /* Log the particle when leaving a rank. */
      csds_log_spart(e->csds, &s->sparts[offset_sparts + k], e,
                     /* log_all_fields */ 1, csds_flag_mpi_exit, node_id);
    }
#endif
  }

  /* Put the bparts into the corresponding proxies. */
  for (size_t k = 0; k < *Nbpart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_bpart[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_bpart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->bparts[offset_bparts + k].id,
          s->bparts[offset_bparts + k].x[0], s->bparts[offset_bparts + k].x[1],
          s->bparts[offset_bparts + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the bpart. */
    if (s->bparts[offset_bparts + k].gpart != NULL) {
      s->bparts[offset_bparts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_bparts_out;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->bparts[offset_bparts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the bpart into the proxy */
    proxy_bparts_load(&e->proxies[pid], &s->bparts[offset_bparts + k], 1);

#ifdef WITH_CSDS
    if (e->policy & engine_policy_csds) {
      error("Not yet implemented.");
    }
#endif
  }

  /* Put the sinks into the corresponding proxies. */
  for (size_t k = 0; k < *Nsink; k++) {
    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_sink[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_sink[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->sinks[offset_sinks + k].id,
          s->sinks[offset_sinks + k].x[0], s->sinks[offset_sinks + k].x[1],
          s->sinks[offset_sinks + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the sink. */
    if (s->sinks[offset_sinks + k].gpart != NULL) {
      s->sinks[offset_sinks + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_sinks_out;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->sinks[offset_sinks + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the sink into the proxy */
    proxy_sinks_load(&e->proxies[pid], &s->sinks[offset_sinks + k], 1);

#ifdef WITH_CSDS
    if (e->policy & engine_policy_csds) {
      error("Not yet implemented.");
    }
#endif
  }

  /* Put the gparts into the corresponding proxies. */
  for (size_t k = 0; k < *Ngpart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_gpart[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_gpart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lli, x=[%e,%e,%e].",
          node_id, s->gparts[offset_gparts + k].id_or_neg_offset,
          s->gparts[offset_gparts + k].x[0], s->gparts[offset_gparts + k].x[1],
          s->gparts[offset_gparts + k].x[2]);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->gparts[offset_gparts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the gpart into the proxy */
    proxy_gparts_load(&e->proxies[pid], &s->gparts[offset_gparts + k], 1);

#ifdef WITH_CSDS
    /* Write only the dark matter particles */
    if ((e->policy & engine_policy_csds) &&
        s->gparts[offset_gparts + k].type == swift_type_dark_matter) {

      /* Log the particle when leaving a rank. */
      csds_log_gpart(e->csds, &s->gparts[offset_gparts + k], e,
                     /* log_all_fields */ 1, csds_flag_mpi_exit, node_id);
    }
#endif
  }

  /* Launch the proxies. */
  MPI_Request reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * engine_maxproxies];
  MPI_Request reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * engine_maxproxies];
  for (int k = 0; k < e->nr_proxies; k++) {
    proxy_parts_exchange_first(&e->proxies[k]);
    reqs_in[k] = e->proxies[k].req_parts_count_in;
    reqs_out[k] = e->proxies[k].req_parts_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (int k = 0; k < e->nr_proxies; k++) {
    int pid = MPI_UNDEFINED;
    if (MPI_Waitany(e->nr_proxies, reqs_in, &pid, MPI_STATUS_IGNORE) !=
            MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "request from proxy %i has arrived." , pid );
    proxy_parts_exchange_second(&e->proxies[pid]);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(e->nr_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Count the total number of incoming particles and make sure we have
     enough space to accommodate them. */
  int count_parts_in = 0;
  int count_gparts_in = 0;
  int count_sparts_in = 0;
  int count_bparts_in = 0;
  int count_sinks_in = 0;
  for (int k = 0; k < e->nr_proxies; k++) {
    count_parts_in += e->proxies[k].nr_parts_in;
    count_gparts_in += e->proxies[k].nr_gparts_in;
    count_sparts_in += e->proxies[k].nr_sparts_in;
    count_bparts_in += e->proxies[k].nr_bparts_in;
    count_sinks_in += e->proxies[k].nr_sinks_in;
  }
  if (e->verbose) {
    message(
        "sent out %zu/%zu/%zu/%zu/%zu parts/gparts/sparts/bparts/sinks, got "
        "%i/%i/%i/%i/%i "
        "back.",
        *Npart, *Ngpart, *Nspart, *Nbpart, *Nsink, count_parts_in,
        count_gparts_in, count_sparts_in, count_bparts_in, count_sinks_in);
  }

  /* Reallocate the particle arrays if necessary */
  if (offset_parts + count_parts_in > s->size_parts) {
    s->size_parts = (offset_parts + count_parts_in) * engine_parts_size_grow;
    struct part* parts_new = NULL;
    struct xpart* xparts_new = NULL;
    if (swift_memalign("parts", (void**)&parts_new, part_align,
                       sizeof(struct part) * s->size_parts) != 0 ||
        swift_memalign("xparts", (void**)&xparts_new, xpart_align,
                       sizeof(struct xpart) * s->size_parts) != 0)
      error("Failed to allocate new part data.");
    memcpy(parts_new, s->parts, sizeof(struct part) * offset_parts);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * offset_parts);
    swift_free("parts", s->parts);
    swift_free("xparts", s->xparts);
    s->parts = parts_new;
    s->xparts = xparts_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_parts; k++) {
      if (s->parts[k].gpart != NULL) {
        s->parts[k].gpart->id_or_neg_offset = -k;
      }
    }
  }

  if (offset_sparts + count_sparts_in > s->size_sparts) {
    s->size_sparts = (offset_sparts + count_sparts_in) * engine_parts_size_grow;
    struct spart* sparts_new = NULL;
    if (swift_memalign("sparts", (void**)&sparts_new, spart_align,
                       sizeof(struct spart) * s->size_sparts) != 0)
      error("Failed to allocate new spart data.");
    memcpy(sparts_new, s->sparts, sizeof(struct spart) * offset_sparts);
    swift_free("sparts", s->sparts);
    s->sparts = sparts_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_sparts; k++) {
      if (s->sparts[k].gpart != NULL) {
        s->sparts[k].gpart->id_or_neg_offset = -k;
      }
    }
  }

  if (offset_bparts + count_bparts_in > s->size_bparts) {
    s->size_bparts = (offset_bparts + count_bparts_in) * engine_parts_size_grow;
    struct bpart* bparts_new = NULL;
    if (swift_memalign("bparts", (void**)&bparts_new, bpart_align,
                       sizeof(struct bpart) * s->size_bparts) != 0)
      error("Failed to allocate new bpart data.");
    memcpy(bparts_new, s->bparts, sizeof(struct bpart) * offset_bparts);
    swift_free("bparts", s->bparts);
    s->bparts = bparts_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_bparts; k++) {
      if (s->bparts[k].gpart != NULL) {
        s->bparts[k].gpart->id_or_neg_offset = -k;
      }
    }
  }

  if (offset_sinks + count_sinks_in > s->size_sinks) {
    s->size_sinks = (offset_sinks + count_sinks_in) * engine_parts_size_grow;
    struct sink* sinks_new = NULL;
    if (swift_memalign("sinks", (void**)&sinks_new, sink_align,
                       sizeof(struct sink) * s->size_sinks) != 0)
      error("Failed to allocate new sink data.");
    memcpy(sinks_new, s->sinks, sizeof(struct sink) * offset_sinks);
    swift_free("sinks", s->sinks);
    s->sinks = sinks_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_sinks; k++) {
      if (s->sinks[k].gpart != NULL) {
        s->sinks[k].gpart->id_or_neg_offset = -k;
      }
    }
  }

  if (offset_gparts + count_gparts_in > s->size_gparts) {
    s->size_gparts = (offset_gparts + count_gparts_in) * engine_parts_size_grow;
    struct gpart* gparts_new = NULL;
    if (swift_memalign("gparts", (void**)&gparts_new, gpart_align,
                       sizeof(struct gpart) * s->size_gparts) != 0)
      error("Failed to allocate new gpart data.");
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * offset_gparts);
    swift_free("gparts", s->gparts);
    s->gparts = gparts_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_gparts; k++) {
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_stars) {
        s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_black_hole) {
        s->bparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_sink) {
        s->sinks[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }
    }
  }

  /* Collect the requests for the particle data from the proxies. */
  int nr_in = 0, nr_out = 0;
  for (int k = 0; k < e->nr_proxies; k++) {
    if (e->proxies[k].nr_parts_in > 0) {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k] =
          e->proxies[k].req_parts_in;
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 1] =
          e->proxies[k].req_xparts_in;
      nr_in += 2;
    } else {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k] =
          reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 1] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_gparts_in > 0) {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 2] =
          e->proxies[k].req_gparts_in;
      nr_in += 1;
    } else {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 2] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sparts_in > 0) {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 3] =
          e->proxies[k].req_sparts_in;
      nr_in += 1;
    } else {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 3] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_bparts_in > 0) {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 4] =
          e->proxies[k].req_bparts_in;
      nr_in += 1;
    } else {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 4] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sinks_in > 0) {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 5] =
          e->proxies[k].req_sinks_in;
      nr_in += 1;
    } else {
      reqs_in[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 5] = MPI_REQUEST_NULL;
    }

    if (e->proxies[k].nr_parts_out > 0) {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k] =
          e->proxies[k].req_parts_out;
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 1] =
          e->proxies[k].req_xparts_out;
      nr_out += 2;
    } else {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k] =
          reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 1] =
              MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_gparts_out > 0) {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 2] =
          e->proxies[k].req_gparts_out;
      nr_out += 1;
    } else {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 2] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sparts_out > 0) {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 3] =
          e->proxies[k].req_sparts_out;
      nr_out += 1;
    } else {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 3] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_bparts_out > 0) {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 4] =
          e->proxies[k].req_bparts_out;
      nr_out += 1;
    } else {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 4] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sinks_out > 0) {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 5] =
          e->proxies[k].req_sinks_out;
      nr_out += 1;
    } else {
      reqs_out[MPI_REQUEST_NUMBER_PARTICLE_TYPES * k + 5] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for each part array to come in and collect the new
     parts from the proxies. */
  int count_parts = 0, count_gparts = 0, count_sparts = 0, count_bparts = 0,
      count_sinks = 0;
  for (int k = 0; k < nr_in; k++) {
    int err, pid;
    if ((err = MPI_Waitany(MPI_REQUEST_NUMBER_PARTICLE_TYPES * e->nr_proxies,
                           reqs_in, &pid, MPI_STATUS_IGNORE)) != MPI_SUCCESS) {
      char buff[MPI_MAX_ERROR_STRING];
      int res;
      MPI_Error_string(err, buff, &res);
      error("MPI_Waitany failed (%s).", buff);
    }
    if (pid == MPI_UNDEFINED) break;
    // message( "request from proxy %i has arrived." , pid /
    // MPI_REQUEST_NUMBER_PARTICLE_TYPES );
    pid = MPI_REQUEST_NUMBER_PARTICLE_TYPES *
          (pid / MPI_REQUEST_NUMBER_PARTICLE_TYPES);

    /* If all the requests for a given proxy have arrived... */
    if (reqs_in[pid + 0] == MPI_REQUEST_NULL &&
        reqs_in[pid + 1] == MPI_REQUEST_NULL &&
        reqs_in[pid + 2] == MPI_REQUEST_NULL &&
        reqs_in[pid + 3] == MPI_REQUEST_NULL &&
        reqs_in[pid + 4] == MPI_REQUEST_NULL &&
        reqs_in[pid + 5] == MPI_REQUEST_NULL) {
      /* Copy the particle data to the part/xpart/gpart arrays. */
      struct proxy* prox = &e->proxies[pid / MPI_REQUEST_NUMBER_PARTICLE_TYPES];
      memcpy(&s->parts[offset_parts + count_parts], prox->parts_in,
             sizeof(struct part) * prox->nr_parts_in);
      memcpy(&s->xparts[offset_parts + count_parts], prox->xparts_in,
             sizeof(struct xpart) * prox->nr_parts_in);
      memcpy(&s->gparts[offset_gparts + count_gparts], prox->gparts_in,
             sizeof(struct gpart) * prox->nr_gparts_in);
      memcpy(&s->sparts[offset_sparts + count_sparts], prox->sparts_in,
             sizeof(struct spart) * prox->nr_sparts_in);
      memcpy(&s->bparts[offset_bparts + count_bparts], prox->bparts_in,
             sizeof(struct bpart) * prox->nr_bparts_in);
      memcpy(&s->sinks[offset_sinks + count_sinks], prox->sinks_in,
             sizeof(struct sink) * prox->nr_sinks_in);

#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds) {
        struct part* parts = &s->parts[offset_parts + count_parts];
        struct xpart* xparts = &s->xparts[offset_parts + count_parts];
        struct spart* sparts = &s->sparts[offset_sparts + count_sparts];
        struct gpart* gparts = &s->gparts[offset_gparts + count_gparts];

        /* Log the gas particles */
        csds_log_parts(e->csds, parts, xparts, prox->nr_parts_in, e,
                       /* log_all_fields */ 1, csds_flag_mpi_enter,
                       prox->nodeID);

        /* Log the stellar particles */
        csds_log_sparts(e->csds, sparts, prox->nr_sparts_in, e,
                        /* log_all_fields */ 1, csds_flag_mpi_enter,
                        prox->nodeID);

        /* Log the gparts */
        csds_log_gparts(e->csds, gparts, prox->nr_gparts_in, e,
                        /* log_all_fields */ 1, csds_flag_mpi_enter,
                        prox->nodeID);

        /* Log the bparts */
        if (prox->nr_bparts_in > 0) {
          error("TODO");
        }

        /* Log the sinks */
        if (prox->nr_sinks_in > 0) {
          /* Not implemented yet */
          error("TODO");
        }
      }
#endif
      /* for (int k = offset; k < offset + count; k++)
         message(
            "received particle %lli, x=[%.3e %.3e %.3e], h=%.3e, from node %i.",
            s->parts[k].id, s->parts[k].x[0], s->parts[k].x[1],
            s->parts[k].x[2], s->parts[k].h, p->nodeID); */

      /* Re-link the gparts. */
      for (int kk = 0; kk < prox->nr_gparts_in; kk++) {
        struct gpart* gp = &s->gparts[offset_gparts + count_gparts + kk];

        if (gp->type == swift_type_gas) {
          struct part* p =
              &s->parts[offset_parts + count_parts - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->parts - p;
          p->gpart = gp;
        } else if (gp->type == swift_type_stars) {
          struct spart* sp =
              &s->sparts[offset_sparts + count_sparts - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->sparts - sp;
          sp->gpart = gp;
        } else if (gp->type == swift_type_black_hole) {
          struct bpart* bp =
              &s->bparts[offset_bparts + count_bparts - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->bparts - bp;
          bp->gpart = gp;
        } else if (gp->type == swift_type_sink) {
          struct sink* sink =
              &s->sinks[offset_sinks + count_sinks - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->sinks - sink;
          sink->gpart = gp;
        }
      }

      /* Advance the counters. */
      count_parts += prox->nr_parts_in;
      count_gparts += prox->nr_gparts_in;
      count_sparts += prox->nr_sparts_in;
      count_bparts += prox->nr_bparts_in;
      count_sinks += prox->nr_sinks_in;
    }
  }

  /* Wait for all the sends to have finished too. */
  if (nr_out > 0)
    if (MPI_Waitall(MPI_REQUEST_NUMBER_PARTICLE_TYPES * e->nr_proxies, reqs_out,
                    MPI_STATUSES_IGNORE) != MPI_SUCCESS)
      error("MPI_Waitall on sends failed.");

  /* Free the proxy memory */
  for (int k = 0; k < e->nr_proxies; k++) {
    proxy_free_particle_buffers(&e->proxies[k]);
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Return the number of harvested parts. */
  *Npart = count_parts;
  *Ngpart = count_gparts;
  *Nspart = count_sparts;
  *Nbpart = count_bparts;
  *Nsink = count_sinks;

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
