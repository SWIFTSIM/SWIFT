/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Peter W. Draper
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

/**
 * @brief simple C wrapper for the C++ infinity library, only provides an
 * interface to meta capabilities we use. Note still exposed to C++ linkage
 * from the infinity library so we must use a C++ compiler.
 */
/* Config parameters. */
//#include "../config.h"
#define HAVE_INFINITY

/* Standard includes. */
#include <arpa/inet.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Infinity C++ headers. */
#ifdef HAVE_INFINITY
#include <infinity/core/Context.h>
#include <infinity/memory/Buffer.h>
#include <infinity/memory/RegionToken.h>
#include <infinity/queues/QueuePair.h>
#include <infinity/queues/QueuePairFactory.h>
#include <infinity/requests/RequestToken.h>
#endif

/* Local defines. */
#include "infinity_wrapper.h"

/* Local includes. */
#include "error.h"

/* Size of a block of memory. MESSAGE_SIZE needs to be a multiple of this as
 * as we need to align in memory. */
#define BLOCKTYPE size_t
#define MPI_BLOCKTYPE MPI_AINT
static const int BYTESINBLOCK = sizeof(BLOCKTYPE);

/* Flags for controlling access. High end of size_t. */
static size_t UNLOCKED = (((size_t)2 << 63) - 1);

/* Struct of QPs and associated data. */
struct qps_data {
  int nr_qps;
  infinity::core::Context *context;
  infinity::queues::QueuePairFactory *factory;
  infinity::queues::QueuePair **qps;
  infinity::memory::Buffer **receive_buffers;
  infinity::memory::RegionToken **remote_buffers;
  infinity::memory::Buffer **readwrite_buffers;
  infinity::memory::RegionToken **token_buffers;
};

/**
 * @brief Find an IP address for the given hostname.
 *
 * @param hostname the hostname
 *
 * @result the IP address, note copy away to keep.
 */

static char *toipaddr(char *hostname) {

  struct hostent *hostent = gethostbyname(hostname);
  if (hostent == NULL) {
    error("Failed to convert hostname '%s' to an IP address", hostname);
  }
  struct in_addr **addr_list = (struct in_addr **)hostent->h_addr_list;
  return inet_ntoa(*addr_list[0]);
}

/**
 * @brief Create a QPs to connect a group of clients to a group of servers.
 *
 * Requires that infinity_create_servers is also running, otherwise we
 * block waiting for the connections.
 *
 * @param servers a #mpi_servers struct with the server details.
 * @param nr_servers the number of servers expected to connect.
 * @param myrank the MPI rank of this process.
 * @param verbose if 1 then report the connections made.
 *
 * @return handle for the QPs and related data.
 */
void *infinity_connect_clients(struct mpi_servers *servers, int nr_servers,
                               int myrank, int verbose) {
#ifdef HAVE_INFINITY

  /* Struct to hold all the persistent data. */
  struct qps_data *cqps = (struct qps_data *)calloc(1, sizeof(struct qps_data));

  /* Need a factory to create QPs. */
  cqps->context = new infinity::core::Context();
  cqps->factory = new infinity::queues::QueuePairFactory(cqps->context);

  /* Create the QPs connecting to all the other ranks. */
  cqps->qps = (infinity::queues::QueuePair **)
    calloc(nr_servers, sizeof(infinity::queues::QueuePair *));
  cqps->nr_qps = nr_servers;

  /* Space for the pointers to the remote memory. */
  cqps->remote_buffers = (infinity::memory::RegionToken **)
    calloc(nr_servers, sizeof(infinity::memory::RegionToken *));

  /* We need to listen for messages from the other rank servers that we can
   * connect to them as they need to be up first. */
  int buf[1];
  MPI_Request reqs[nr_servers];
  for (int k = 0; k < nr_servers; k++) {
    if (k != myrank) {
      MPI_Irecv(buf, 1, MPI_INT, k, k, MPI_COMM_WORLD, &reqs[k]);
    } else {
      reqs[myrank] = MPI_REQUEST_NULL;
    }
  }

  /* Now we poll for any servers that are ready to connect. */
  int index;
  MPI_Status status;
  while (1) {
    MPI_Waitany(nr_servers, reqs, &index, &status);

    /* All done when all requests have completed. */
    if (index == MPI_UNDEFINED) break;

    /*  Got one, so connect. */
    char *ip = &servers->ip[index * infinity_max_server_ip];
    if (verbose)
      message("%d waiting for connection to remote server %s %d on %d", myrank,
              ip, index, BASE_PORT + myrank);
    cqps->qps[index] = cqps->factory->connectToRemoteHost(ip, BASE_PORT + myrank);
    if (verbose)
      message("%d connected to remote server %s %d on %d", myrank, ip, index,
              BASE_PORT + myrank);

    /* Remote buffer access */
    cqps->remote_buffers[index] =
      (infinity::memory::RegionToken *)cqps->qps[index]->getUserData();
  }

  /* Result is opaque. */
  return (void *)cqps;

#else
  return NULL;
#endif
}

/**
 * @brief Send a buffer to a server listening on a QP.
 *
 * @param qphandle the handle from infinity_connect_clients.
 * @param index index of the server to send to.
 * @param buffer the buffer to send, should be block aligned.
 * @param size the size of the buffer in bytes.
 * @param offset the offset into the remote buffer, note in bytes not blocks.
 */
void infinity_send_data(void *qphandle, int index, void *buffer, size_t size,
                        size_t offset) {
#ifdef HAVE_INFINITY
  struct qps_data *cqps = (struct qps_data *)qphandle;

  /* Need to assign to a buffer to register memory. XXX make this as big as
   * necessary per server and reuse. */
  auto *sendBuffer =
    new infinity::memory::Buffer(cqps->context, buffer, size);

  /* And send. */
  infinity::requests::RequestToken requestToken(cqps->context);
  cqps->qps[index]->write(sendBuffer,
                          0,                           // localOffset
                          cqps->remote_buffers[index], // destination
                          offset,                      // remoteOffset
                          size,                        // sizeInBytes
                          infinity::queues::OperationFlags(),
                          &requestToken);
  requestToken.waitUntilCompleted();
  requestToken.reset();

  /* Now we update the unlock field. */
  ((BLOCKTYPE *)sendBuffer->getData())[0] = UNLOCKED;
  cqps->qps[index]->write(sendBuffer,
                          0,                           // localOffset
                          cqps->remote_buffers[index], // destination
                          offset,                      // remoteOffset
                          BYTESINBLOCK,                // sizeInBytes
                          infinity::queues::OperationFlags(),
                          &requestToken);
  requestToken.waitUntilCompleted();  // Since we reuse the sendBuffer.

  delete sendBuffer;

#endif
  return;
}

/* @brief Free the resource associated with handle.
 *
 * @param qphandle the handle from infinity_connect_clients.
 */
void infinity_clients_free(void *qphandle) {

#ifdef HAVE_INFINITY
  struct qps_data *cqps = (struct qps_data *)qphandle;
  for (int k = 0; k < cqps->nr_qps; k++) delete cqps->qps[k];
  free(cqps->qps);
  delete cqps->factory;
  delete cqps->context;
  free(cqps->receive_buffers);
  free(cqps->remote_buffers);
  if (cqps->readwrite_buffers != NULL) {
    for (int k = 0; k < cqps->nr_qps; k++) delete cqps->readwrite_buffers[k];
    free(cqps->readwrite_buffers);
  }
  if (cqps->token_buffers != NULL) {
    for (int k = 0; k < cqps->nr_qps; k++) delete cqps->token_buffers[k];
    free(cqps->token_buffers);
  }
  free(cqps);
#endif
  return;
}

/**
 * @brief Create QPs for server to receive data from our clients.
 *
 * Requires that infinity_connect_clients is also ran, otherwise we
 * block waiting for the connections.
 *
 * @param servers a #mpi_servers struct with the server details.
 * @param nr_servers the number of servers we will create.
 * @param sizes the sizes, in bytes, of the various windows needed to receive
 *              all the remote data from a client. Array size of nr_servers.
 * @param myrank the MPI rank of this process.
 * @param verbose if 1 then report the connections made.
 *
 * @return handle for the QPs and related data.
 */
void *infinity_create_servers(struct mpi_servers *servers, int nr_servers,
                              size_t *sizes, int myrank, int verbose) {

#ifdef HAVE_INFINITY
  /* Struct to hold all the persistent data. */
  struct qps_data *cqps = (struct qps_data *)calloc(1, sizeof(struct qps_data));

  /* Need a factory to create QPs. */
  cqps->context = new infinity::core::Context();
  cqps->factory = new infinity::queues::QueuePairFactory(cqps->context);

  /* Create the QPs connecting to all the other ranks. */
  cqps->qps = (infinity::queues::QueuePair **)
    calloc(nr_servers, sizeof(infinity::queues::QueuePair *));
  cqps->nr_qps = nr_servers;

  /* Create buffers to receive all the remote data. */
  cqps->readwrite_buffers = (infinity::memory::Buffer **)
    calloc(nr_servers, sizeof(infinity::memory::Buffer *));
  cqps->token_buffers = (infinity::memory::RegionToken **)
    calloc(nr_servers, sizeof(infinity::memory::RegionToken *));

  for (int k = 0; k < nr_servers; k++) {
    if (sizes[k] > 0) {
      cqps->readwrite_buffers[k] =
        new infinity::memory::Buffer(cqps->context, sizes[k]);
    } else {
      /* Dummy: not expecting any data, but setup anyway. */
      cqps->readwrite_buffers[k] =
        new infinity::memory::Buffer(cqps->context, BYTESINBLOCK);
    }
    cqps->token_buffers[k] = cqps->readwrite_buffers[k]->createRegionToken();
  }

  /* Do the port binding for each other rank. */
  int buf[1];
  MPI_Request req;
  for (int k = 0; k < nr_servers; k++) {
    if (k != myrank) {
      if (verbose)
        message("%d binding to %d on port %d", myrank, k, BASE_PORT + k);
      cqps->factory->bindToPort(BASE_PORT + k);

      /* Send message this port is about to block for a connection. */
      if (verbose) message("Blocking for first message on %d", BASE_PORT + k);
      MPI_Isend(buf, 1, MPI_INT, k, myrank, MPI_COMM_WORLD, &req);
      cqps->qps[k] = cqps->factory->acceptIncomingConnection
        (cqps->token_buffers[k], sizeof(infinity::memory::RegionToken));
      if (verbose)
        message("Accepting incoming connections on %d", BASE_PORT + k);
    }
  }

  return (void *)cqps;
#else
  return NULL;
#endif
}

/**
 * @brief Check if data is ready, that is has arrived.
 *
 * @param qphandle the handle from infinity_create_servers.
 * @param index index of the client we are checking.
 * @param offset the offset of this data in the RDMA buffer,
 *               note in blocks not bytes.
 *
 * @result pointer to the start of the data, otherwise NULL.
 */
void *infinity_check_ready(void *qphandle, int index, size_t offset) {

  void *result = NULL;
#ifdef HAVE_INFINITY
  struct qps_data *cqps = (struct qps_data *)qphandle;

  /* Get the data location. */
  BLOCKTYPE *dataptr = &((BLOCKTYPE *)cqps->readwrite_buffers[index]->getData())[offset];

  /* Check if this has been unlocked. */
  BLOCKTYPE volatile lock = dataptr[0];
  if (lock == UNLOCKED) result = (void *)dataptr;

#endif
  return result;
}
