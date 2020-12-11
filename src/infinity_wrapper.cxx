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
#include "../config.h"

/* Standard includes. */
#include <arpa/inet.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <netdb.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Infinity C++ headers. */
#if defined(HAVE_INFINITY) && defined(WITH_MPI)
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
extern "C" {
#include "clocks.h"
#include "error.h"
}

/* The base port. */
int infinity_base_port = 27771;

/* Size of a block of memory. MESSAGE_SIZE needs to be a multiple of this as
 * as we need to align in memory. */
#if defined(HAVE_INFINITY) && defined(WITH_MPI)
#define BLOCKTYPE size_t
#define MPI_BLOCKTYPE MPI_AINT
static const int BYTESINBLOCK = sizeof(BLOCKTYPE);

/* Flags for controlling access. High end of size_t. */
static size_t UNLOCKED = (((size_t)2 << 63) - 1);
#endif

/* Struct of QPs and associated data. */
struct qps_data {
  int nr_qps;
#if defined(HAVE_INFINITY) && defined(WITH_MPI)
  infinity::core::Context *context;
  infinity::queues::QueuePairFactory *factory;
  infinity::queues::QueuePair **qps;
  infinity::memory::Buffer **receive_buffers;
  infinity::memory::RegionToken **remote_buffers;
  infinity::memory::Buffer **readwrite_buffers;
  infinity::memory::RegionToken **token_buffers;
#endif
};

struct server_data {
  struct mpi_servers *servers;
  int nr_servers;
  size_t *sizes;
  int myrank;
  int verbose;
  void *handle;
};

/**
 * @brief Find an IP address for the given hostname.
 *
 * @param hostname the hostname
 *
 * @result the IP address, note copy away to keep.
 */
static char *infinity_toipaddr(char *hostname) {

  struct hostent *hostent = gethostbyname(hostname);
  if (hostent == NULL) {
    error("Failed to convert hostname '%s' to an IP address", hostname);
  }
  struct in_addr **addr_list = (struct in_addr **)hostent->h_addr_list;
  return inet_ntoa(*addr_list[0]);
}

/**
 * @brief Send side thread to connect RDMA clients.
 */
static void *infinity_send_thread(void *arg) {

  /* Connect QPs to the remote servers. */
  struct server_data *data = (struct server_data *)arg;
  data->handle = infinity_connect_clients(data->servers, 
                                          data->nr_servers,
                                          data->myrank,
                                          data->verbose);
}

/**
 * @brief recv thread, listens for remote sends from another rank.
 */
static void *infinity_recv_thread(void *arg) {

  struct server_data *data = (struct server_data *)arg;
  data->handle = infinity_create_servers(data->servers, 
                                         data->nr_servers,
                                         data->sizes,
                                         data->myrank, 
                                         data->verbose);
}

/**
 * @brief Open up for RDMA communications between all the servers and
 * clients.
 *
 * @param nr_servers number of MPI ranks.
 * @param sizes the size needed to receive messages.
 * @param recv_handle handle for the recv QPs.
 * @param send_handle handle for the send QPs.
 */
void infinity_open_communications(int nr_servers, size_t *sizes,
                                  void **recv_handle, void **send_handle) {

#if defined(HAVE_INFINITY) && defined(WITH_MPI)

  /* Increment base port for next time. XXX to avoid collisions, must do
   * better. */
  infinity_base_port += nr_servers;

  /* Get the IPs of all the ranks. */
  char name[MPI_MAX_PROCESSOR_NAME];
  int namelen = 0;
  MPI_Get_processor_name(name, &namelen);
  char ip[infinity_max_server_ip];
  strncpy(ip, infinity_toipaddr(name), infinity_max_server_ip);

  /* And distribute, so we all know everyone's IPs. */
  struct mpi_servers servers;
  servers.ip = (char *)malloc(sizeof(char) * nr_servers * infinity_max_server_ip);
  MPI_Allgather(ip, infinity_max_server_ip, MPI_BYTE, servers.ip,
                infinity_max_server_ip, MPI_BYTE, MPI_COMM_WORLD);

  if (engine_rank == 0) {
    message("RDMA servers will listen on:");
    for (int j = 0; j < nr_servers; j++) {
      for (int k = 0; k < nr_servers; k++) {
        if (k != j) {
          message("  %d: %s on port %d", j,
                  &servers.ip[j * infinity_max_server_ip], infinity_base_port + k);
        }
      }
    }
  }

  /* Now we need to set up the RDMA connections, client and servers block when
   * this starts up, so we need a couple of threads to do the work. */
  struct server_data recv_data;
  recv_data.servers = &servers;
  recv_data.nr_servers = nr_servers;
  recv_data.sizes = sizes;
  recv_data.myrank = engine_rank;
  recv_data.verbose = 1;
  recv_data.handle = NULL;

  pthread_t recvthread;
  if (pthread_create(&recvthread, NULL, &infinity_recv_thread, &recv_data) != 0)
    error("Failed to create recv thread to start RDMA servers.");

  struct server_data send_data;
  send_data.servers = &servers;
  send_data.nr_servers = nr_servers;
  send_data.myrank = engine_rank;
  send_data.verbose = 1;
  send_data.handle = NULL;

  pthread_t sendthread;
  if (pthread_create(&sendthread, NULL, &infinity_send_thread, &send_data) != 0)
    error("Failed to create send thread to connect RDMA clients.");

  /* Wait until the threads complete. */
  pthread_join(sendthread, NULL);
  pthread_join(recvthread, NULL);

  /* Pass out handles for accessing the QPs. */
  *recv_handle = recv_data.handle;
  *send_handle = send_data.handle;

#endif
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
#if defined(HAVE_INFINITY) && defined(WITH_MPI)

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
              ip, index, infinity_base_port + myrank);
    cqps->qps[index] = cqps->factory->connectToRemoteHost(ip, infinity_base_port + myrank);
    if (verbose)
      message("%d connected to remote server %s %d on %d", myrank, ip, index,
              infinity_base_port + myrank);

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
#if defined(HAVE_INFINITY) && defined(WITH_MPI)
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
void infinity_free_handle(void *qphandle) {

#if defined(HAVE_INFINITY) && defined(WITH_MPI)
  struct qps_data *cqps = (struct qps_data *)qphandle;
  if (cqps == NULL) return;
  if (cqps->qps != NULL) {
    for (int k = 0; k < cqps->nr_qps; k++) delete cqps->qps[k];
    free(cqps->qps);
  }
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
  delete cqps->factory;
  delete cqps->context;
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

#if defined(HAVE_INFINITY) && defined(WITH_MPI)
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
        message("%d binding to %d on port %d", myrank, k, infinity_base_port + k);
      cqps->factory->bindToPort(infinity_base_port + k);

      /* Send message this port is about to block for a connection. */
      if (verbose) message("Blocking for first message on %d", infinity_base_port + k);
      MPI_Isend(buf, 1, MPI_INT, k, myrank, MPI_COMM_WORLD, &req);
      cqps->qps[k] = cqps->factory->acceptIncomingConnection
        (cqps->token_buffers[k], sizeof(infinity::memory::RegionToken));
      if (verbose)
        message("Accepting incoming connections on %d", infinity_base_port + k);
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
#if defined(HAVE_INFINITY) && defined(WITH_MPI)
  struct qps_data *cqps = (struct qps_data *)qphandle;

  /* Get the data location. */
  BLOCKTYPE *dataptr = &((BLOCKTYPE *)cqps->readwrite_buffers[index]->getData())[offset];

  /* Check if this has been unlocked. */
  BLOCKTYPE volatile lock = dataptr[0];
  if (lock == UNLOCKED) result = (void *)dataptr;

#endif
  return result;
}
