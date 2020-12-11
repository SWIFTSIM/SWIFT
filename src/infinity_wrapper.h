/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef INFINITY_WRAPPER_H
#define INFINITY_WRAPPER_H

/* Config parameters. */
#include "../config.h"

/* Base port no. Ranks use +rank. */
extern int infinity_base_port;

/* Maximum length of formatted server IP address. */
#define infinity_max_server_ip 24

/* Struct of MPI server ip addresses as formatted strings.*/
struct mpi_servers {
  char *ip;
};

#ifdef __cplusplus
extern "C" {
#endif
void infinity_open_communications(int nr_servers, size_t *sizes,
                                  void **recv_handle, void **send_handle);
void infinity_send_data(void *qphandle, int index, void *buffer, size_t size,
                        size_t offset);
void infinity_free_handle(void *qphandle);
void *infinity_check_ready(void *qphandle, int index, size_t offset);
#ifdef __cplusplus
}
#endif

void *infinity_connect_clients(struct mpi_servers *servers, int nr_servers,
                               int myrank, int verbose);
void *infinity_create_servers(struct mpi_servers *servers, int nr_servers,
                              size_t *sizes, int myrank, int verbose);

#endif /* INFINITY_WRAPPER_H */
