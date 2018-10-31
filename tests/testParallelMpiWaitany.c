// Test MPI_Waitany in parallel.
//
// This test has each rank send each other rank an array of random ints with
// asynchronous sends and receives. The program spawns as many threads
// as it has incoming receives, which all call MPI_Waitany concurrently.
//
// Compile with:
//   mpicc testParallelMpiWaitany.c
//
// Run with:
//    mpirun -n 10 ./a.out

#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define BUFFER_SIZE 1000

struct thread_func_data {
  int num_ranks;
  MPI_Request *requests;
};

void *thread_func(void *data_in) {
  const struct thread_func_data *data = (struct thread_func_data *)data_in;

  // Wait on any of the requests.
  int index = MPI_UNDEFINED;
  MPI_Status status;
  const int res = MPI_Waitany(data->num_ranks, data->requests, &index, &status);

  // Make sure all went well.
  if (res != MPI_SUCCESS) {
    int len = 1024;
    char error_message[len];
    MPI_Error_string(res, error_message, &len);
    fprintf(stderr, "MPI_Waitany failed: %s\n", error_message);
  }
  if (index == MPI_UNDEFINED) {
    fprintf(stderr, "MPI_Waitany didn't return a valid index.\n");
  }

  return NULL;
}

int main(int argc, char **argv) {
  // Initialize the MPI environment.
  MPI_Init(NULL, NULL);

  // Get the number of processes.
  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  // Get the rank of the process.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Rank %02i is ready.\n", rank);

  // Create a data buffers, fill the outgoing one with random numbers.
  int *data_out = malloc(sizeof(int) * BUFFER_SIZE);
  for (int k = 0; k < BUFFER_SIZE; k++) {
    data_out[k] = rand();
  }
  int *data_in = malloc(sizeof(int) * BUFFER_SIZE);

  // Create the request objects.
  MPI_Request requests_in[num_ranks];
  MPI_Request requests_out[num_ranks];

  // Asynchronously send the buffer to all other ranks.
  for (int rid = 0; rid < num_ranks; rid++) {
    if (rid != rank) {
      MPI_Isend(data_out, BUFFER_SIZE, MPI_INT, rid, /*tag=*/0, MPI_COMM_WORLD,
                &requests_out[rid]);
    } else {
      requests_out[rid] = MPI_REQUEST_NULL;
    }
  }

  // Launch the asynchronous receives.
  for (int rid = 0; rid < num_ranks; rid++) {
    if (rid != rank) {
      MPI_Irecv(data_in, BUFFER_SIZE, MPI_INT, rid, /*tag=*/0, MPI_COMM_WORLD,
                &requests_in[rid]);
    } else {
      requests_in[rid] = MPI_REQUEST_NULL;
    }
  }

  // Launch the threads that will call MPI_Waitany.
  struct thread_func_data thread_data = {num_ranks, requests_in};
  pthread_t threads[num_ranks - 1];
  for (int tid = 0; tid < num_ranks - 1; tid++) {
    pthread_create(&threads[tid], /*attr=*/NULL, thread_func, &thread_data);
  }

  // Wait for all the threads to be done.
  for (int tid = 0; tid < num_ranks - 1; tid++) {
    pthread_join(threads[tid], /*value_ptr=*/NULL);
  }

  // Just to be sure, make sure the sends are all done.
  MPI_Waitall(num_ranks, requests_out, MPI_STATUSES_IGNORE);

  // Clean up allocated memory.
  free(data_in);
  free(data_out);

  // Finalize the MPI environment.
  MPI_Finalize();
  printf("Rank %02i is done.\n", rank);
}
