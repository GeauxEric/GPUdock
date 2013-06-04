#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "threads.h"
#include "mpiprocess.h"


void
mpiprocess (float beta_low, float beta_high, char *mydir)
{
  printf ("mpi_wrapper: hello\n");

  int rank, size;

  MPI_Init (NULL, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  printf ("hello from %02d/%02d\n", rank, size);
  thread_parent (beta_low, beta_high, mydir, rank);

  MPI_Finalize();
}


