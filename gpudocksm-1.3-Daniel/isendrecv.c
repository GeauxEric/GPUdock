// mpicc -std=c99 isendrecv.c
// mpirun -np 9 ./a.out

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ROWS 7
#define COLS 4


int
minimal_int (int a, int b)
{
  return a < b ? a : b;
}


int
main (int argc, char *argv[])
{
  int *a = (int *) malloc (sizeof (int) * ROWS * COLS);

  int rank, size;
  MPI_Init (NULL, NULL);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  const int receiver_rank = size - 1;

  MPI_Request req[size];
  MPI_Status st[size];
  const int mytag = 0;

  int rows_per_p = (int) ceilf ((float) ROWS / size);
  int n_active_p = (int) ceilf ((float) ROWS / rows_per_p);
  int begin_row = rows_per_p * rank;
  int end_row = minimal_int (begin_row + rows_per_p, ROWS) - 1;

  // compute
  if (rank < n_active_p) {
    printf ("rank %d: row %d - %d\n", rank, begin_row, end_row);

    for (int i = begin_row; i <= end_row; ++i)
      for (int j = 0; j < COLS; ++j)
	a[COLS * i + j] = rank;
  }

  // reduce
  if (rank < n_active_p) {     // sender
    int msg_size = COLS * (end_row - begin_row + 1);
    MPI_Isend (&a[COLS * begin_row], msg_size, MPI_INT, receiver_rank, mytag, MPI_COMM_WORLD, &req[rank]);
    //printf ("%d sending\n", rank);
    MPI_Wait (&req[rank], &st[rank]);
    //printf ("%d sent\n", rank);
  }
  if (rank == receiver_rank) { // receiver
    for (int i = 0; i < n_active_p; ++i) {
      int begin_row_i = rows_per_p * i;
      int end_row_i = minimal_int (begin_row_i + rows_per_p, ROWS) - 1;
      int msg_size = COLS * (end_row_i - begin_row_i + 1);
      MPI_Irecv (&a[COLS * begin_row_i], msg_size, MPI_INT, i, mytag, MPI_COMM_WORLD, &req[i]);
      //printf ("%d receiving\n", rank);
    }
    MPI_Waitall (n_active_p, req, st);
    //printf ("%d all received\n", rank);
  }



  // print
  if (rank == receiver_rank)
    for (int i = 0; i < ROWS; ++i) {
      for (int j = 0; j < COLS; ++j)
	printf ("%2d\t", a[COLS * i + j]);
      printf ("\n");
    }


  MPI_Finalize ();
  free (a);

  return 0;
}
