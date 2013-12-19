#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <cstdio>

#include "debug.h"

using namespace std;

#define     MP3	 1
#define     TUNE 0

int main(int argc, char *argv[])
{
	int ROWS = 16;
	int COLUMNS = 16;
	int *a;
	int total_confs = ROWS;


	// to make sure the matrix elements are stored in a consequential mem
	// since the MPI_Send and MPI_Recv are exchanging consequential mem buffer
	a = new int[ROWS * COLUMNS];

#if MP3
	int rank;
	int size;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	int half_total = total_confs / 2;
	int half_rows_per_process = (int)ceil((float)half_total / size);
	int begin_row = half_rows_per_process * rank;
	int end_row = min(begin_row + half_rows_per_process - 1, half_total - 1);

	for (int i = begin_row; i <= end_row; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			a[COLUMNS * i + j] = i * 100 + j;
		}
	}

#if TUNE
	DEBUG_1_("rank:	", rank);
	DEBUG_1_("begin_row: ", begin_row);
	DEBUG_1_("end_row: ", end_row);
	DEBUG_2_("[][]", a[COLUMNS * begin_row]);
#endif

	//============================================================================
	//reduce
	MPI_Status send_status[size], recv_status[size];
	MPI_Request send_req[size], recv_req[size];
	const int tag_head = 1, tag_tail = 2;
	// sender
	if (rank != 0) {
		int msg_size = COLUMNS * (end_row - begin_row + 1);
		MPI_Isend(&a[COLUMNS * begin_row + 0], msg_size, MPI_FLOAT, 0, tag_head, MPI_COMM_WORLD, &send_req[rank]);

		//DEBUG_2_("sending...", rank);
		MPI_Wait(&send_req[rank], &send_status[rank]);
		//DEBUG_2_("sent...", rank);

	}
	// receiver
	else {
		for (int sender = 1; sender < size; ++sender) {
			int begin_row_sender = half_rows_per_process * sender;
			int end_row_sender = min(begin_row_sender + half_rows_per_process - 1, half_total - 1);
			int msg_size = COLUMNS * (end_row_sender - begin_row_sender + 1);
			MPI_Irecv(&a[COLUMNS * begin_row_sender + 0], msg_size, MPI_FLOAT,
				  sender, tag_head, MPI_COMM_WORLD, &recv_req[sender]);
		}

                              DEBUG_2_("recving...", 888);
                            MPI_Waitall(size-1, &recv_req[1], &recv_status[1]);
                              DEBUG_2_("recved...", 888);

		for (int sender = 1; sender < size; ++sender) {
			//DEBUG_2_("recving...", sender);
			// MPI_Wait(&recv_req[sender], &recv_status[sender]);
			//DEBUG_2_("received...", sender);
		}
	}

#if 1
	if (rank == 0) {
		// DEBUG_3_("[]", a[end_row][0]);
		for (int j = 0; j < ROWS; j++) {
			for (int i = 0; i < COLUMNS; i++) {
				printf("%04d\t", a[COLUMNS * j + i]);
				// DEBUG_4_("[]", a[end_row+1][i]);
			}
			cout << endl;
		}

	}
#endif

#if MP3
	MPI_Finalize();
#endif

	delete[]a;

	return 0;
}
