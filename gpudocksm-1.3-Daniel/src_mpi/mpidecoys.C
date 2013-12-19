#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <mpi.h>
#include "mpidecoys.h"
#include "../lib/debug.h"
#include "assert.h"

#define     MPI		1
#define     TUNE	1
#define     ORIGINAL	0


using namespace std;

void Decoy(Complex * d_complex, std::string iname, std::string oname)
{
	ifstream infile;
	ofstream outfile;
	outfile.open(oname.c_str());

	vector < string > data;	// contains the input file content
	vector < float >mcc;

	string inputline;
	long int number_of_lines = 0;
	infile.open(iname.c_str());

	if (infile.is_open()) {
		while (!infile.eof()) {
			getline(infile, inputline);
			data.push_back(inputline);
			number_of_lines++;
		}
	}
	else{
	    cout << iname << " not exists" << endl;
	    exit(EXIT_FAILURE);
	}

	long int atoms, newpens, newlens;
	infile.clear();
	infile.seekg(0, ios::beg);
	atoms = d_complex->getLigandAtomsTotal();

	long int total_confs = ((number_of_lines - 1) / (atoms + 3));

	// -----------------------------
	// to tune the program, reduced total_confs
	// -----------------------------
	// total_confs = 3000;
	
	// allocate the 1D mcc matrx in sequential mem
	// long int is used given the large size of the matrix
	long int matrix_size = total_confs;
	long int row_size = matrix_size;
	long int col_size = matrix_size;
	double *mcc_matrix = new double[matrix_size * matrix_size];
	assert(mcc_matrix != NULL);

	double newcoords[MAXLIG][3];
	long int conf = 0;
	long int n1;			// index for fetching new native parameters
	double elapsed_time;

#if MPI
	int rank;
	int size;
	MPI_Init(NULL, NULL);

	MPI_Barrier ( MPI_COMM_WORLD );
	elapsed_time = - MPI_Wtime();

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	//==========================================================================
	// if odd number, last row computed individually
	if (total_confs % 2 == 1) {
		// compute the last conf by rank 0
		if (rank == 0) {
			mcc_matrix[matrix_size * matrix_size - 1] = 1.f;
		}
		row_size--;
	}
	//==========================================================================
	// dispatch the work and compute
	long int half_rows = row_size / 2;
	long int half_rows_per_process = (long int)ceil((float)half_rows / size);
	long int begin_row = min(half_rows_per_process * rank, half_rows - 1);
	long int end_row = min(begin_row + half_rows_per_process - 1, half_rows - 1);

#if 0
	DEBUG_2_("half_rows_per_process: ", half_rows_per_process);
	DEBUG_3_("rank: ", rank);
	DEBUG_3_("begin_row: ", begin_row);
	DEBUG_3_("end_row: ", end_row);
	DEBUG_3_("tail_begin: ", row_size - 1 - end_row);
	DEBUG_3_("tail_end: ", row_size - 1 - begin_row);
#endif

	//==========================================================================
	//compute the heads of the upper triangular matrix
	for (conf = begin_row; conf <= end_row; ++conf) {
/*{{{*/
		n1 = ((atoms + 3) * conf);

		newpens = atoi((data[n1]).c_str());
		newlens = atoi((data[n1 + 1]).c_str());

#if TUNE
		d_complex->setProteinEnsembleCurrent(newpens);
		d_complex->setLigandEnsembleCurrent(newlens);
#endif

		for (long int j = n1 + 2, k = 0; k < atoms; j++, k++) {
			stringstream vec1_entry((data[j]).c_str());
			vec1_entry >> newcoords[k][0] >> newcoords[k][1] >> newcoords[k][2];
		}

#if TUNE
		d_complex->setConfCoords(newcoords);	// set lig conf coords
		d_complex->clearConfusionMatrix();
		d_complex->createConfusionMatrix();
		// d_complex->clearContactMatrix();
		// d_complex->createContactMatrix();
#endif

		long int pens_p;
		long int lens_p;
		Complex complex_p = *d_complex;
		double coords_p[MAXLIG][3];

		long int conf2_p;
		long int n;		//used for vector entries

		for (conf2_p = conf; conf2_p < col_size; conf2_p++) {

#if 0
			DEBUG_3_("rank:", rank);
			DEBUG_3_("cycle: ", cycle);
			DEBUG_3_("conf2_p: ", conf2_p);
#endif

			//Calculate mcc for each conformation relative the the current "native"
			n = ((atoms + 3) * conf2_p);

			pens_p = atoi((data[n]).c_str());
			lens_p = atoi((data[n + 1]).c_str());

#if TUNE
			complex_p.setProteinEnsembleCurrent(pens_p);
			complex_p.setLigandEnsembleCurrent(lens_p);
#endif

			for (long int j = n + 2, k = 0; k < atoms; j++, k++) {
				stringstream vec_entry((data[j]).c_str());
				vec_entry >> coords_p[k][0] >> coords_p[k][1] >> coords_p[k][2];
			}

#if TUNE
			complex_p.setConfCoords(coords_p);
			complex_p.calculateMccMtrx();
#endif

			long int mcc_index = 0;

			if (conf == 0) {
				mcc_index = conf2_p;	
			} else {
				for (long int i = 1; i <= conf; i++) {
					mcc_index += total_confs - (i - 1);
				}
				mcc_index += conf2_p - conf;
			}


#if 0
			cout << mcc_index << "	" << complex_p.getEnergy(1) << endl;
			DEBUG_4_("conf: ", conf);
			DEBUG_4_("mcc_index: ", mcc_index);
			DEBUG_4_("cmcc: ", complex_p.getEnergy(1));
			DEBUG_4_("rank: ", rank);
			DEBUG_4_("conf2_p: ", conf2_p);
#endif

			mcc_matrix[col_size * conf + conf2_p] = complex_p.getAmcc();
		}

		// compute_a_row (row, mcc);
		// compute_a_row (total_confs - row - 1, mcc);/*}}}*/
	}

	//==========================================================================
	// compute the tails of the upper triangular matrix
	for (conf = row_size - 1 - end_row; conf <= row_size - 1 - begin_row; ++conf) {
/*{{{*/
		// cout << conf << endl;
		// DEBUG_4_("rank: ", rank);
		// DEBUG_4_("tails", conf);
		n1 = ((atoms + 3) * conf);

		newpens = atoi((data[n1]).c_str());
		newlens = atoi((data[n1 + 1]).c_str());

#if TUNE
		d_complex->setProteinEnsembleCurrent(newpens);
		d_complex->setLigandEnsembleCurrent(newlens);
#endif

		for (long int j = n1 + 2, k = 0; k < atoms; j++, k++) {
			stringstream vec1_entry((data[j]).c_str());
			vec1_entry >> newcoords[k][0] >> newcoords[k][1] >> newcoords[k][2];
		}

#if TUNE
		d_complex->setConfCoords(newcoords);	// set lig conf coords
		d_complex->clearConfusionMatrix();
		d_complex->createConfusionMatrix();
		// d_complex->clearContactMatrix();
		// d_complex->createContactMatrix();
#endif

		long int pens_p;
		long int lens_p;
		Complex complex_p = *d_complex;
		double coords_p[MAXLIG][3];

		long int conf2_p;
		long int n;		//used for vector entries

		for (conf2_p = conf; conf2_p < col_size; conf2_p++) {

#if 0
			DEBUG_3_("rank:", rank);
			DEBUG_3_("cycle: ", cycle);
			DEBUG_3_("conf2_p: ", conf2_p);
#endif
			//Calculate mcc for each conformation relative the the current "native"

			n = ((atoms + 3) * conf2_p);

			pens_p = atoi((data[n]).c_str());
			lens_p = atoi((data[n + 1]).c_str());

#if TUNE
			complex_p.setProteinEnsembleCurrent(pens_p);
			complex_p.setLigandEnsembleCurrent(lens_p);
#endif

			for (long int j = n + 2, k = 0; k < atoms; j++, k++) {
				stringstream vec_entry((data[j]).c_str());
				vec_entry >> coords_p[k][0] >> coords_p[k][1] >> coords_p[k][2];
			}

#if TUNE
			complex_p.setConfCoords(coords_p);
			complex_p.calculateMccMtrx();
#endif

			long int mcc_index = 0;

			if (conf == 0) {
				mcc_index = conf2_p;	// does rank start at 0 or 1?
			} else {
				for (long int i = 1; i <= conf; i++) {
					mcc_index += total_confs - (i - 1);
				}
				mcc_index += conf2_p - conf;
			}

			// cout <<  rank << conf << conf2_p << " mcc_index: " << mcc_index <<" cmcc " << complex_p.getEnergy(1) << endl;

#if 0
			cout << mcc_index << "	" << complex_p.getEnergy(1) << endl;
			DEBUG_4_("conf: ", conf);
			DEBUG_4_("mcc_index: ", mcc_index);
			DEBUG_4_("cmcc: ", complex_p.getEnergy(1));
			DEBUG_4_("rank: ", rank);
			DEBUG_4_("conf2_p: ", conf2_p);
#endif

			mcc_matrix[col_size * conf + conf2_p] = complex_p.getAmcc();

		}

		// compute_a_row (row, mcc);
		// compute_a_row (total_confs - row - 1, mcc);/*}}}*/
	}

	//==========================================================================
	// reduce
	MPI_Status send_status[size], recv_status[size];
	MPI_Request send_req[size], recv_req[size];
	const long int tag_head = 1, tag_tail = 2;
	// sender
	if (rank != 0) {
		long int msg_size = matrix_size * (end_row - begin_row + 1);	    // send a block of sequential mem/*{{{*/

		MPI_Isend(&mcc_matrix[matrix_size * begin_row], msg_size, MPI_DOUBLE,	// heads
			  0, tag_head, MPI_COMM_WORLD, &send_req[rank]);
		MPI_Isend(&mcc_matrix[matrix_size * (row_size - 1 - end_row)], msg_size, MPI_DOUBLE,	// tails
			  0, tag_tail, MPI_COMM_WORLD, &send_req[rank]);
		MPI_Wait(&send_req[rank], &send_status[rank]);/*}}}*/
	}
	// receiver
	else {
		for (long int sender = 1; sender < size; ++sender) {/*{{{*/
			long int begin_row_sender = min(half_rows_per_process * sender, half_rows - 1);
			long int end_row_sender = min(begin_row_sender + half_rows_per_process - 1, half_rows - 1);
			long int msg_size = matrix_size * (end_row_sender - begin_row_sender + 1);
			MPI_Irecv(&mcc_matrix[matrix_size * begin_row_sender], msg_size, MPI_DOUBLE,	// heads
				  sender, tag_head, MPI_COMM_WORLD, &recv_req[sender]);
			MPI_Irecv(&mcc_matrix[matrix_size * (row_size - 1 - end_row_sender)], msg_size, MPI_DOUBLE,	// tails
				  sender, tag_tail, MPI_COMM_WORLD, &recv_req[sender]);
		}
		MPI_Waitall(size - 1, &recv_req[1], &recv_status[1]);/*}}}*/
	}

#if  MPI
	if (rank == 0) {
		/* copy from upper triangular part to lower */
		for (long int i = 0; i < matrix_size; i++)
			for (long int j = 0; j < matrix_size; j++)
				mcc_matrix[matrix_size * j + i] = mcc_matrix[matrix_size * i + j];

		/* write to the output file */
		outfile << matrix_size << endl;
		for (long int i = 0; i < matrix_size; i++) {
			for (long int j = 0; j < matrix_size; j++) {
				outfile << mcc_matrix[matrix_size * i + j];
				outfile << " ";
			}
			outfile << endl;
		}
		outfile.close();

		elapsed_time += MPI_Wtime();
		DEBUG_1_("parallel running time: ", elapsed_time);
	}

	MPI_Finalize();
#endif

	infile.close();

	delete[]mcc_matrix;
}

std::ifstream & GotoLine(std::ifstream & file, int num)
{
	file.seekg(std::ios::beg);
	for (int i = 0; i < num - 1; ++i) {
		file.ignore(std::numeric_limits < std::streamsize >::max(), '\n');
	}
	return file;
}
