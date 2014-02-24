/*
==============================================================================================
     __________________ ____ ___              .___             __      _________   _____   
    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
           \/                                  \/           \/     \/        \/         \/ 

      GPU-accelerated hybrid-resolution ligand docking using Replica Exchange Monte Carlo

==============================================================================================
*/

/*
 *  this the codes to generate the native and decoys for machine learning
 *
 *  while the system being perturbed:
 *  1. the coords of the lig
 *  2. the conf number of the lig
 *  3. the conf number of the prt 
 *  4. value of each energy term 
 *  are recorded
 *
 *  <lig_id>MCcoefs.txt
 *  <lig_id>Decoys.txt 
 *  <lig_id>DecoyEner.txt
 *  as the output files
 */

#include "decoy_gen.h"
#include "debug.h"

using namespace std;

int main(int argc, char *argv[])
{
	Param *para = new Param[1];
//---------------------------------------------------------
// default settings of parameters 
	para->remc_steps = 122;
	para->remc_cycles = 122;
	para->remc_replicas = 6;	// total replica numbers
	para->mc_temp = 0.1;		// temperature in each replica
	para->step_t = 0.5;
	para->step_r = 5.0;
	para->dsphere = 6.0;
	para->boltzmann = 0.05;
	para->data_pt_total = 1000000000;

//---------------------------------------------------------
	string molid = "MOLID";

	srand(time(0));

	time_t t_start1, t_end1, t_start2, t_end2;

	time(&t_start1);

	cout << "------------------------------------------------------------" << endl
		<< "                         Decoys Generator" << endl
		<< "		    Generating decoys in the mcc space " << endl
		<< "		by running Replica Exchange Monte Carlo " << endl << endl
		<< "------------------------------------------------------------" << endl << endl;

	string protein_name;		//name of protein file/*{{{*/
	bool protein_opt = false;

	string compounds_name;
	bool compounds_opt = false;

	string lhm_name;
	bool lhm_opt = false;

	string output_name;
	bool output_opt = false;
	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-p") && i < argc) {
			protein_name = string(argv[i + 1]);
			protein_opt = true;
		}
		if (!strcmp(argv[i], "-l") && i < argc) {
			compounds_name = string(argv[i + 1]);
			compounds_opt = true;
		}
		if (!strcmp(argv[i], "-s") && i < argc) {
			lhm_name = string(argv[i + 1]);
			lhm_opt = true;
		}
		if (!strcmp(argv[i], "-dp") && i < argc) {
			para->data_pt_total = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-o") && i < argc) {
			output_name = string(argv[i + 1]);
			output_opt = true;
		}
		if (!strcmp(argv[i], "-i") && i < argc) {
			molid = string(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-nr") && i < argc) {
			para->remc_replicas = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-ns") && i < argc) {
			para->remc_steps = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-nc") && i < argc) {
			para->remc_cycles = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-t") && i < argc) {
			para->mc_temp = atof(argv[i + 1]);
		}
		if (!strcmp(argv[i], "-bz") && i < argc) {
			para->boltzmann = atof(argv[i + 1]);
		}
	}

	char *path1;
	path1 = getenv("GPUDOCKSMDAT_FF");
	if (path1 == NULL) {
		cout << "GPUDOCKSMDAT_FF is not set" << endl;
		exit(EXIT_FAILURE);
	}
	string data_path;
	data_path = getenv("GPUDOCKSMDAT_FF");
	if (!protein_opt) {
		cout << "Provide target protein structure" << endl;
		exit(EXIT_FAILURE);
	}

	if (!compounds_opt) {
		cout << "Provide compound library in SD format" << endl;
		exit(EXIT_FAILURE);
	}

	if (!lhm_opt) {
		cout << "Provide LHM potentials" << endl;
		exit(EXIT_FAILURE);
	}

	list < string > l1_sdf;
	list < string >::iterator i1_sdf;

	string line1;

	ifstream compounds_file(compounds_name.c_str());

	if (!compounds_file.is_open()) {
		cout << "Cannot open " << compounds_name << endl;
		exit(EXIT_FAILURE);
	}

	while (getline(compounds_file, line1))
		l1_sdf.push_back(line1);

	compounds_file.close();

	int gpu_n = 0;

	string lib1[MAXSDF];
	int lib2 = 0;

	Complex *gpu_complex[MAXLIB];

	for (i1_sdf = l1_sdf.begin(); i1_sdf != l1_sdf.end(); i1_sdf++) {	/*{{{ */
		lib1[lib2++] = (*i1_sdf);

		if ((*i1_sdf) == "$$$$") {
			if (lib2 > 10) {
				gpu_complex[gpu_n] = new Complex(0, 0);
				gpu_complex[gpu_n]->clearMoves();

				bool load1 = gpu_complex[gpu_n]->loadProtein(protein_name);

				if (load1) {
					cout << "Cannot read target protein structure" << endl;
					exit(EXIT_FAILURE);
				}

				bool load2 = gpu_complex[gpu_n]->loadParams(data_path);

				if (load2) {
					cout << "Cannot read parameter file" << endl;
					exit(EXIT_FAILURE);
				}

				bool load3 = gpu_complex[gpu_n]->loadLigand(lib1, lib2, molid);

				if (load3) {
					cout << "Cannot load ligand data" << endl;
					exit(EXIT_FAILURE);
				}

				bool load4 = gpu_complex[gpu_n]->loadLHM(lhm_name);

				if (load4) {
					cout << "Cannot read LHM potentials" << endl;
					exit(EXIT_FAILURE);
				}
				gpu_complex[gpu_n]->createTrialCoords();
				gpu_complex[gpu_n]->calculateEnergy();

				gpu_n++;
			}

			lib2 = 0;
		}
	}

	cout << "Ligand ID: " << gpu_complex[0]->getLigandID() << endl;
	cout << "MC steps: " << para->remc_steps << endl;
	cout << "remc cycles: " << para->remc_cycles << endl;
	cout << "remc replicas: " << para->remc_replicas << endl;
	cout << "MC temperatures: " << para->mc_temp << endl;
	cout << "boltzmann constant: " << para->boltzmann << endl << endl;

	OutStream *outfiles = new OutStream[1];
	outfiles->all_stop = 1;
	outfiles->low_stop = 1;
	outfiles->high_stop = 1;

	string lig_id = gpu_complex[0]->getLigandID();
	initOutStream(outfiles, lig_id);

	for (int in = 0; in < gpu_n; in++) {
		time(&t_start2);
		gpu_complex[in]->clearMoves();
		gpu_complex[in]->createTrialCoords();

		gpu_complex[in]->restoreCoords(1);

		// run remc
		REMC(gpu_complex[in], para, false, outfiles);

		gpu_complex[in]->clearMoves();
		gpu_complex[in]->setBest();
		gpu_complex[in]->createTrialCoords();
		gpu_complex[in]->calculateEnergy();

		for (int rn = 1; rn <= para->remc_replicas; rn++)
			cout << "Replica" << setw(4) << rn << " acceptance ratio:" << fixed << setprecision(6) <<
				setw(10) << gpu_complex[in]->getAcceptanceENE(rn) << endl;

		time(&t_end2);

		double dif2 = difftime(t_end2, t_start2);

		cout << "Docking time: " << fixed << setprecision(0) << dif2 << " s" << endl;
	}

	cout << "total decoys points generated:	" << outfiles->all << endl;
	// cout << "low decoys points generated:    " << outfiles->low << endl;
	// cout << "high decoys points generated:   " << outfiles->high << endl;

	time(&t_end1);

	cout << endl;

	printTime(difftime(t_end1, t_start1));

#if WRITE
	closeOutStream(outfiles);
#endif

	delete[]outfiles;
	delete[]para;

	return 0;
}
