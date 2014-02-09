
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "mpi_mcc.h"
#include "../lib/debug.h"

using namespace std;


int main(int argc, char *argv[])
{

	string molid = "MOLID";

	if (argc < 2) {
		cout << " gpudocksm -p  <target protein structure, PDB>" << endl
		     << "           -l  <compound library in SD format>" << endl;
		exit(EXIT_SUCCESS);
	}

	string protein_name;
	bool protein_opt = false;

	string compounds_name;
	bool compounds_opt = false;

	string output_name;
	bool output_opt = false;

	string data_name;

	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-p") && i < argc) {
			protein_name = string(argv[i + 1]);
			protein_opt = true;
		}
		if (!strcmp(argv[i], "-l") && i < argc) {
			compounds_name = string(argv[i + 1]);
			compounds_opt = true;
		}
		if (!strcmp(argv[i], "-o") && i < argc) {
			output_name = string(argv[i + 1]);
			output_opt = true;
		}
		if (!strcmp(argv[i], "-d") && i < argc) {
			data_name = string(argv[i + 1]);
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
	} else {
		ofstream outprot(output_name.c_str());
		outprot.close();
	}

	list < string > l1_sdf;
	list < string >::iterator i1_sdf;

	string line1;

	ifstream compounds_file(compounds_name.c_str());

	if (!compounds_file.is_open()) {
		cout << "Cannot open " << compounds_name << endl;
		exit(EXIT_FAILURE);
	}

	while (getline(compounds_file, line1)) {
		l1_sdf.push_back(line1);
	}
	compounds_file.close();

	int gpu_n = 0;

	string lib1[MAXSDF];
	int lib2 = 0;

	Complex *gpu_complex[MAXLIB];
	for (i1_sdf = l1_sdf.begin(); i1_sdf != l1_sdf.end(); i1_sdf++) {
		//cout<<"test5.5"; 
		lib1[lib2++] = (*i1_sdf);
		//cout< cout<endl;
		//cout<<lib1[lib2]<<endl;
		if ((*i1_sdf) == "$$$$") {
			if (lib2 > 10) {
				gpu_complex[gpu_n] = new Complex(0, 0);

				bool load1 = gpu_complex[gpu_n]->loadProtein(protein_name);
				if (load1) {
					cout << "Cannot read target protein structure" << endl;
					exit(EXIT_FAILURE);
				}

				bool load2 = false;
				load2 = gpu_complex[gpu_n]->loadParams(data_path);
				if (load2) {
					cout << "Cannot read parameter file" << endl;
					exit(EXIT_FAILURE);
				}

				bool load3 = gpu_complex[gpu_n]->loadLigand(lib1, lib2, molid);
				if (load3) {
					cout << "Cannot load ligand data" << endl;
					exit(EXIT_FAILURE);
				}
				gpu_n++;
			}

			lib2 = 0;

		}
	}

	// gpu_n = 1
	for (int in = 0; in < gpu_n; in++) {
		// gpu_complex[in]->createContactMatrix();
		// Decoy(gpu_complex[in], data_name.c_str(), output_name.c_str());
	  string ligand_id = gpu_complex[in]->getLigandID();
	  int pnp = gpu_complex[in]->getProteinPointsTotal();
	  int pnr = gpu_complex[in]->getProteinResiduesTotal();
	  int pens_total = gpu_complex[in]->getProteinEnsembleTotal();
	  int lna = gpu_complex[in]->getLigandAtomsTotal();
	  int lnb = gpu_complex[in]->getLigandBondsTotal();
	  int lens_total = gpu_complex[in]->getLigandEnsembleTotal();
	  printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", ligand_id, pnp, pnr, pens_total, lna, lnb, lens_total);

	}


	return 0;
}

/*
std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
*/
