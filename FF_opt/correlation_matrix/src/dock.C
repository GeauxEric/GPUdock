#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <stdio.h>
#include <time.h>

#include "size.h"
#include "dock.h"
#include "run.h"
#include "util.h"
#include "load.h"
#include "corr_matx_gen.h"

using namespace std;

int
main (int argc, char **argv)
{
  cout << "------------------------------------------------------------" << endl
       << "                         GPU-dockSM-linear" << endl
       << "                         version 1.2" << endl << endl
       << "   GPU-accelerated mixed-resolution ligand docking using" << endl
       << "                Replica Exchange Monte Carlo" << endl
       << "------------------------------------------------------------" << endl << endl;

  srand (time (0));


  McPara *mcpara = new McPara;
  McLog *mclog = new McLog;
  ExchgPara *exchgpara = new ExchgPara;
  InputFiles *inputfiles = new InputFiles[1];
  OutputFiles * outfile = new OutputFiles[1];

  ParseArguments (argc, argv, mcpara, exchgpara, inputfiles, outfile);
  
  ofstream of;
  of.open(outfile->corr_mat_path.c_str());


  // load into preliminary data structures
  Ligand0 *lig0 = new Ligand0[MAXEN2];
  Protein0 *prt0 = new Protein0[MAXEN1];
  Psp0 *psp0 = new Psp0;
  Kde0 *kde0 = new Kde0;
  Mcs0 *mcs0 = new Mcs0[MAXPOS];
  EnePara0 *enepara0 = new EnePara0;

  // corellation matrix
  int total_rows = inputfiles->track_file.total_rows;
  int total_cols = TOTAL_COL;
  float *track_mat = new float[total_rows * total_cols];
  float *correlation_mat = new float[total_rows * total_rows];
  

  // loading
  loadTrack(&inputfiles->track_file, track_mat);
  loadLigand (&inputfiles->lig_file, lig0);
  loadProtein (&inputfiles->prt_file, prt0);
  loadLHM (&inputfiles->lhm_file, psp0, kde0, mcs0);
  loadEnePara (&inputfiles->enepara_file, enepara0);
  loadWeight(&inputfiles->weight_file, enepara0);
  loadNorPara(&inputfiles->norpara_file, enepara0);


  // sizes
  ComplexSize complexsize;
  complexsize.n_prt = inputfiles->prt_file.conf_total;	// number of protein conf
  complexsize.n_tmp = exchgpara->num_temp;	// number of temperature
  complexsize.n_lig = inputfiles->lig_file.conf_total;	// number of ligand conf
  complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;
  complexsize.lna = inputfiles->lig_file.lna;
  complexsize.pnp = inputfiles->prt_file.pnp;
  complexsize.pnk = kde0->pnk;
  complexsize.pos = inputfiles->lhm_file.pos;	// number of MCS positions

  // data structure optimizations 
  Ligand *lig = new Ligand[complexsize.n_rep];
  Protein *prt = new Protein[complexsize.n_prt];
  Psp *psp = new Psp;
  Kde *kde = new Kde;
  Mcs *mcs = new Mcs[complexsize.pos];
  EnePara *enepara = new EnePara;
  Temp *temp = new Temp[complexsize.n_tmp];
  Replica *replica = new Replica[complexsize.n_rep];


  OptimizeLigand (lig0, lig, complexsize);

  OptimizeProtein (prt0, prt, enepara0, lig0, complexsize);

  OptimizePsp (psp0, psp, lig, prt);
  OptimizeKde (kde0, kde);
  OptimizeMcs (mcs0, mcs, complexsize);
  OptimizeEnepara (enepara0, enepara);


  delete[]lig0;
  delete[]prt0;
  delete[]psp0;
  delete[]kde0;
  delete[]mcs0;
  delete[]enepara0;


  // initialize system
  InitLigCoord (lig, complexsize);
  SetTemperature (temp, exchgpara);
  // SetTemperature (temp, mcpara, complexsize);
  SetReplica (replica, lig, complexsize);
  SetMcLog (mclog);


  // debug
  // PrintDataSize (lig, prt, psp, kde, mcs, enepara);
  // PrintLigand (lig);
  // PrintProtein (prt);

  clock_t begin, end;
  double time_spent;
  begin = clock();

  cout << "begin correlation matrix generation" << endl;

  // generate correllation matrix
  GenCorrMat (correlation_mat, track_mat, total_rows, lig, prt, enepara);
  
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  // flip the matrix
  for (int i = 0; i < total_rows; i++)
    for (int j = 0; j < i; j++) 
      correlation_mat[i*total_rows +j] = correlation_mat[j*total_rows +i];
  
  // write to file
  of << total_rows << endl;

  for (int i = 0; i < total_rows; i++) {
    for (int j = 0; j < total_rows; j++) {
      of << correlation_mat[i*total_rows + j] << " ";
    }
    of << endl;
  }
  
  of.close();

  // clean up

  delete[]correlation_mat;
  delete[]outfile;
  delete[]track_mat;
  delete[]mcpara;
  delete[]mclog;
  delete[]inputfiles;
  delete[]lig;
  delete[]prt;
  delete[]psp;
  delete[]kde;
  delete[]mcs;
  delete[]enepara;
  delete[]temp;
  delete[]replica;
  delete[]exchgpara;

  return 0;
}
