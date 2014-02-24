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

#include <string.h>
#include "debug.h"
#include "montecarlo.h"
#include "remc.h"

using namespace std;

// ==================================================================================   Monte Carlo

double
MonteCarlo(gsl_rng * &r, Complex * mc_complex, Param * para, int mc_n, bool mc_prt, int data_part, OutStream * outfiles)
{
  int mc_steps = para->remc_steps;
  double mc_t = para->step_t;
  double mc_r = para->step_r;
  double mc_b = para->boltzmann;
  // int data_pt_total = para->data_pt_total;

  double dx[3];
  size_t n = 3;

  double scaled_t = mc_t * sqrt(mc_complex->getTemperature());
  double scaled_r = mc_r * sqrt(mc_complex->getTemperature());

  double m_conf1[MAXLIG][3];
  double m_conf2[22];

  // double low_bound = 0.4;
  // double high_bound = 0.6;

  for (int i = 0; i < mc_steps; i++) {
    mc_complex->getConfCoordsPerm(m_conf1);
    mc_complex->getConfParams(m_conf2);
    mc_complex->restoreCoords(1);
    mc_complex->clearMoves();

    // if (outfiles->all > data_pt_total) {
    // 	outfiles->all_stop = 0;
    // }
    // if (outfiles->high > data_pt_total) {
    // 	outfiles->high_stop = 0;
    // }
    // if (outfiles->low > data_pt_total) {
    // 	outfiles->low_stop = 0;
    // }

    double ini_mcc = mc_complex->getMCC();

    // record the native pose
    if (ini_mcc == 1.0) {
      // mc_complex->printLigandInfo(outfiles->all_decoys, outfiles->all_MCcoef, outfiles->all_decoyenergies);
      mc_complex->printEnerInfo(outfiles->all_decoyenergies);
      // mc_complex->printLigandInfo(outfiles->H_decoys, outfiles->H_MCcoef, outfiles->H_decoyenergies);
    }

    mc_complex->createTrialCoords();
    mc_complex->calculateEnergy();

    double ene1 = mc_complex->getEnergy(1);

    gsl_ran_dir_nd(r, n, dx);

    mc_complex->ligTranslate(1, dx[0] * scaled_t);
    mc_complex->ligTranslate(2, dx[1] * scaled_t);
    mc_complex->ligTranslate(3, dx[2] * scaled_t);

    if (gsl_rng_uniform(r) < 0.5 && mc_complex->getLigandEnsembleTotal() > 1)
      mc_complex->ligPerturb(gsl_rng_uniform_int(r, mc_complex->getLigandEnsembleTotal()));

    gsl_ran_dir_nd(r, n, dx);

    mc_complex->ligRotate(dx[0] * scaled_r, dx[1] * scaled_r, dx[2] * scaled_r);

    if (gsl_rng_uniform(r) < 0.2)
      mc_complex->setProteinEnsembleCurrent(gsl_rng_uniform_int(r, mc_complex->getProteinEnsembleTotal()));

    if (gsl_rng_uniform(r) < 0.40)
      mc_complex->setLigandEnsembleCurrent(gsl_rng_uniform_int(r, mc_complex->getLigandEnsembleTotal()));

    mc_complex->createTrialCoords();
    mc_complex->calculateEnergy();
    double mcc;
    mcc = mc_complex->getMCC();
    double ene2 = mc_complex->getEnergy(1);
    bool w = false;			// keep track of flipping or not

    if (ene2 > ene1) {
      double pr1 = exp((-1.0 * (ene2 - ene1)) / (((double)mc_b * (double)(mc_complex->getTemperature()))));
      if (gsl_rng_uniform(r) > pr1) {
	mc_complex->setConfParams(m_conf2);
	mc_complex->restoreCoords(1);
	ene2 = ene1;
      } else {
	mc_complex->getConfCoords(m_conf1);
	mc_complex->getConfParams(m_conf2);
	mc_complex->setConfCoordsPerm(m_conf1);
	w = true;
	// push to the buffer 
	if (outfiles->all_stop == 1) {
	  outfiles->all += 1;
	  // mc_complex->printLigandInfo(outfiles->all_decoys, outfiles->all_MCcoef, outfiles->all_decoyenergies);
	  mc_complex->printEnerInfo(outfiles->all_decoyenergies);
	}
	// if (mcc < low_bound && (outfiles->low_stop == 1)) {
	// 	outfiles->low += 1;
	// 	mc_complex->printLigandInfo(outfiles->L_decoys, outfiles->L_MCcoef, outfiles->L_decoyenergies);
	// }
	// if (mcc >= high_bound && (outfiles->high_stop == 1)) {
	// 	outfiles->high += 1;
	// 	mc_complex->printLigandInfo(outfiles->H_decoys, outfiles->H_MCcoef, outfiles->H_decoyenergies);
	// }

      }
    } else {

      mc_complex->getConfCoords(m_conf1);
      mc_complex->getConfParams(m_conf2);
      mc_complex->setConfCoordsPerm(m_conf1);
      w = true;
      // push to the buffer 
      if (outfiles->all_stop == 1) {
	outfiles->all += 1;
	// mc_complex->printLigandInfo(outfiles->all_decoys, outfiles->all_MCcoef, outfiles->all_decoyenergies);
	mc_complex->printEnerInfo(outfiles->all_decoyenergies);
      }
      // if (mcc < low_bound && (outfiles->low_stop == 1)) {
      // 	outfiles->low += 1;
      // 	mc_complex->printLigandInfo(outfiles->L_decoys, outfiles->L_MCcoef, outfiles->L_decoyenergies);
      // }
      // if (mcc >= high_bound && (outfiles->high_stop == 1)) {
      // 	outfiles->high += 1;
      // 	mc_complex->printLigandInfo(outfiles->H_decoys, outfiles->H_MCcoef, outfiles->H_decoyenergies);
      // }
    }

    mc_complex->restoreCoords(1);
    mc_complex->addAcceptanceENE(mc_n, w);

  }

  return mc_complex->getEnergy(1);
}
