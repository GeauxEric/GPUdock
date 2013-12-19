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


#include "montecarlo.h"
#include "remc.h"

using namespace std;

// ==================================================================================   Monte Carlo

double MonteCarlo( gsl_rng*& r, Complex * mc_complex, int mc_steps, double mc_t, double mc_r, double mc_d, double mc_b, int mc_n, bool mc_prt )
{
 double dx[3];
 size_t n = 3;

 double scaled_t = mc_t * sqrt( mc_complex->getTemperature() );
 double scaled_r = mc_r * sqrt( mc_complex->getTemperature() );

 double m_conf1[MAXLIG][3];
 double m_conf2[22];

 //Output File
 ofstream decoysH, MCcoefH, decoyenergiesH, decoysL, MCcoefL, decoyenergiesL;
 MCcoefH.open("MCcoefsH.txt", fstream::app);
 decoysH.open("decoysH.txt", fstream::app);
 decoyenergiesH.open("decoyeneH.txt", fstream::app);
 MCcoefL.open("MCcoefsL.txt", fstream::app);
 decoysL.open("decoysL.txt", fstream::app);
 decoyenergiesL.open("decoyeneL.txt", fstream::app);

 for ( int i = 0; i < mc_steps; i++ )
 {
  mc_complex->getConfCoordsPerm(m_conf1);
  mc_complex->getConfParams(m_conf2);
  mc_complex->restoreCoords(1);
  mc_complex->clearMoves();
  mc_complex->createTrialCoords();
  mc_complex->calculateEnergy();

  double ene1 = mc_complex->getEnergy(1);
  gsl_ran_dir_nd (r, n, dx);

  mc_complex->ligTranslate(1, dx[0] * scaled_t );
  mc_complex->ligTranslate(2, dx[1] * scaled_t );
  mc_complex->ligTranslate(3, dx[2] * scaled_t );

  if ( gsl_rng_uniform(r) < 0.5 && mc_complex->getLigandEnsembleTotal() > 1 )
   mc_complex->ligPerturb( gsl_rng_uniform_int(r, mc_complex->getLigandEnsembleTotal()) );

  gsl_ran_dir_nd (r, n, dx);

  mc_complex->ligRotate(dx[0] * scaled_r, dx[1]*scaled_r, dx[2]*scaled_r );
  
  if ( gsl_rng_uniform(r) < 0.2)
   mc_complex->setProteinEnsembleCurrent( gsl_rng_uniform_int(r, mc_complex->getProteinEnsembleTotal()));

  if ( gsl_rng_uniform(r) < 0.40)
   mc_complex->setLigandEnsembleCurrent( gsl_rng_uniform_int(r, mc_complex->getLigandEnsembleTotal()));

  mc_complex->createTrialCoords();
  mc_complex->calculateEnergy();
  double mcc;
  mcc=mc_complex->getMCC();
  if (mcc<0.5){
    mc_complex->printLigandInfo(decoysL, MCcoefL, decoyenergiesL);
  }
  if (mcc>=0.8){
    mc_complex->printLigandInfo(decoysH, MCcoefH, decoyenergiesH);
  }
  double ene2 = mc_complex->getEnergy(1);
  bool w = false;

  if ( ene2 > ene1 )
  {
   double pr1 = exp( ( -1.0 * ( ene2 - ene1 ) ) / (((double) mc_b * (double)(mc_complex->getTemperature()) ) ));
   if ( gsl_rng_uniform(r) > pr1 )
   {
   mc_complex->setConfParams(m_conf2);
   mc_complex->restoreCoords(1);
    ene2 = ene1;
    }
   else
   {

     mc_complex->getConfCoords(m_conf1);
     mc_complex->getConfParams(m_conf2);
     mc_complex->setConfCoordsPerm(m_conf1);
    w = true;
   }
  }
  else
  {

   mc_complex->getConfCoords(m_conf1);
   mc_complex->getConfParams(m_conf2);
   mc_complex->setConfCoordsPerm(m_conf1);
   w = true;
  }
  mc_complex->restoreCoords(1);
  mc_complex->addAcceptanceENE(mc_n, w);
 
 /*
  if (w = (true) ){
   mc_complex->printLigandInfo(decoys);
  }
 */

 }
 cout<<mc_complex->getEnergy(1)<<endl;
 decoysL.close();
 decoysH.close();
 MCcoefL.close();
 MCcoefH.close();
 decoyenergiesL.close();
 decoyenergiesH.close();
 return mc_complex->getEnergy(1);
}

                       
