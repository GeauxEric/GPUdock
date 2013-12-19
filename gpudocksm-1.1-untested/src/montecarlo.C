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

using namespace std;

// ==================================================================================   Monte Carlo

void MonteCarlo( Complex * mc_complex, int mc_steps, double mc_t, double mc_r, double mc_b, bool mc_prt, double mc_temp, double mc_conf1[], int mc_conf2[], double &mc_energy, int mc_number, int &mc_tot, int &mc_acc )
{
 for ( int step = 0; step < mc_steps; step++ )
 {
  double energy_old = mc_energy;
  
  double mc_conf3[6];
  
  for ( int i = 0; i < 3; i++ )
  {
   mc_conf3[i] = mc_conf1[i];
   
   mc_conf1[i] += mc_t * unirand( -1, 1 );
  }
  
  for ( int i = 3; i < 6; i++ )
  {
   mc_conf3[i] = mc_conf1[i];
   
   mc_conf1[i] += ( ( mc_r * unirand( -1, 1 ) ) * PI ) / 180.0;
  }
  
  mc_complex->setConfiguration( mc_conf1 );
  
  mc_complex->setProteinEnsembleCurrent( mc_conf2[0] );
  mc_complex->setLigandEnsembleCurrent( mc_conf2[1] );
  
  mc_complex->setTemperature( mc_temp );
  
  mc_complex->calculateEnergy();
  
  mc_energy = mc_complex->getEnergy(1);
  
  if ( mc_energy > energy_old )
   if ( unirand( 0, 1 ) > exp( ( -1.0 * ( mc_energy - energy_old ) ) / ( mc_b * mc_temp ) ) )
   {
    
    for ( int i = 0; i < 6; i++ )
     mc_conf1[i] = mc_conf3[i];
    
    mc_energy = energy_old;
    
    mc_acc--;
   }
  
  mc_tot++;
  mc_acc++;
 }
}
