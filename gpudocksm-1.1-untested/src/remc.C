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


#include "remc.h"
#include "montecarlo.h"

using namespace std;

// ==================================================================================   Replica Exchange Monte Carlo

void REMC( Complex * re_complex, int re_replicas, int re_steps, int re_cycles, double re_t, double re_r, double re_b, bool re_prt )
{
 Replica replicas[MAXREP];
 
 int replicas_tot = 0;
 
 int pairs[MAXSWP][3];
 
 int pairs_tot = 0;
 
 int swaps_tot = 0;
 int swaps_acc = 0;
 
 double df = 6.0;
 
 cout << "Creating replicas ... " << flush;
 
 for ( int ii1 = 0; ii1 < re_complex->getProteinEnsembleTotal(); ii1++ )
  for ( int ii2 = 0; ii2 < re_complex->getLigandEnsembleTotal(); ii2++ )
   for ( int ii3 = 0; ii3 < re_replicas; ii3++ )
   {
    replicas[replicas_tot].temp = 1.1;
    
    if ( ii3 )
     replicas[replicas_tot].temp = exp( ( (double) ii3 ) / sqrt(df) );
    
    for ( int j = 0; j < 3; j++ )
     replicas[replicas_tot].conf1[j] = re_t * unirand( -1, 1 );
    
    for ( int j = 3; j < 6; j++ )
     replicas[replicas_tot].conf1[j] = ( ( re_r * unirand( -1, 1 ) ) * PI ) / 180.0;
    
    replicas[replicas_tot].conf2[0] = ii1;
    replicas[replicas_tot].conf2[1] = ii2;
    
    /* create temporary configuration */
    
    re_complex->setConfiguration( replicas[replicas_tot].conf1 );
    
    re_complex->setProteinEnsembleCurrent( replicas[replicas_tot].conf2[0] );
    re_complex->setLigandEnsembleCurrent( replicas[replicas_tot].conf2[1] );
    
    re_complex->setTemperature( replicas[replicas_tot].temp );
    
    re_complex->calculateEnergy();
    
    replicas[replicas_tot].energy = re_complex->getEnergy(1);
    
    replicas[replicas_tot].number = replicas_tot;
    
    replicas[replicas_tot].steps_tot = 0;
    replicas[replicas_tot].steps_acc = 0;
    
    if ( ii3 < re_replicas - 1 )
    {
     pairs[pairs_tot][0] = replicas_tot;
     pairs[pairs_tot][1] = replicas_tot + 1;
     pairs[pairs_tot][2] = ii3%2;
     
     pairs_tot++;
    }
    
    replicas_tot++;
   }
 
 cout << "done" << endl;
 
 cout << "System comprises " << replicas_tot << " replicas and " << pairs_tot << " swapping pairs" << endl;
 
 cout << "Running REMC ... 0.0%" << flush;
 
 for (int c = 0; c < re_cycles; c++)
 {
  for ( int r = 0; r < replicas_tot; r++ )
   MonteCarlo( re_complex, re_steps, re_t, re_r, re_b, re_prt, replicas[r].temp, replicas[r].conf1, replicas[r].conf2, replicas[r].energy, replicas[r].number, replicas[r].steps_tot, replicas[r].steps_acc );
  
  for ( int p = 0; p < pairs_tot; p++ )
   if ( pairs[p][2] == c%2 )
   {
    swaps_tot++;
    
    bool swap = false;
    
    if ( replicas[pairs[p][0]].energy > replicas[pairs[p][1]].energy )
     swap = true;
    else if ( unirand( 0, 1 ) < exp( -1.0 * ( ( ( 1.0 / ( re_b * replicas[pairs[p][0]].temp ) ) - ( 1.0 / ( re_b * replicas[pairs[p][1]].temp ) ) ) * ( replicas[pairs[p][1]].energy - replicas[pairs[p][0]].energy ) ) ) )
     swap = true;
    
    if ( swap )
    {
     for ( int j = 0; j < 6; j++ )
     {
      double swap_tmp = replicas[pairs[p][0]].conf1[j];
      
      replicas[pairs[p][0]].conf1[j] = replicas[pairs[p][1]].conf1[j];
      
      replicas[pairs[p][1]].conf1[j] = swap_tmp;
     }
     
     double ene_tmp = replicas[pairs[p][0]].energy;
     
     replicas[pairs[p][0]].energy = replicas[pairs[p][1]].energy;
     
     replicas[pairs[p][1]].energy = ene_tmp;
     
     swaps_acc++;
    }
   }
  
  cout << "\r" << "Running REMC ... " << setprecision(1) << ( (double) c / (double) re_cycles ) * 100 << "%" << flush;
 }
 
 cout << "\r" << "Running REMC ... 100.0%" << endl;
}
