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
#include "complex.h"
using namespace std;

// ==================================================================================   Replica Exchange Monte Carlo

void REMC( Complex * re_complex, int re_replicas, int re_steps, int re_cycles, double re_t, double re_r, double re_d, double re_b, bool re_prt )
{
 vector<replica> replicas;

 ofstream pcoords;
 pcoords.open("pcoords.txt");
// re_complex->printProteinInfo(pcoords);
 pcoords.close();


 double df = 6.0;

 const gsl_rng_type * T;
 gsl_rng * rn;
 gsl_rng_env_setup();
 T  = gsl_rng_default;
 rn = gsl_rng_alloc(T);

 cout << "Creating replicas ..... " << flush;

 for (int r1 = 0; r1 < re_replicas; r1++)
 {
  double t_conf1[MAXLIG][3];
  double t_conf2[22];

   re_complex->clearMoves();

  if ( r1 > 0 && re_complex->getLigandEnsembleTotal() > 1 ){
   re_complex->setTemperature(exp(4.0/((double)re_replicas+1.0-(double)r1)));
   re_complex->ligPerturb( gsl_rng_uniform_int ( rn, re_complex->getLigandEnsembleTotal() ) );
  }
  re_complex->createTrialCoords();
  re_complex->calculateEnergy();

  re_complex->getConfCoords(t_conf1);
  re_complex->getConfParams(t_conf2);

  if ( r1 == 0  )
   t_conf2[0] = 1.1;
  else
   t_conf2[0] = exp( ( (double) r1 ) / sqrt(df) );

 replica re_new;

  for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
   for ( int j = 0; j < 3; j++ )
    re_new.conf1[i][j] = t_conf1[i][j];

  for ( int i = 0; i < 22; i++ )
   re_new.conf2[i] = t_conf2[i];

  re_new.num = r1 + 1;

  replicas.push_back(re_new);

  re_complex->restoreCoords(1);

  cout << re_new.num << " " << flush;

 }

 cout << "done" << endl;

 cout << "Minimizing replicas ... " << flush;

 vector<replica>::iterator r2;

 for ( r2 = replicas.begin() ; r2 != replicas.end(); r2++ )
 {
  double t_conf1[MAXLIG][3];
  double t_conf2[22];

  for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
   for ( int j = 0; j < 3; j++ )
    t_conf1[i][j] = (*r2).conf1[i][j];

  for ( int i = 0; i < 22; i++ )
   t_conf2[i] = (*r2).conf2[i];

  re_complex->setConfParams(t_conf2);
  re_complex->setConfCoords(t_conf1);

  //Simplex(re_complex, si_steps, re_t, re_r, false); //////////////// MAY NOT NEED THIS SECTION IF NO SIMPLEX/////////////

  re_complex->getConfCoords(t_conf1);
  re_complex->getConfParams(t_conf2);

  for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
   for ( int j = 0; j < 3; j++ )
    (*r2).conf1[i][j] = t_conf1[i][j];

  for ( int i = 0; i < 22; i++ )
   (*r2).conf2[i] = t_conf2[i];

  re_complex->restoreCoords(1);

  cout << (*r2).num << " " << flush;
 }

 cout << "done" << endl;

 cout << "REMC .................. " << flush;

 for (int c = 0; c < re_cycles; c++)
 {
  cout << c + 1 << " " << flush<<endl;
  
  for ( r2 = replicas.begin() ; r2 != replicas.end(); r2++ )
  {
   double t_conf1[MAXLIG][3];
   double t_conf2[22];

   for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
    for ( int j = 0; j < 3; j++ )
     t_conf1[i][j] = (*r2).conf1[i][j];

   for ( int i = 0; i < 22; i++ )
    t_conf2[i] = (*r2).conf2[i];

   re_complex->setConfParams(t_conf2);
   re_complex->setConfCoordsPerm(t_conf1);

   MonteCarlo(rn, re_complex, re_steps, re_t, re_r, re_d, re_b, (*r2).num, true);

   re_complex->getConfCoordsPerm(t_conf1);
   re_complex->getConfParams(t_conf2);

   for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
    for ( int j = 0; j < 3; j++ )
     (*r2).conf1[i][j] = t_conf1[i][j];

   for ( int i = 0; i < 22; i++ )
    (*r2).conf2[i] = t_conf2[i];

   re_complex->restoreCoords(1);
  }

  int swp0 = 0;

  if ( c%2 == 0 )
  {
   for ( int i = 0; i < re_replicas - 1; i += 2)
   {
    unsigned int i1 = i;
    unsigned int i2 = i + 1;

    bool swp1 = false;

    if ( replicas[i1].conf2[1] > replicas[i2].conf2[1] )
    {
     swp1 = true;
    }
    else
    {
     double pr1 = exp( -1.0 * ( ( ( 1.0 / ( re_b * replicas[i1].conf2[0] ) ) - ( 1.0 / ( re_b * replicas[i2].conf2[0] ) ) ) * ( replicas[i2].conf2[1] - replicas[i1].conf2[1] ) ) );

     if ( gsl_rng_uniform(rn) < pr1 )
      swp1 = true;
    }
    
    if ( swp1 )
    {
     double tmp1t = replicas[i1].conf2[0];
     double tmp2t = replicas[i2].conf2[0];

     int tmp1e = replicas[i1].num;
     int tmp2e = replicas[i2].num;

     replica swp2;

     swp2 = replicas[i1];

     replicas[i1] = replicas[i2];

     replicas[i2] = swp2;

     replicas[i1].conf2[0] = tmp1t;
     replicas[i2].conf2[0] = tmp2t;

     replicas[i1].num = tmp1e;
     replicas[i2].num = tmp2e;

     swp0++;
    }

   re_complex->addAcceptanceSWP(swp1);
   }
  }
  else
  {
   for ( int i = 1; i < re_replicas - 1; i += 2)
   {
    unsigned int i1 = i;
    unsigned int i2 = i + 1;

    bool swp1 = false;

    if ( replicas[i1].conf2[1] > replicas[i2].conf2[1] )
    {
     swp1 = true;
    }
    else
    {
     double pr1 = exp( -1.0 * ( ( ( 1.0 / ( re_b * replicas[i1].conf2[0] ) ) - ( 1.0 / ( re_b * replicas[i2].conf2[0] ) ) ) * ( replicas[i2].conf2[1] - replicas[i1].conf2[1] ) ) );

     if ( gsl_rng_uniform(rn) < pr1 )
      swp1 = true;
    }
   
        if ( swp1 )
    {
     double tmp1t = replicas[i1].conf2[0];
     double tmp2t = replicas[i2].conf2[0];

     int tmp1e = replicas[i1].num;
     int tmp2e = replicas[i2].num;

     replica swp2;

     swp2 = replicas[i1];

     replicas[i1] = replicas[i2];

     replicas[i2] = swp2;

     replicas[i1].conf2[0] = tmp1t;
     replicas[i2].conf2[0] = tmp2t;

     replicas[i1].num = tmp1e;
     replicas[i2].num = tmp2e;

     swp0++;
    }

    re_complex->addAcceptanceSWP(swp1);
   }
  }
for(int q=1;q<8;q++)
 cout<<re_complex->getEnergy(q)<<" ";
 cout<<endl;
 
 }

 cout << "done" << endl;
 /*
 cout << "Minimizing replicas ... " << flush;

 for ( r2 = replicas.begin() ; r2 != replicas.end(); r2++ )
 {
  double t_conf1[MAXLIG][3];
  double t_conf2[22];

  //re_complex->backupCoords(1);

  for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
   for ( int j = 0; j < 3; j++ )
    t_conf1[i][j] = (*r2).conf1[i][j];

  for ( int i = 0; i < 22; i++ )
   t_conf2[i] = (*r2).conf2[i];

  re_complex->setConfParams(t_conf2);
  re_complex->setConfCoords(t_conf1);
  //re_complex->setConfParams(t_conf2);

  //Simplex(re_complex, si_steps, re_t, re_r, false);

  re_complex->getConfCoords(t_conf1);
  re_complex->getConfParams(t_conf2);

  for ( int i = 0; i < re_complex->getLigandAtomsTotal(); i++ )
   for ( int j = 0; j < 3; j++ )
    (*r2).conf1[i][j] = t_conf1[i][j];

  for ( int i = 0; i < 22; i++ )
   (*r2).conf2[i] = t_conf2[i];

  re_complex->restoreCoords(1);

  cout << (*r2).num << " " << flush;
 }

 cout << "done" << endl;
 
*/
 
 

 cout << "Stop 1 " << endl;
 
// exit(0);
}
