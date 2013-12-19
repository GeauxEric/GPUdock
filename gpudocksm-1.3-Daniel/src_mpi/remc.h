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


#ifndef __REMC_H_
#define __REMC_H_

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<gsl/gsl_randist.h>

#include "size.h"
#include "complex.h"
#include "montecarlo.h"

using namespace std;

////////////////////////////////////////////////////ADDED BY DANIEL/////////
struct replica
{ 
 double conf1[MAXLIG][3];
 double conf2[22];
 int num;
};

////////////////////////////////////////////////////////////////////////////

void REMC( Complex *, int, int, int, double, double, double, double, bool );

#endif

