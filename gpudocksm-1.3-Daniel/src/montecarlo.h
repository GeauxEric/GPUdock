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


#ifndef __MONTECARLO_H_
#define __MONTECARLO_H_

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include <gsl/gsl_randist.h>

#include "complex.h"
#include "montecarlo.h"
#include "remc.h"

using namespace std;
 
double MonteCarlo( gsl_rng*&, Complex *, int, double, double, double, double, int, bool );


#endif
