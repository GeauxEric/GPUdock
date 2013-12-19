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

#include "size.h"
#include "complex.h"
#include "montecarlo.h"
#include "random.h"

using namespace std;

struct Replica
{
 double conf1[6]; /* translation + rotation */
 int    conf2[2]; /* protein + ligand ensemble */
 double     temp; /* temperature */
 double   energy; /* total energy */
 int      number; /* replica number */
 int   steps_tot; /* total MC steps */
 int   steps_acc; /* accepted MC steps */
};

void REMC( Complex *, int, int, int, double, double, double, bool );

#endif

