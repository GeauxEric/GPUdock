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


#ifndef __SIZE_H_
#define __SIZE_H_

const int MAXPRO = 3000;  /* protein residues */
const int MAXLIG = 100;   /* ligand heavy atoms */
const int MAXEN1 = 3;     /* protein confs */
const int MAXEN2 = 20;    /* ligand confs */
const int MAXREP = 1000;  /* max replicas */
const int MAXSWP = 2000;  /* max swapping pairs */
const int MAXLIB = 100;   /* library cmps */
const int MAXSDF = 1000;  /* sdf length */
const int MAXTP1 = 30;    /* point types */
const int MAXTP2 = 24;    /* atom types */
const int MAXTP3 = 50;    /* point types (for ele) */
const int MAXTP4 = 20;    /* residue types */
const int MAXFP1 = 1024;  /* smiles */
const int MAXFP2 = 168;   /* maccs */
const int MAXWEI = 9;     /* energy terms */
const int MAXKDE = 10000; /* kde points */
const int MAXMCS = 500;   /* mcs fields */
const int MAXPOS = 1000;  /* position restraints */

const double PI = 3.14159265;

#endif
