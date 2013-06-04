#include <stdio.h>
#include <stdlib.h>

#include "sim.h"
/*
#define SZ 6
#define SZ_CUBE (SZ*SZ*SZ)
#define ITER_SWAP 4000000
#define ITER_SWAP_KERN 1000
#define NBETA 24
#define STR_LENG 64
*/

// analysis realization_[REALIZATION]_e.txt

//#define RUNS 40
//#define SAMPLES (RUNS*GD/2)
//#define SAMPLES 35840
#define SAMPLES 30720


int
main (int argc, char **argv)
{
  char myfile[STR_LENG]={0};
  int rec_start=1000;
  int rec_max = ITER_SWAP / ITER_SWAP_KERN;
  printf("rec_max=%d\n",rec_max);
  double tmp_q;
  double qq;
  //int e[GD][2000][NBETA];
  /*
  float *q;
  q=(float *)malloc(sizeof(float)*SAMPLES*rec_max*NBETA);
  if(q==NULL){
    printf("allocation failed");
    return(0);
  }
  */
  //free(e);
  double q1[NBETA];
  double q2[NBETA];
  double q3[NBETA];
  double q4[NBETA];
  int pair, rec, beta;


  for(beta = 0 ; beta< NBETA; beta++){
    q1[beta]=0.0f;
    q2[beta]=0.0f;
    q3[beta]=0.0f;
    q4[beta]=0.0f;
  }

  /// read from file to e[GD][rec_max][NBETA]
  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"q_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    


    for (rec = 0; rec< rec_start; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf", &tmp_q);
      }
      fscanf (fp, "\n");
    }

    for (rec = rec_start; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf", &tmp_q);
	qq=tmp_q/SZ_CUBE;
	q1[beta]+=qq;
	q2[beta]+=qq*qq;
	q3[beta]+=qq*qq*qq;
	q4[beta]+=qq*qq*qq*qq;

	//	q[pair*rec_max*NBETA+rec*NBETA+beta]=tmp_q;
      }
      fscanf (fp, "\n");
    }
    
    fclose (fp);
  }
    
    
    
  /// compute summation
  /*  
  for (rec = rec_start; rec < rec_max; rec++)
    for (pair = 0; pair < SAMPLES; pair++)
      for(beta = 0 ; beta< NBETA; beta++)
      {
	tmp_q= q[pair * rec_max * NBETA + rec * NBETA + beta];

	//	q4_avrg[rec]+=qq*qq*qq*qq/SAMPLES;	
      }
  free(q);
  */
  for(beta=0;beta<NBETA;beta++){
    double aq1=q1[beta]/SAMPLES/(rec_max-rec_start);
    double aq2=q2[beta]/SAMPLES/(rec_max-rec_start);
    double aq3=q3[beta]/SAMPLES/(rec_max-rec_start);
    double aq4=q4[beta]/SAMPLES/(rec_max-rec_start);
    double binder=1.5-0.5*(aq4-4.0*aq3*aq1+6.0*aq2*aq1*aq1-3.0*aq1*aq1*aq1*aq1)/(aq2-aq1*aq1)/(aq2-aq1*aq1);//0.5*(3-q4_avrg[rec]/(q2_avrg[rec]*q2_avrg[rec]));
    printf("%1.3f\n",binder);
  }

  return 0;
}


