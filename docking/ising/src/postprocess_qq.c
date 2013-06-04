#include <stdio.h>
#include <stdlib.h>

#include "sim.h"
/*
#define SZ 16
#define SZ_CUBE (SZ*SZ*SZ)
#define ITER_SWAP 20000000
#define ITER_SWAP_KERN 1000
#define NBETA 24
#define STR_LENG 64
*/


// analysis realization_[REALIZATION]_e.txt

//#define RUNS 60
//#define SAMPLES (RUNS*GD/2)
#define SAMPLES 3000

int
main (int argc, char **argv)
{
  char myfile[STR_LENG]={0};
  int rec_start=0;
  int rec_max = ITER_SWAP / ITER_SWAP_KERN;
  printf("rec_max=%d\n",rec_max);
  float tmp_q;
  double qq;
  //int e[GD][2000][NBETA];
  float *q;
  double *q1,*q2,*q3,*q4;
  q1=(double *)malloc(sizeof(double)*rec_max);
  q2=(double *)malloc(sizeof(double)*rec_max);
  q3=(double *)malloc(sizeof(double)*rec_max);
  q4=(double *)malloc(sizeof(double)*rec_max);
  q=(float *)malloc(sizeof(float)*SAMPLES*NBETA*rec_max);
  if(q==NULL){
    printf("allocation failed");
    return(0);
  }
  //free(e);
  /*
  double q1[rec_max];
  double q2[rec_max];
  double q3[rec_max];
  double q4[rec_max];
  */
  int pair, rec, beta;


  /// read from file to e[GD][rec_max][NBETA]
  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"q_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    
    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%f", &tmp_q);
	q[pair*rec_max*NBETA+rec*NBETA+beta]=tmp_q;
      }
      fscanf (fp, "\n");
    }      
    fclose (fp);
  }



  /// compute summation
    
  for (rec = 0; rec < rec_max; rec++)
    {
      q1[rec]=0.0f;
      q2[rec]=0.0f;
      q3[rec]=0.0f;
      q4[rec]=0.0f;
    }
  

  for (pair = 0; pair < SAMPLES; pair++)
    for (rec = 0; rec < rec_max; rec++)
      {
	tmp_q= q[pair*rec_max*NBETA+rec*NBETA+NBETA-1];
	qq=(double)tmp_q/SZ_CUBE;
	q1[rec]+=qq;
	q2[rec]+=qq*qq;
	q3[rec]+=qq*qq*qq;
	q4[rec]+=qq*qq*qq*qq;
	//	q4_avrg[rec]+=qq*qq*qq*qq/SAMPLES
      }
  free(q);
  
  for(rec=0;rec<rec_max;rec++){
    double aq1=q1[rec]/SAMPLES;
    double aq2=q2[rec]/SAMPLES;
    double aq3=q3[rec]/SAMPLES;
    double aq4=q4[rec]/SAMPLES;

    double binder=1.5-0.5*(aq4-4*aq3*aq1+6*aq2*aq1*aq1-3*aq1*aq1*aq1*aq1)/(aq2-aq1*aq1)/(aq2-aq1*aq1);//0.5*(3-q4_avrg[rec]/(q2_avrg[rec]*q2_avrg[rec]));
    printf("%7d\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n",(rec+1)*ITER_SWAP_KERN,binder,aq1,aq2,aq2,aq4);
  }

  free(q1);
  free(q2);
  free(q3);
  free(q4);

//printf("\n");
  /// print
  
  /*
  for (beta = 0; beta < NBETA; beta++)
    printf ("   b%02d ", beta);
  printf ("\n");

  for (beta = 0; beta < NBETA; beta++)
    printf ("-------");
  printf ("\n");
  
  for (rec = 0; rec < rec_max; rec++) {
    printf ("%02d\t", rec*ITER_SWAP_KERN);
    for (beta = 0; beta < NBETA; beta++) {
      e_avrg[rec][beta] /= GD;
      printf ("%7d", (int) e_avrg[rec][beta]);
    }
    printf ("\n");
    // printf ("\t iteration %02d\n", rec*ITER_SWAP_KERN);
  }
  
  for (beta = 0; beta < NBETA; beta++)
    printf ("-------");
  printf ("\n");
  */

  return 0;
}


