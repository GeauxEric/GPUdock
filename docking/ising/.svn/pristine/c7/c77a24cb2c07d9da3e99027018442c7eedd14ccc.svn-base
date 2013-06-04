#include <stdio.h>
#include <stdlib.h>

#include "sim.h"


// analysis realization_[REALIZATION]_e.txt



int
main (int argc, char **argv)
{
  char myfile[STR_LENG]={0};

  int rec_max = ITER_SWAP / ITER_SWAP_KERN;
  printf("rec_max=%d\n",rec_max);
  float tmp_e;
  //int e[GD][2000][NBETA];
   float *e;
   e=(float *)malloc(sizeof(float)*GD*rec_max*NBETA);
   if(e==NULL){
    printf("allocation failed");
    return(0);
   }
  //free(e);
   //double e_avrg[rec_max][NBETA];
   double e_avrg[NBETA];
  int realization, rec, beta;


  /// read from file to e[GD][rec_max][NBETA]
    for (realization = 0; realization < GD; realization++) {
      if(myfile!=NULL)
	int charcheck = sprintf ( myfile ,"realization_%04d_e.txt", realization);
        FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
    exit (1);
    }
    
    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%f", &tmp_e);
	 e[realization*rec_max*NBETA+rec*NBETA+beta]=tmp_e;
      }
      fscanf (fp, "\n");
    }

    fclose (fp);
  }



  /// compute summation

    //  for (rec = 0; rec < rec_max; rec++)
    for (beta = 0; beta < NBETA; beta++)
      e_avrg[beta] = 0;
    
    for (rec = 10; rec < rec_max; rec++)
    for (realization = 0; realization < GD; realization++)
      for (beta = 0; beta < NBETA; beta++)
	e_avrg[beta] += (double) e[realization*rec_max*NBETA+rec*NBETA+beta];

  free(e);

  /// print
  
  /*
  for (beta = 0; beta < NBETA; beta++)
    printf ("   b%02d ", beta);
  printf ("\n");

  for (beta = 0; beta < NBETA; beta++)
    printf ("-------");
  printf ("\n");
  */
  //  for (rec = 0; rec < rec_max; rec++) {
  //  printf ("%02d\t", rec*ITER_SWAP_KERN);
    for (beta = 0; beta < NBETA; beta++) {
      e_avrg[beta] /= (GD*1000);
      printf ("%f\n", e_avrg[beta]);
    }
    //printf ("\n");
    // printf ("\t iteration %02d\n", rec*ITER_SWAP_KERN);
    //  }
  /*
  for (beta = 0; beta < NBETA; beta++)
    printf ("-------");
  printf ("\n");
  */

  return 0;
}


