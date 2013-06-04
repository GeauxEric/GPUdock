#include <stdio.h>
#include <stdlib.h>

#include "../sim.h"


// analysis realization_[REALIZATION]_e.txt

int
main (int argc, char **argv)
{
  char myfile[STR_LENG];

  int rec_max = ITER_SWAP / ITER_SWAP_KERN;

  int tmp_e;
  int e[GD][rec_max][NBETA];
  double e_avrg[GD][NBETA];
  int realization, rec, beta;


  /// read from file to e[GD][rec_max][NBETA]
  for (realization = 0; realization < GD; realization++) {

    sprintf (myfile, "realization_%04d_e.txt", realization);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%d", &tmp_e);
	e[realization][rec][beta] = tmp_e;
      }
      fscanf (fp, "\n");
    }

    fclose (fp);
  }


  /*
  // print out original data for verification
  for (realization = 23; realization < 24; realization++) {
    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	printf ("%7d", e[realization][rec][beta]);
      }
      printf ("\t %d\n", rec);
      //printf ("\n");
    }
  }
  */



  /// compute summation

  for (realization = 0; realization < GD; realization++)
    for (beta = 0; beta < NBETA; beta++)
      e_avrg[realization][beta] = 0;

  for (realization = 0; realization < GD; realization++)
    for (rec = 0; rec < rec_max; rec++)
      for (beta = 0; beta < NBETA; beta++)
	e_avrg[realization][beta] += (double) e[realization][rec][beta];



  /// print

  for (beta = 0; beta < NBETA; beta++)
    printf ("   b%02d   ", beta);
  printf ("\n");

  for (beta = 0; beta < NBETA; beta++)
    printf ("---------");
  printf ("\n");

  for (realization = 0; realization < GD; realization++) {
    for (beta = 0; beta < NBETA; beta++) {
      e_avrg[realization][beta] /= rec_max;
      printf ("%9.1f", e_avrg[realization][beta]);
    }
    //printf ("\n");
    printf ("\t realization %02d\n", realization);
  }

  for (beta = 0; beta < NBETA; beta++)
    printf ("---------");
  printf ("\n");


  return 0;
}


