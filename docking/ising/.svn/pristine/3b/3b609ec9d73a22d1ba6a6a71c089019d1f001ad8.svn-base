#include <stdio.h>
#include <stdlib.h>

#include "../sim.h"


// analysis q[PAIR].txt


int
main (int argc, char **argv)
{
  char myfile[STR_LENG];

  int rec_max = ITER_SWAP / ITER_SWAP_KERN;

  int tmp_q;
  int q[GD/2][rec_max][NBETA];
  double q_avrg[GD/2][NBETA];
  int pair, rec, beta;



  /// read from file to q[GD/2][rec_max][NBETA]
  for (pair = 0; pair < GD/2; pair++) {

    sprintf (myfile, "pair_%04d_q.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%d", &tmp_q);
	q[pair][rec][beta] = tmp_q;
      }
      fscanf (fp, "\n");
    }

    fclose (fp);
  }




  /// compute summation

  for (pair = 0; pair < GD/2; pair++)
    for (beta = 0; beta < NBETA; beta++)
      q_avrg[pair][beta] = 0;

  for (pair = 0; pair < GD/2; pair++)
    for (rec = 0; rec < rec_max; rec++)
      for (beta = 0; beta < NBETA; beta++)
	q_avrg[pair][beta] += (double) abs(q[pair][rec][beta]);




  /// print

  for (beta = 0; beta < NBETA; beta++)
    printf ("   b%02d  ", beta);
  printf ("\n");

  for (beta = 0; beta < NBETA; beta++)
    printf ("--------");
  printf ("\n");

  for (pair = 0; pair < GD/2; pair++) {
    for (beta = 0; beta < NBETA; beta++) {
      q_avrg[pair][beta] /= rec_max;
      printf ("%7.1f\t",  q_avrg[pair][beta]);
    }
    //printf ("\n");
    printf ("\t pair %02d\n", pair);
  }

  for (beta = 0; beta < NBETA; beta++)
    printf ("--------");
  printf ("\n");


  return 0;
}


