#include "COPYING"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>

#include "sim.h"




void
host_timing_init (Timing * t, int n)
{
  for (int i = 0; i < n; i++) {
    t[i].start = 0.0;
    t[i].span = 0.0;
    t[i].n = 0;
  }
}


/*
  execution mode
  0 - marking start time stamp
  1 - marking end time stamp, and accumulated time span
  2 - culculate average
*/
void
host_timing (Timing * t, int idx, int mode)
{
  switch (mode) {
  case 0:
    t[idx].start = host_time_now ();
    break;
  case 1:
    t[idx].span += t[idx].start - host_time_now ();
    t[idx].n += 1;
    break;
  default:
    t[idx].avrg = t[idx].span / t[idx].n;
  }
}



double
host_time_now ()
{
  struct timeval mytime;
  gettimeofday (&mytime, NULL);
  double mytime_second =
    (double) mytime.tv_sec + (double) mytime.tv_usec / 1.0e6;

  return mytime_second;
}





/*
// convert lattice allocation, shared -> seperate
for (int run = 0; run < 2; run++) {
  for (z = 0; z < SZ; z++) {
    for (y = 0; y < SZ; y++) {
      for (x = 0; x < SZ_HF; x++) {
	idx_src = (SZ * SZ_HF * z) + (SZ_HF * y) + (x * 2 + run);
	idx_dst = (SZ * SZ_HF * z) + (SZ_HF * y) + (x + SZ_CUBE_HF * run);
	l1[idx_dst] = l0[idx_src];
      }
    }   
  }   
 }
*/




// initiate J
void
host_init_J (MSC_DATATYPE * l)
{
  MSC_DATATYPE jx[SZ_CUBE];
  MSC_DATATYPE jy[SZ_CUBE];
  MSC_DATATYPE jz[SZ_CUBE];

  MSC_DATATYPE j;
  int i;
  int x, xa, xb, y, ya, yb, z, za, zb;
  int idx;

#if DENSE == SPARSE
  MSC_DATATYPE mask_s = MASK_S;
#elif DENSE == COMPACT
  MSC_DATATYPE mask_s = MASK_S0;
#endif

  for (i = 0; i < SZ_CUBE; i++) {
    jx[i] = 1;
    jy[i] = 1;
    jz[i] = 1;
#ifdef RANDJ
    jx[i] = rand () & 1;
    jy[i] = rand () & 1;
    jz[i] = rand () & 1;
#endif
  }

  for (z = 0; z < SZz; z++) {
    za = (z + SZz - 1) % SZz;
    zb = z;			// zb != z + 1
    for (y = 0; y < SZ; y++) {
      ya = (y + SZ - 1) % SZ;
      yb = y;
      for (x = 0; x < SZ; x++) {
	xa = (x + SZ - 1) % SZ;
	xb = x;
	idx = SZ * SZ * z + SZ * y + x;

	j =			// integrated
	  MASK_J0 * jx[SZ * SZ * z + SZ * y + xa] |	// left
	  MASK_J1 * jx[SZ * SZ * z + SZ * y + xb] |	// right
	  MASK_J2 * jy[SZ * SZ * z + SZ * ya + x] |	// up
	  MASK_J3 * jy[SZ * SZ * z + SZ * yb + x] |	// down
	  MASK_J4 * jz[SZ * SZ * za + SZ * y + x] |	// front
	  MASK_J5 * jz[SZ * SZ * zb + SZ * y + x];	// back
	l[idx] = j | (l[idx] & mask_s);
      }
    }
  }

}


void
host_init_S (MSC_DATATYPE * l)
{

  MSC_DATATYPE s;

#if DENSE == SPARSE
  MSC_DATATYPE mask_s = MASK_S;
#elif DENSE == COMPACT
  MSC_DATATYPE mask_s = MASK_S0;
#endif

  for (int i = 0; i < SZ_CUBE; i++) {
    s = 1;
#ifdef RANDS
    s = rand ();
#endif
    s = s & mask_s;
    l[i] = (l[i] & MASK_J) | s;
  }
}



// logically, lattice[GD][SZ_CUBE]
void
host_init_lattice (MSC_DATATYPE * lattice)
{
  MSC_DATATYPE *l;			// local lattice
  l = (MSC_DATATYPE *) malloc (sizeof (MSC_DATATYPE) * SZ_CUBE);
  for (int i = 0; i < SZ_CUBE; i++) {
    l[i] = 0;
  }

  for (int cube = 0; cube < GD; cube++) {
    host_init_S (l);
    if ((cube & 1) == 0)	// common J for adjacent realization
      host_init_J (l);

    for (int word = 0; word < NWORD; word++) {
      int offset = (NWORD * SZ_CUBE * cube) + (SZ_CUBE * word);
      for (int i = 0; i < SZ_CUBE; i++)
	lattice[offset + i] = l[i];
    }
  }

  free (l);
}








/*
  output format of
  e[REALIZATION].txt
  q[REALIZATION_PAIR].txt

  X axis: beta, from low to high
  Y axis: iteration space
*/


void
host_save_st (St * st, char *mydir, int node, int device)
{
  //printf ("node%03d device%d write outputs to %s/\n", node, device, mydir);

  const int rec_max = ITER_SWAP / ITER_SWAP_KERN;

  /*
     for (int realization = 0; realization < GD; realization++) {
     char myfile[STR_LENG];
     sprintf (myfile, "%s/node%03d_dev%d_e%04d.txt", mydir, node, device, realization);
     FILE *fp = fopen (myfile, "a");
     if (fp == NULL) {
     fprintf (stderr, "failed to open %s\n", myfile);
     exit (1);
     }

     for (int rec = 0; rec < rec_max; rec++) {
     for (int b = 0; b < NBETA; b++)
     fprintf (fp, "%8.2f ", st[rec].e[realization][b]);
     fprintf (fp, "\n");
     }

     fclose (fp);

     }
   */

  //output q(0)
  for (int pair = 0; pair < GD / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/q_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].q[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }

  //output q(k1)
  for (int pair = 0; pair < GD / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/qk1r_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].qk_real[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }

  for (int pair = 0; pair < GD / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/qk1i_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].qk_imag[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }



  //output q(k2)
  for (int pair = 0; pair < GD / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/qk2r_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].qk2_real[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }
  for (int pair = 0; pair < GD / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/qk2i_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].qk2_imag[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }

}







// initialize temperatures
void
host_init_temp (Temp * temp, float beta_low, float beta_high)
{
  double a = beta_low;
  Temp *tmp = (Temp *) malloc (sizeof (Temp) * NBETA_MAX);
  const double beta_ratio = exp (log (beta_high / beta_low) / (NBETA - 1));
  //const float beta_delta = (beta_high - beta_low) / NBETA;

  // generating temperatures for one lattice
  for (int b = 0; b < NBETA; b++) {
    tmp[b].idx = b;
    tmp[b].beta = a;
    //a += beta_delta;
    a *= beta_ratio;
  }

#if 0
  printf("List of Betas used in this run:\n");
  for (int b = 0; b < NBETA; b++)
    printf("%f\n",b,tmp[b].beta);
#endif

  for (int b = NBETA; b < NBETA_MAX; b++) {
    tmp[b].idx = b;
    tmp[b].beta = 0;
  }

  // duplication for all lattices
  for (int lattice = 0; lattice < GD; lattice++) {
    for (int b = 0; b < NBETA_MAX; b++) {
      int i = NBETA_MAX * lattice + b;
      //temp_idx[i] = tmp[b].idx;
      //temp_beta[i] = tmp[b].beta;
      temp[i] = tmp[b];
    }
  }

  free (tmp);
}




void
host_report_speed_title ()
{
  printf ("\t\t\t\t     second   iter/s    spin/s    s/spin    ");
  printf ("ps/spin   ");
  //printf ("SM_R/W    ");
  //printf ("GFLOPS  ");
  putchar ('\n');
}



void
host_report_speed (double start, double stop, int iter, char *event)
{
  double s = stop - start;
  float iter_per_s = iter / s;
  float spin_per_s = NBETA * SZ_CUBE * GD * iter_per_s;
  float s_per_spin = 1 / spin_per_s;
  //int sm_w_bw = NWORD * sizeof (MSC_DATATYPE) * SZ_CUBE * iter_per_s / 1048576;
  //int sm_r_bw = sm_w_bw * 7;  // average per TB, GB/s

  printf ("%s   %9.6f, ", event, s);
  printf ("%.3e %.3e %.3e ", iter_per_s, spin_per_s, s_per_spin);
  printf ("%06.3f    ", s_per_spin * 1.0e12);
  //printf ("%04d/%04d ", sm_r_bw, sm_w_bw);
  //printf ("%s", "???");
  putchar ('\n');
}





void
host_usage (char *bin)
{
  fprintf (stderr, "usage: %s [options]\n", bin);
  fprintf (stderr, " -l <beta_lower_bound>     default %f\n", BETA_LOW);
  fprintf (stderr, " -u <beta_upper_bound>     default %f\n", BETA_HIGH);
  fprintf (stderr, " where lower bound must less than upper bound\n");
  exit (1);
}



void
host_summary (float beta_low, float beta_high, char *mydir)
{

  char prob_gen[2][16] = {"DIRECT_COMPUTE", "TABLE_LOOKUP"};
  char allocation[3][16] = {"SHARED", "SEPARATED", "INTEGRATED"};
  char dense[2][16] = {"SPARSE", "COMPACT"};

  printf ("\n");
  printf ("lattice size:\t\t %d*%d*%d\n", SZ, SZ, SZz);
  printf ("realizations per GPU:\t %d\n", GD);
  printf ("parallel betas:\t\t %.5f-%.5f,", beta_low, beta_high);
  printf (" %d samples\n", NBETA);
  printf ("external field:\t\t %.5f\n", H);
  printf ("ITER_WARMUP:\t\t %d\n", ITER_WARMUP);
  printf ("ITER_WARMUP_KERN:\t %d\n", ITER_WARMUP_KERN);
  printf ("ITER_WARMUP_KERNFUNC:\t %d\n", ITER_WARMUP_KERNFUNC);
  printf ("ITER_SWAP:\t\t %d\n", ITER_SWAP);
  printf ("ITER_SWAP_KERN:\t\t %d\n", ITER_SWAP_KERN);
  printf ("ITER_SWAP_KERNFUNC:\t %d\n", ITER_SWAP_KERNFUNC);
  printf ("output dir:\t\t %s\n", mydir);
  printf ("\n");

  printf ("PROB_GEN:\t\t %s\n", prob_gen[PROB_GEN]);
  printf ("ALLOCATION:\t\t %s\n", allocation[ALLOCATION]);
  printf ("DENSE:\t\t\t %s\n", dense[DENSE]);
  printf ("MSCT:\t\t\t %d-%d\n", MSCT, MSC_FORMAT);
  printf ("threads per block:\t %d\n", TperB);
  printf ("blocks per GPU:\t\t %d\n", GD);
  printf ("\n");
}



void
host_makedir (char *mydir)
{
  struct stat st;

  // there is no such a dir
  if (stat (mydir, &st) != 0) {
    printf ("mkdir %s\n", mydir);
    mkdir (mydir, 0777);
  }

  // the dir exists
  else { 
	printf ("program terminated: output dir conflict\n");
    exit (1);
  }
/*
{
    printf ("the output directory %s exists, force overwriting? y/n ", mydir);
    char input;
    scanf ("%c", &input);
    if (input == 'n') {
      printf ("program terminated\n");
      exit (1);
    }
  }
*/
}
