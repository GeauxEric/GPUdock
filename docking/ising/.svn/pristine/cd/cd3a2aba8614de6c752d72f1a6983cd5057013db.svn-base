#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sim.h"


// analysis realization_[REALIZATION]_e.txt

#define RUNS 3
#define SAMPLES (GD*RUNS/2)
//#define SAMPLES 10
#define rec_start 100
#define NBINS 6
#define rec_max (ITER_SWAP / ITER_SWAP_KERN - rec_start)
#define BIN_SIZE (SAMPLES / NBINS)



int
main (int argc, char **argv)
{
  char myfile[STR_LENG]={0};
  //  int rec_start=5000;
  //  int rec_max = ITER_SWAP / ITER_SWAP_KERN - rec_start;
  //  int rec_max=5;
  
  //  printf("rec_max=%d\n",rec_max);
  int tmp_q;
  double qq;
  double tmp_qk;
  //int e[GD][2000][NBETA];
  int *q;
  int bin;

   q=(int *)malloc(sizeof(int)*SAMPLES*rec_max*NBETA);

   if(q==NULL){
     printf("allocation failed");
    return(0);
   }
   //free(e);

   
  double q2_avrg[NBINS][NBETA];
  double qk_avrg[NBINS][NBETA];
  double qk2_avrg[NBINS][NBETA];
  double q4_avrg[NBINS][NBETA];
  double q22_avrg[NBINS][NBETA];
  double q22[NBINS][NBETA];
  int pair, rec, beta;
  //  double q2[BINS],q4[BINS];

  double r12_1[SAMPLES][NBETA];
  double r12_2[SAMPLES][NBETA];
  
  /// read from file to e[GD][rec_max][NBETA]
  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"q_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for(rec=0; rec<rec_start;rec++)
      for (beta = 0; beta < NBETA; beta++) 
	fscanf (fp, "%d", &tmp_q);//data not to be used, simply throw away~      
     if(pair==0)
      printf("%d steps of warm up data thrown away\n",rec);

    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%d", &tmp_q);
	q[pair*rec_max*NBETA+rec*NBETA+beta]=tmp_q;
      }
      fscanf (fp, "\n");
    }
    if(pair==0)
      printf("%d steps of data read in for calculation\n",rec);

    fclose (fp);
    printf("\r finished reading q_pair_%05d.txt",pair);
  }

  printf("\n");
  for(bin=0; bin<NBINS; bin++)
    for (beta = 0; beta < NBETA; beta++)
      qk_avrg[bin][beta]=0;


  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk1_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    
    for (rec = 0; rec < rec_start; rec++) 
      for (beta = 0; beta < NBETA; beta++) 
      	fscanf (fp, "%lf\t", &tmp_qk);

    bin= pair / BIN_SIZE;
    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf\t", &tmp_qk);
	qk_avrg[bin][beta]+=tmp_qk*tmp_qk/(SZ_CUBE*SZ_CUBE)/(BIN_SIZE)/(SAMPLES);//[<qk^2>]
      }
      fscanf (fp, "\n");
    }
    
    fclose (fp);
    printf("\r finished reading qk1_pair_%05d.txt",pair);
  }


  printf("\n");
    /*read q(k2) */
  for(bin=0; bin<NBINS; bin++)
    for (beta = 0; beta < NBETA; beta++)
      qk2_avrg[bin][beta]=0;


  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk2_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    
    for (rec = 0; rec < rec_start; rec++) 
      for (beta = 0; beta < NBETA; beta++) 
      	fscanf (fp, "%lf\t", &tmp_qk);

    bin= pair / BIN_SIZE;
    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf\t", &tmp_qk);
	qk2_avrg[bin][beta]+=tmp_qk*tmp_qk/(SZ_CUBE*SZ_CUBE)/(BIN_SIZE)/(SAMPLES);//[<qk^2>]
      }
      fscanf (fp, "\n");
    }
    
    fclose (fp);
    printf("\r finished reading qk2_pair_%05d.txt",pair);
  }



  printf("\n Computing...\n");
    
  /// compute summation
  
  //  for (rec = 0; rec < rec_max; rec++)
  for(bin=0;bin<NBINS;bin++){
    for (beta = 0; beta < NBETA; beta++){
      q2_avrg[bin][beta] = 0;
      q4_avrg[bin][beta] = 0;
      q22_avrg[bin][beta] = 0;
      q22[bin][beta] = 0;
    }
  }    
  for (pair = 0; pair < (SAMPLES); pair++){
    bin=pair/BIN_SIZE;
    for (rec = 0; rec < rec_max; rec++){
      for (beta = 0; beta < NBETA; beta++){
	//	  printf("\rbin:%d,pair:%d,rec:%d,beta:%d",nbin,pair,rec,beta);
	//	q_avrg[rec][beta] += (double) q[realization*rec_max*NBETA+rec*NBETA+beta];
	tmp_q= q[pair*rec_max*NBETA+rec*NBETA+beta];
	qq=(double)tmp_q/SZ_CUBE;
	q2_avrg[bin][beta]+=qq*qq/(BIN_SIZE)/rec_max;//[<q^2>]
	q4_avrg[bin][beta]+=qq*qq*qq*qq/(BIN_SIZE)/rec_max;//[<q^4>]
	q22[bin][beta]+=qq*qq/(rec_max);//<q^2>
      }
    }
    for(beta=0;beta<NBETA;beta++){
      q22_avrg[bin][beta]+=q22[bin][beta]*q22[bin][beta]/BIN_SIZE;//[<q^2>^2]
      q22[bin][beta]=0;
    }
  }
  //free(q);
    
    
 
  
  
  free(q);

  double measure[5][NBINS];
  /*
  double cl[BINS];
  double r12[BINS];
  double binder[BINS];
  double A[BINS];
  double G[BINS];
  */

  
  double average[5][NBETA];
  double error[5][NBETA];
  
  for(beta=0;beta<NBETA;beta++){
    for(bin=0;bin<NBINS;bin++){
      measure[0][bin]=0.5f/sin(3.141592654f/SZ)*sqrt(q2_avrg[bin][beta]/qk_avrg[bin][beta]-1)/SZ;
      measure[1][bin]=qk_avrg[bin][beta]/qk2_avrg[bin][beta];
      
      
      //      for(beta=0;beta<NBETA;beta++)
      q2_avrg[bin][beta]=q2_avrg[bin][beta]*q2_avrg[bin][beta];//[<q^2>]^2
      
      
      
      measure[2][bin]=  0.5f*(3-q4_avrg[bin][beta]/(q2_avrg[bin][beta]));
      measure[3][bin]=(q22_avrg[bin][beta]-q2_avrg[bin][beta])/q2_avrg[bin][beta];
      measure[4][bin]=(q22_avrg[bin][beta]-q2_avrg[bin][beta])/(q4_avrg[bin][beta]-q2_avrg[bin][beta]);
    }

    
     
    for(int n=0;n<5;n++){
      average[n][beta]=0;
      error[n][beta]=0;
      for(bin=0;bin<NBINS;bin++){
	average[n][beta]+=measure[n][bin];
	error[n][beta]+=measure[n][bin]*measure[n][bin];
      }
     average[n][beta]/=NBINS;
     error[n][beta]/=NBINS;
     error[n][beta]=sqrt((error[n][beta]-average[n][beta]*average[n][beta])/(NBINS-1));
    }

  }

  printf("Correlation Length/L:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",average[0][beta]);
  printf("\n");
  
  
  printf("R12:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",average[1][beta]);
  printf("\n");
  
  
  printf("B:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",average[2][beta]);
  printf("\n");
  
    
   
  printf("A:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",average[3][beta]);
  printf("\n");
  
  printf("G:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",average[4][beta]);
  printf("\n");
  
  printf("error analysis:\n");
  printf("Correlation Length/L:\n");
    for(beta=0;beta<NBETA;beta++)
      printf("%f\n",error[0][beta]);
    printf("\n");
    
    
    printf("R12:\n");
    for(beta=0;beta<NBETA;beta++)
      printf("%f\n",error[1][beta]);
    printf("\n");
    
    
    printf("B:\n");
    for(beta=0;beta<NBETA;beta++)
      printf("%f\n",error[2][beta]);
    printf("\n");
    
    
   
  printf("A:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",error[3][beta]);
  printf("\n");
  
  printf("G:\n");
  for(beta=0;beta<NBETA;beta++)
    printf("%f\n",error[4][beta]);
  printf("\n");
  


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


