#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
#define SAMPLES 32


int
main (int argc, char **argv)
{
  clock_t t_start,t1,t2;

  char myfile[STR_LENG]={0};
  int rec_start=0;
  int rec_max = ITER_SWAP / ITER_SWAP_KERN;
  int rec_av=rec_max-rec_start;
  //printf("rec_max=%d\n",rec_max);
  double tmp_q;

  double q1[NBETA];
  double sx2[NBETA];
  double sx4[NBETA];

  double x2[NBETA][SAMPLES];
  double x4[NBETA][SAMPLES];

  double q2[NBETA][SAMPLES];
  double q3[NBETA][SAMPLES];
  double q4[NBETA][SAMPLES];

  
  double qk1_real[NBETA];
  double qk1_imag[NBETA];
  double qk2[NBETA][SAMPLES];
  double sxk2[NBETA];


  double qkk1_real[NBETA];
  double qkk1_imag[NBETA];
  double qkk2[NBETA][SAMPLES];
  double sxkk2[NBETA];


  double rj1[NBETA];
  double rj2[NBETA];

  double corrj1[NBETA];
  double corrj2[NBETA];
  //  double binderj[NBETA][SAMPLES];
  double binderj1[NBETA];
  double binderj2[NBETA];

  int pair, rec, beta;

  double C=1.0/(2*sin(PI/SZ))/SZ;

  t_start=clock();
  for(beta = 0 ; beta< NBETA; beta++){
    q1[beta]=0.0f;
    qk1_real[beta]=0.0f;
    qk1_imag[beta]=0.0f;

    qkk1_real[beta]=0.0f;
    qkk1_imag[beta]=0.0f;

    sx2[beta]=0.0f;
    sx4[beta]=0.0f;    
    binderj1[beta]=0.0f;
    binderj2[beta]=0.0f;
    sxk2[beta]=0.0f;
    sxkk2[beta]=0.0f;
    corrj1[beta]=0.0f;
    corrj2[beta]=0.0f;
    rj1[beta]=0.0f;
    rj2[beta]=0.0f;
  }


  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"q_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
            

    for(beta=0;beta<NBETA;beta++){
      q2[beta][pair]=0.0f;
      q3[beta][pair]=0.0f;
      q4[beta][pair]=0.0f;
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
	tmp_q/=SZ_CUBE;
        q1[beta]+=tmp_q;
	q2[beta][pair]+=tmp_q*tmp_q;
	q3[beta][pair]+=tmp_q*tmp_q*tmp_q;
	q4[beta][pair]+=tmp_q*tmp_q*tmp_q*tmp_q;
      }
      fscanf (fp, "\n");
    }
    fclose (fp);


    
    for(beta=0;beta<NBETA;beta++){
      q2[beta][pair]/=(rec_av);
      q3[beta][pair]/=(rec_av);
      q4[beta][pair]/=(rec_av);
    }
  } 
    

  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk1r_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
            
    for(beta=0;beta<NBETA;beta++){
      qk2[beta][pair]=0.0f;
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
	tmp_q/=SZ_CUBE;
        qk1_real[beta]+=tmp_q;
	qk2[beta][pair]+=tmp_q*tmp_q;
      }
      fscanf (fp, "\n");
    }
    fclose (fp);
    
    
  } 



  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk1i_pair_%05d.txt", pair);
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
	tmp_q/=SZ_CUBE;
        qk1_imag[beta]+=tmp_q;
	qk2[beta][pair]+=tmp_q*tmp_q;
      }
      fscanf (fp, "\n");
    }
    fclose (fp);
    
    for(beta=0;beta<NBETA;beta++){
      qk2[beta][pair]/=(rec_av);
    }
  } 


  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk2r_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    

    for(beta=0;beta<NBETA;beta++){
      qkk2[beta][pair]=0.0f;
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
	tmp_q/=SZ_CUBE;
        qkk1_real[beta]+=tmp_q;
	qkk2[beta][pair]+=tmp_q*tmp_q;
      }
      fscanf (fp, "\n");
    }
    fclose (fp);

    
    for(beta=0;beta<NBETA;beta++){
      qkk2[beta][pair]/=(rec_av);
    }
  } 
  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk2i_pair_%05d.txt", pair);
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
	tmp_q/=SZ_CUBE;
        qkk1_imag[beta]+=tmp_q;
	qkk2[beta][pair]+=tmp_q*tmp_q;
      }
      fscanf (fp, "\n");
    }
    fclose (fp);

    
    for(beta=0;beta<NBETA;beta++){
      qkk2[beta][pair]/=(rec_av);
    }
  } 

  t1=clock();
  
   
  for(beta=0;beta<NBETA;beta++){
      double aq1=q1[beta]/(SAMPLES*(rec_av));
      double aqk1=sqrt(qk1_real[beta]*qk1_real[beta]+qk1_imag[beta]*qk1_imag[beta])/(SAMPLES*(rec_av));
      double aqkk1=sqrt(qkk1_real[beta]*qkk1_real[beta]+qkk1_imag[beta]*qkk1_imag[beta])/(SAMPLES*(rec_av));
      //      double aqkk1=qkk1[beta]/(SAMPLES*(rec_av));

    for(pair=0;pair<SAMPLES;pair++){

      x2[beta][pair]=q2[beta][pair]-aq1*aq1;
      x4[beta][pair]=q4[beta][pair]-4.0*q3[beta][pair]*aq1+6.0*q2[beta][pair]*aq1*aq1-3.0*aq1*aq1*aq1*aq1;
      sx2[beta]+=x2[beta][pair];
      sx4[beta]+=x4[beta][pair];

      qk2[beta][pair]=qk2[beta][pair]-aqk1*aqk1;
      sxk2[beta]+=qk2[beta][pair];

      qkk2[beta][pair]=qkk2[beta][pair]-aqkk1*aqkk1;
      sxkk2[beta]+=qkk2[beta][pair];

    }
    
    for(pair=0;pair<SAMPLES;pair++){
      double binder=1.5-0.5*(SAMPLES-1)*(sx4[beta]-x4[beta][pair])/(sx2[beta]-x2[beta][pair])/(sx2[beta]-x2[beta][pair]);
      double corr=C*sqrt((sx2[beta]-x2[beta][pair])/(sxk2[beta]-qk2[beta][pair])-1.0);
      double r=(sxk2[beta]-qk2[beta][pair])/(sxkk2[beta]-qkk2[beta][pair]);
      binderj1[beta]+=binder;
      binderj2[beta]+=binder*binder;
      corrj1[beta]+=corr;
      corrj2[beta]+=corr*corr;
      rj1[beta]+=r;
      rj2[beta]+=r*r;
    }
    binderj2[beta]/=SAMPLES;
    binderj1[beta]/=SAMPLES;
    rj2[beta]/=SAMPLES;
    rj1[beta]/=SAMPLES;
    corrj1[beta]/=SAMPLES;
    corrj2[beta]/=SAMPLES;
    printf("%1.7f\t%1.7f\t",binderj1[beta],sqrt((binderj2[beta]-binderj1[beta]*binderj1[beta])*(SAMPLES-1)));
    printf("%1.7f\t%1.7f\t",corrj1[beta],sqrt((corrj2[beta]-corrj1[beta]*corrj1[beta])*(SAMPLES-1)));
    printf("%1.7f\t%1.7f\n",rj1[beta],sqrt((rj2[beta]-rj1[beta]*rj1[beta])*(SAMPLES-1)));
  }



  t2=clock();
  printf("t1=%f\nt2=%f\n",((float)t1-t_start)/CLOCKS_PER_SEC,((float)(t2-t1))/CLOCKS_PER_SEC);

  return 0;
}


