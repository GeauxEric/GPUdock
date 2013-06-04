#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sim.h"




//#define RUNS 80
//#define SAMPLES (GD*RUNS/2*3)
#define SAMPLES 51840

#define rec_start 1000
#define rec_max (ITER_SWAP / ITER_SWAP_KERN - rec_start )


int
main (int argc, char **argv)
{
  char myfile[STR_LENG]={0};
  double tmp_qk;
  int pair, rec, beta;

   
  double r12_1[SAMPLES][NBETA];
  double r12_2[SAMPLES][NBETA];
  double r12[SAMPLES];

  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk1_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    
    for(beta=0;beta<NBETA;beta++)
      r12_1[pair][beta]=0;


    for(rec=0;rec<rec_start;rec++){
      for(beta=0;beta<NBETA;beta++)
	fscanf (fp, "%lf\t", &tmp_qk);
      fscanf (fp, "\n");
    }

    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf\t", &tmp_qk);
	tmp_qk=tmp_qk*tmp_qk/(SZ_CUBE*SZ_CUBE)/(rec_max);//[<qk^2>];
	//	qk_avrg[beta]+=tmp_qk;
	r12_1[pair][beta]+=tmp_qk;
      }
      

      fscanf (fp, "\n");
    }
    
    fclose (fp);
    //printf("\rfinished reading qk1_pair_%05d.txt",pair);
  }

  //printf("\n");
  for (pair = 0; pair < SAMPLES; pair++) {
    if(myfile!=NULL)
      int charcheck = sprintf ( myfile ,"qk2_pair_%05d.txt", pair);
    FILE *fp = fopen (myfile, "r");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

        
    for(beta=0;beta<NBETA;beta++)
      r12_2[pair][beta]=0;

    
    for(rec=0;rec<rec_start;rec++){
      for(beta=0;beta<NBETA;beta++)
	fscanf (fp, "%lf\t", &tmp_qk);
      fscanf (fp, "\n");
    }


    for (rec = 0; rec < rec_max; rec++) {
      for (beta = 0; beta < NBETA; beta++) {
	fscanf (fp, "%lf\t", &tmp_qk);
	tmp_qk=tmp_qk*tmp_qk/(SZ_CUBE*SZ_CUBE)/(rec_max);//[<q^2>]
	//qk2_avrg[beta]+=tmp_qk;
	r12_2[pair][beta]+=tmp_qk;
      }
      fscanf (fp, "\n");
    }
    
    fclose (fp);
    //printf("\rfinished reading qk2_pair_%05d.txt",pair);
  }


  /*calculating r12 with jack knife*/
  printf("\nR12:\n");
  for(beta=0;beta<NBETA;beta++){
    double sum_1=0;
    double sum_2=0;
    double error_1=0;
    double error_2=0;
    for(pair=0;pair<SAMPLES;pair++){
	sum_1+=r12_1[pair][beta];
	sum_2+=r12_2[pair][beta];
    }
    for(pair=0;pair<SAMPLES;pair++){
      double temp_r=(sum_1-r12_1[pair][beta])/(SAMPLES-1);
      error_1+=temp_r*temp_r;
      temp_r=(sum_2-r12_2[pair][beta])/(SAMPLES-1);
      error_2+=temp_r*temp_r;
      //	r12[pair]=(sum_1-r12_1[pair][beta])/(sum_2-r12_2[pair][beta]);
    }
    sum_1/=SAMPLES;
    //    error_1=sqrt((error_1/SAMPLES-sum_1*sum_1)/(SAMPLES-1))/sum_1;
    error_1=sqrt(((double)SAMPLES-1)*(error_1/SAMPLES-sum_1*sum_1))/sum_1;
    sum_2/=SAMPLES;
    //    error_2=sqrt((error_2/SAMPLES-sum_2*sum_2)/(SAMPLES-1))/sum_2;
    error_2=sqrt(((double)SAMPLES-1)*(error_2/SAMPLES-sum_2*sum_2))/sum_2;


    //printf("%f,%f\n",sum_1,sum_2);
    //printf("%f,%f\n",error_1,error_2);
    double avrg=0;
    double error=0;
    /* for(pair=0;pair<SAMPLES;pair++){
      avrg+=r12[pair];
      error+=r12[pair]*r12[pair];
    }
    */
    avrg=sum_1/sum_2;
    error=avrg*sqrt(error_1*error_1+error_2*error_2);
    printf("%f\t%f\n",avrg,error);
    
  }

  return 0;
}


