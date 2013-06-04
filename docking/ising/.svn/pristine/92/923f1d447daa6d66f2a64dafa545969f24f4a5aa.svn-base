// nvcc -O2 -gencode arch=compute_20,code=sm_20 -o rand rand.cu lcg.cu
// nvcc -O2 -gencode arch=compute_20,code=sm_20 -o rand rand.cu


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <sys/time.h>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include "lcg.h"




// number of random number to generate per thread
#define NRAND 1000000

#define GD 64;
#define BDx 8;
#define BDy 16;
#define BDz 16;
#define TperB 256;

#define PROB_DATATYPE float
//#define PROB_DATATYPE u_int32_t



__constant__ int seed;
__constant__ curandState *pgpuseed;
__constant__ LCG_DATATYPE *pgpuseed1;
__constant__ int *pcnt;
__constant__ PROB_DATATYPE *pa;


/*
// initilize seeds for curand
// CURAND_Library.pdf, pp21
__global__ void
kernel_init_seed ()
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;
  const int gidx = TperB * blockIdx.x + bidx;
  curand_init (seed, gidx, 0, &pgpuseed[gidx]);
}
*/






#if 0
__device__ void
kernel0 ()
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;

  curandState seed0 = pgpuseed[TperB * blockIdx.x + bidx];       // curand sequence
  //LCG_DATATYPE seed1 = pgpuseed1[TperB * blockIdx.x + bidx];    // LCG PRNG sequence

  /*
  int rand123_cnt = pcnt[TperB * blockIdx.x + bidx];
  philox4x32_key_t rand123_k = {{bidx, seed1}};
  philox4x32_ctr_t rand123_c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};
  union {
    philox4x32_ctr_t c;
    uint4 i;
  } rand123_u;
  */


  for (int i = 0; i < NRAND; i++) {
    PROB_DATATYPE myrand = curand (&seed0);
  }

  pgpuseed[TperB * blockIdx.x + bidx] = seed0;
  //pgpuseed1[TperB * blockIdx.x + bidx] = seed1;
  //pcnt[TperB * blockIdx.x + bidx] = rand123_cnt;
}
#endif








double
host_time_now ()
{
  struct timeval mytime;
  gettimeofday (&mytime, NULL);
  double mytime_second =
    (double) mytime.tv_sec + (double) mytime.tv_usec / 1.0e6;

  return mytime_second;
}



void
launcher (int ver)
{
  int dim_grid = GD;
  dim3 dim_block (BDx, BDy, BDz);

  printf ("hello\n");

#if 0
  switch (ver) {
  case 0:
    //kernel0 <<< dim_grid, dim_block >>> ();
    printf ("hello\n");
    break;
    /*
  case 1:
    kernel1 <<< dim_grid, dim_block >>> ();
    break;
  case 2:
    kernel2 <<< dim_grid, dim_block >>> ();
    break;
  case 3:
    kernel3 <<< dim_grid, dim_block >>> ();
    break;
    */
  default:
    printf ("branch %d not available\n", ver);
    break;
  }
#endif
}


void
report_speed (double start, double stop, int ver)
{
  double s = stop - start;
  long nrand = NRAND * TperB * GD;
  float rand_per_s = (float) nrand / s;
  float ps_per_rand = 1.0e12 / rand_per_s;
  printf ("PRNG version %d: \t%10.5f rand/s \t%10.5f ps/rand\n",
	  rand_per_s,
	  ps_per_rand);
}



int
main (int argc, char **argv)
{
  srand (time (NULL));

  int dim_grid = GD;
  dim3 dim_block (BDx, BDy, BDz);



  // gpuseed - curand seed
  int myrand = rand ();
  cudaMemcpyToSymbol ("seed", &myrand, sizeof (int), 0, cudaMemcpyHostToDevice);

  curandState *gpuseed_dev;
  size_t gpuseed_sz = sizeof (curandState) * TperB * GD;
  cudaMalloc ((void **) &gpuseed_dev, gpuseed_sz);
  cudaMemcpyToSymbol ("pgpuseed", &gpuseed_dev, sizeof (curandState *), 0, cudaMemcpyHostToDevice);
  //kernel_init_seed <<< dim_grid, dim_block >>> ();
  cudaThreadSynchronize ();


  // seeds for LCG PRNG
  LCG_DATATYPE *gpuseed1, *gpuseed1_dev;
  size_t gpuseed1_sz = sizeof (LCG_DATATYPE) * TperB * GD;
  gpuseed1 = (LCG_DATATYPE *) malloc (gpuseed1_sz);
  cudaMalloc ((void **) &gpuseed1_dev, gpuseed1_sz);
  cudaMemcpyToSymbol ("pgpuseed1", &gpuseed1_dev, sizeof (LCG_DATATYPE *), 0,
                      cudaMemcpyHostToDevice);


  // cnt - counter for rand123 PRNG
  int *cnt, *cnt_dev;
  size_t cnt_sz = sizeof (int) * TperB * GD;
  cnt = (int *) malloc (cnt_sz);
  cudaMalloc ((void **) &cnt_dev, cnt_sz);
  cudaMemcpyToSymbol ("pcnt", &cnt_dev, sizeof (int *), 0, cudaMemcpyHostToDevice);
  for (int i = 0; i < TperB * GD; i++)
    cnt[i] = 0;
  cudaMemcpy (cnt_dev, cnt, cnt_sz, cudaMemcpyHostToDevice);




  PROB_DATATYPE a_dev;
  size_t a_sz = sizeof (PROB_DATATYPE) * TperB * GD;
  cudaMalloc ((void **) &a_dev, a_sz);
  cudaMemcpyToSymbol ("pa", &a_dev, sizeof (PROB_DATATYPE *), 0, cudaMemcpyHostToDevice);




  for (int i = 0; i < 4; i++) {
    double time_start = host_time_now ();
    launcher (i);
    cudaThreadSynchronize ();
    double time_stop = host_time_now ();
    report_speed (time_start, time_stop, i);
  }


  cudaFree (gpuseed_dev);
  cudaFree (gpuseed1_dev);
  cudaFree (cnt_dev);
  cudaFree (a_dev);
  free (gpuseed1);
  free (cnt);

  return 0;
}



