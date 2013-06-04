#include "COPYING"


#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include "cuda_util.cuh"
#include "sim.h"
#include "sim.cuh"
#include "lcg.h"
#include "host_kernel.h"


extern texture <PROB_DATATYPE, 2, cudaReadModeElementType> mytexture;


void
host_launcher (float beta_low, float beta_high, char *mydir, int node,
	       int device)
{
  // select a device
  cudaSetDevice (device);

  // configure the GPU SRAM
  //cudaFuncSetCacheConfig (kernel_unified, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig (kernel_warmup, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig (kernel_swap, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig (kernel_rearrange, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig (kernel_compute_q, cudaFuncCachePreferShared);


  const int dim_grid0 = GD;
  const dim3 dim_block0 (BDx0, BDy0, BDz0);

  const int dim_grid3 = GD;
  const int dim_block3 = TperB;

  const int dim_grid4 = GD_HF;
  const int dim_block4 = TperB;


  // initilize random sequence
  srand (time (NULL) + 10 * node + device);


  // gpuseed - curand seed
  int myrand = rand ();
  cudaMemcpyToSymbol ("seed", &myrand, sizeof (int), 0,
		      cudaMemcpyHostToDevice);

  curandState *gpuseed_dev;
  size_t gpuseed_sz = sizeof (curandState) * TperB * GD;
  cudaMalloc ((void **) &gpuseed_dev, gpuseed_sz);
  cudaMemcpyToSymbol ("pgpuseed", &gpuseed_dev, sizeof (curandState *), 0,
		      cudaMemcpyHostToDevice);

  // how often should we re-initialize gpuseed???
  kernel_init_seed <<< dim_grid3, dim_block3 >>> ();
  cudaThreadSynchronize ();



  // seeds for LCG PRNG
  LCG_DATATYPE *gpuseed1, *gpuseed1_dev;
  size_t gpuseed1_sz = sizeof (LCG_DATATYPE) * TperB * GD;
  gpuseed1 = (LCG_DATATYPE *) malloc (gpuseed1_sz);
  cudaMalloc ((void **) &gpuseed1_dev, gpuseed1_sz);
  cudaMemcpyToSymbol ("pgpuseed1", &gpuseed1_dev, sizeof (LCG_DATATYPE *), 0, cudaMemcpyHostToDevice);

  host_init_seed (gpuseed1);
  CudaSafeCall (cudaMemcpy (gpuseed1_dev, gpuseed1, gpuseed1_sz, cudaMemcpyHostToDevice));



  // cnt - counter for rand123 PRNG
  // ziter * 2 * NBETA * ITER ~= 540 * ITER ~= 1 * 10^9
  int *cnt, *cnt_dev;
  size_t cnt_sz = sizeof (int) * TperB * GD;
  cnt = (int *) malloc (cnt_sz);
  cudaMalloc ((void **) &cnt_dev, cnt_sz);
  cudaMemcpyToSymbol ("pcnt", &cnt_dev, sizeof (int *), 0,
		      cudaMemcpyHostToDevice);
  for (int i = 0; i < TperB * GD; i++)
    cnt[i] = 0;
  cudaMemcpy (cnt_dev, cnt, cnt_sz, cudaMemcpyHostToDevice);



  // lattice
  MSC_DATATYPE *lattice, *lattice_dev;
  size_t lattice_sz = sizeof (MSC_DATATYPE) * SZ_CUBE * NWORD * GD;
  lattice = (MSC_DATATYPE *) malloc (lattice_sz);
  cudaMalloc ((void **) &lattice_dev, lattice_sz);
  host_init_lattice (lattice);
  cudaMemcpy (lattice_dev, lattice, lattice_sz, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol ("plattice", &lattice_dev, sizeof (MSC_DATATYPE *), 0,
		      cudaMemcpyHostToDevice);

  // lattice1
  // spins have been rearranged to reflect the temperature order
  MSC_DATATYPE *lattice1_dev;
  size_t lattice1_sz = sizeof (MSC_DATATYPE) * SZ_CUBE * GD;
  cudaMalloc ((void **) &lattice1_dev, lattice1_sz);
  cudaMemcpyToSymbol ("plattice1", &lattice1_dev, sizeof (MSC_DATATYPE *), 0,
		      cudaMemcpyHostToDevice);

  // temp - index and beta
  Temp *temp, *temp_dev;
  size_t temp_sz = sizeof (Temp) * NBETA_MAX * GD;
  temp = (Temp *) malloc (temp_sz);
  cudaMalloc ((void **) &temp_dev, temp_sz);
  host_init_temp (temp, beta_low, beta_high);
  cudaMemcpy (temp_dev, temp, temp_sz, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol ("ptemp", &temp_dev, sizeof (Temp *), 0,
		      cudaMemcpyHostToDevice);

#if 0
  // testing 1D texture
  PROB_DATATYPE *temp_texture, *temp_texture_dev;
  size_t temp_texture_sz = sizeof (float) * NBETA_MAX;
  temp_texture = (PROB_DATATYPE *) malloc (temp_texture_sz);
  for (int i = 0; i < NBETA_MAX; i++)
    temp_texture[i] = (PROB_DATATYPE) i;
  cudaMalloc ((void **) &temp_texture_dev, temp_texture_sz);
  cudaMemcpy (temp_texture_dev, temp_texture, temp_texture_sz,
	      cudaMemcpyHostToDevice);
  cudaBindTexture (NULL, mytexture, temp_texture_dev, temp_texture_sz);
#endif

#if 1
  // testing 2D texture
  PROB_DATATYPE *temp_texture, *temp_texture_dev;
  size_t temp_texture_sz = sizeof (float) * NBETA_MAX * NPROB_MAX;
  temp_texture = (PROB_DATATYPE *) malloc (temp_texture_sz);
  for (int i = 0; i < NBETA_MAX * NPROB_MAX; i++)
    temp_texture[i] = (PROB_DATATYPE) i;
  cudaMalloc ((void **) &temp_texture_dev, temp_texture_sz);
  cudaMemcpy (temp_texture_dev, temp_texture, temp_texture_sz,
	      cudaMemcpyHostToDevice);

  cudaChannelFormatDesc desc = cudaCreateChannelDesc <PROB_DATATYPE> ();
  cudaBindTexture2D(NULL, mytexture, temp_texture_dev, desc,
		    NPROB_MAX, NBETA_MAX, sizeof(PROB_DATATYPE) * NPROB_MAX);
#endif


  // st - status records
  St *st, *st_dev;
  size_t st_sz = sizeof (St) * ITER_SWAP / ITER_SWAP_KERN;
  st = (St *) malloc (st_sz);
  cudaMalloc ((void **) &st_dev, st_sz);
  cudaMemcpyToSymbol ("pst", &st_dev, sizeof (St *), 0,
		      cudaMemcpyHostToDevice);



#ifdef DEBUG0
  printf ("st_sz = %f MB\n", (float) st_sz / 1024 / 1024);
#endif



  char event[STR_LENG];

  //Timing t[4];
  //host_timing_init (t, 4);

  double t[4][2];
  double t2[3];
  double t3 = 0;

  putchar ('\n');
  host_report_speed_title ();




#if 1
  /// GPU ising

  // warm up runs
  t2[0] = host_time_now ();
  for (int i = 0; i < ITER_WARMUP; i += ITER_WARMUP_KERN) {
    t[0][0] = host_time_now ();

    kernel_warmup <<< dim_grid0, dim_block0 >>> ();
    //kernel_unified <<< dim_grid0, dim_block0 >>> (0, 0);
    CudaCheckError ();
    cudaThreadSynchronize ();

    t[0][1] = host_time_now ();
    sprintf (event, "n%03d d%d warmup %8d/%08d", node, device, i, ITER_WARMUP);
    host_report_speed (t[0][0], t[0][1], ITER_WARMUP_KERN, event);
  }

  t2[1] = host_time_now ();

  // swap runs
  for (int i = 0; i < ITER_SWAP; i += ITER_SWAP_KERN) {
    t[1][0] = host_time_now ();

    kernel_swap <<< dim_grid0, dim_block0 >>> (i / ITER_SWAP_KERN);
    //kernel_unified <<< dim_grid0, dim_block0 >>> (i / ITER_SWAP_KERN, 1);
    cudaThreadSynchronize ();
    t[1][1] = host_time_now ();
    t3 += t[1][1] - t[1][0];

    kernel_rearrange <<< dim_grid3, dim_block3 >>> ();
    CudaCheckError ();
    kernel_compute_q <<< dim_grid4, dim_block4 >>> (i / ITER_SWAP_KERN);
    CudaCheckError ();
    cudaThreadSynchronize ();
    t[2][1] = host_time_now ();

    sprintf (event, "n%03d d%d PT     %8d/%08d", node, device, i, ITER_SWAP);
    host_report_speed (t[1][0], t[1][1], ITER_SWAP_KERN, event);
  }
  t2[2] = host_time_now ();
#endif

#if 0
  /// CPU ising

  // warm up runs
  t2[0] = host_time_now ();
  for (int i = 0; i < ITER_WARMUP; i += ITER_WARMUP_KERN) {
    t[0][0] = host_time_now ();
    host_kernel_warmup (lattice, temp);
    t[0][1] = host_time_now ();
    sprintf (event, "n%03d d%d warmup %8d/%08d", node, device, i, ITER_WARMUP);
    host_report_speed (t[0][0], t[0][1], ITER_WARMUP_KERN, event);
  }
  
  t2[1] = host_time_now ();
  
  // swap runs
  for (int i = 0; i < ITER_SWAP; i += ITER_SWAP_KERN) {
    t[1][0] = host_time_now ();
    host_kernel_swap (lattice, temp, st, i / ITER_SWAP_KERN);
    t[1][1] = host_time_now ();
    t3 += t[1][1] - t[1][0];
    sprintf (event, "n%03d d%d PT     %8d/%08d", node, device, i, ITER_SWAP);
    host_report_speed (t[1][0], t[1][1], ITER_SWAP_KERN, event);
  }
  t2[1] = host_time_now ();
#endif




#ifndef NO_OUTPUT
  CudaSafeCall (cudaMemcpy (st, st_dev, st_sz, cudaMemcpyDeviceToHost));
  host_save_st (st, mydir, node, device);
#endif



  // report overall speed
  putchar ('\n');
  sprintf (event, "n%03d d%d overall warmup          ", node, device);
  host_report_speed (t2[0], t2[1], ITER_WARMUP, event);
  sprintf (event, "n%03d d%d overall PT (no measure) ", node, device);
  host_report_speed (0, t3, ITER_SWAP, event);
  sprintf (event, "n%03d d%d overall PT              ", node, device);
  host_report_speed (t2[1], t2[2], ITER_SWAP, event);
  sprintf (event, "n%03d d%d overall simulation      ", node, device);
  host_report_speed (t2[0], t2[2], ITER_WARMUP + ITER_SWAP, event);
  putchar ('\n');


  cudaFree (gpuseed_dev);
  cudaFree (gpuseed1_dev);
  cudaFree (cnt_dev);
  cudaFree (lattice_dev);
  cudaFree (lattice1_dev);
  cudaFree (temp_dev);
  cudaFree (st_dev);

  free (gpuseed1);
  free (cnt);
  free (lattice);
  free (temp);
  free (st);
}
