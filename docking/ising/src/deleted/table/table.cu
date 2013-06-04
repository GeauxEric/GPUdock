#include "table.cuh"
#define PRINT


__constant__ float exp_beta_e2_dev[NBETA_ALL];
__constant__ float exp_beta_e4_dev[NBETA_ALL];
__constant__ float exp_beta_e6_dev[NBETA_ALL];


__global__ void
kernel_init_state (curandState * state, int seed)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init (seed ^ tid, tid, 0, &state[tid]);

  // seed, subsequence, offset, state
  // skipahead(100000,&state[tid]);
}



__global__ void
kernel_generate (u_int64_t * table_dev, curandState * state, int size)
{

  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  int nthreads = blockDim.x * gridDim.x;	// 64 * 16
  int iter = (size + nthreads - 1) / nthreads;	// (200 + 1024 - 1) / 1024
  int iter_offset;
  float myrand;
  curandState *mystate = state + tid;
  u_int64_t item, item_seg;
  int nbeta_offset = NBETA * blockIdx.x;
  int k;

  for (int i = 0; i < iter; i++) {
    iter_offset = nthreads * i;

    if (tid < size - iter_offset) {
      item = 0;
      for (int j = 0; j < NBETA; j++) {
	k = nbeta_offset + j;

	myrand = curand_uniform (mystate);	// range: (0,1]
	item_seg = 1;
	item_seg += (myrand < exp_beta_e2_dev[k]);
	item_seg += (myrand < exp_beta_e4_dev[k]);
	item_seg += (myrand < exp_beta_e6_dev[k]);
	item |= (item_seg << (j << 2));
      }
      table_dev[iter_offset + tid] = item;
    }

  }
}





void
check (u_int64_t * table, int size)
{
  u_int64_t item, item_seg;

  int cnt[16];
  for (int i = 0; i < 16; i++)
    cnt[i] = 0;

  for (int i = 0; i < size; i++) {
    item = table[i];
    for (int j = 0; j < NBETA; j++) {
      item_seg = (item >> (j << 2)) & 0xF;
      cnt[item_seg]++;
    }
#ifdef PRINT
//    printf ("%08d: \t %016llx\n", i, item);
#endif
  }

  putchar ('\n');
  printf ("the distribution among %d item segments:\n", size * NBETA);
  for (int i = 0; i < 16; i++) {
    printf ("\t#%02d : %10d\n", i, cnt[i]);
  }
  putchar ('\n');
}





void
init_beta ()
{
  const float beta_begin = BETA_BEGIN;
  const float beta_end = BETA_END;
  const float beta_delta = (beta_end - beta_begin) / NBETA_ALL;

  float *exp_beta_e2 = (float *) malloc (sizeof (float) * NBETA_ALL);
  float *exp_beta_e4 = (float *) malloc (sizeof (float) * NBETA_ALL);
  float *exp_beta_e6 = (float *) malloc (sizeof (float) * NBETA_ALL);

  float betai = beta_begin;
  for (int i = 0; i < NBETA_ALL; i++) {
    exp_beta_e2[i] = exp (-2 * betai);
    exp_beta_e4[i] = exp (-4 * betai);
    exp_beta_e6[i] = exp (-6 * betai);
    betai += beta_delta;
  }

  cudaMemcpyToSymbol (exp_beta_e2_dev, exp_beta_e2,
		      sizeof (float) * NBETA_ALL);
  cudaMemcpyToSymbol (exp_beta_e4_dev, exp_beta_e4,
		      sizeof (float) * NBETA_ALL);
  cudaMemcpyToSymbol (exp_beta_e6_dev, exp_beta_e6,
		      sizeof (float) * NBETA_ALL);

  free (exp_beta_e2);
  free (exp_beta_e4);
  free (exp_beta_e6);
}




void
run (int size)
{

  u_int64_t *table, *table_dev;
  table = (u_int64_t *) malloc (sizeof (u_int64_t) * size);
  cudaMalloc ((void **) &table_dev, sizeof (u_int64_t) * size);

  curandState *state_dev;
  cudaMalloc ((void **) &state_dev, sizeof (curandState) * NT);

  srand (time (NULL));
  int seed = rand ();



  init_beta ();
  kernel_init_state <<< GD, BD >>> (state_dev, seed);

  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);
  cudaEventRecord (start);

  kernel_generate <<< GD, BD >>> (table_dev, state_dev, size);
  cudaMemcpy (table, table_dev, sizeof (u_int64_t) * size,
	      cudaMemcpyDeviceToHost);

  cudaEventRecord (stop);
  cudaEventSynchronize (stop);
  float elapsed_time;
  cudaEventElapsedTime (&elapsed_time, start, stop);

  check (table, size);
  printf ("CUDA Time for generating %d items: \t%8.6f ms\n", size, time);

  cudaEventDestroy (start);
  cudaEventDestroy (stop);
  cudaFree (table_dev);
  cudaFree (state_dev);
  free (table);
}



int
main (int argc, char **argv)
{
  int size = SZ;

  if (argc == 2) {
    size = atoi (argv[1]);
  }

  run (size);

  return 0;
}
