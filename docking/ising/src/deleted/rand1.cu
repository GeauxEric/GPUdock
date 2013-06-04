#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#define GD 16
#define BD 64
#define NT (BD * GD)
#define SZ 2000

# define CUDA_CALL ( x ) do { if (( x ) != cudaSuccess ) { \
      printf (" Error at % s :% d \ n " , __FILE__ , __LINE__ ) ; \
      return EXIT_FAILURE ;}} while (0)

# define CURAND_CALL ( x ) do { if (( x ) != CURAND_STATUS_SUCCESS ) { \
      printf (" Error at % s :% d \ n " , __FILE__ , __LINE__ ) ; \
      return EXIT_FAILURE ;}} while (0)





__global__ void
kernel_init_state (curandState * state, int seed)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init (seed, tid, 0, &state[tid]);
  // seed, subsequence, offset, state
  // skipahead(100000,&state[tid]);
}


__global__ void
kernel_rand (float *table, curandState * state, int size)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int iter = (size + NT - 1) / NT;
  int iter_offset;
  float myrand;

  for (int i = 0; i < iter; i++) {
    iter_offset = NT * iter;
    if (tid < size - iter_offset) {
      myrand = curand_uniform (&state[tid]);	// range: [0,1]
      table[iter_offset + tid] = myrand;
    }
  }
}


void
check_table (float *table, int size)
{
  float item;
  for (int i = 0; i < size; i++) {
    item = table[i];
    printf ("%08d:\t%8.8f\n", i, item);
  }
}



void
run (int seed, int size)
{
  float *table, *table_dev;
  table = (float *) malloc (sizeof (float) * SZ);
  cudaMalloc ((void **) &table_dev, sizeof (float) * SZ);

  curandState *state_dev;
  cudaMalloc ((void **) &state_dev, sizeof (curandState) * NT);

  kernel_init_state <<< GD, BD >>> (state_dev, seed);
  kernel_rand <<< GD, BD >>> (table_dev, state_dev, size);

  cudaMemcpy (table, table_dev, sizeof (float) * SZ, cudaMemcpyDeviceToHost);

  check_table (table, size);

  free (table);
  cudaFree (table_dev);
  cudaFree (state_dev);
}


int
main (int argc, char **argv)
{
  srand (time(NULL));
  int seed = rand ();
  //int seed = 1234234;
  int size = SZ;

  if (argc == 2) {
    seed = atoi (argv[1]);
  }
  if (argc == 3) {
    seed = atoi (argv[1]);
    size = atoi (argv[2]);
  }

  run (seed, size);

  return 0;
}
