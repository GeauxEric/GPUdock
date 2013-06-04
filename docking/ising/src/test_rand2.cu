/*
  ALLOCATION = SHARED
  DENSE      = COMPACT
  MSCT       = 4
*/


__device__ void
stencil (float *temp_beta_shared, int iter)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB


  curandState seed0 = pgpuseed[TperB * blockIdx.x + bidx];
  LCG_DATATYPE seed1 = pgpuseed1[TperB * blockIdx.x + bidx];

  int rand123_cnt = pcnt[TperB * blockIdx.x + bidx];
  philox4x32_key_t rand123_k = {{bidx, seed1}};
  philox4x32_ctr_t rand123_c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};
  union {
    philox4x32_ctr_t c;
    uint4 i;
  } rand123_u;


  PROB_DATATYPE aaa = (PROB_DATATYPE) plattice[TperB * blockIdx.x + bidx];

  for (int i = 0; i < iter; i++) {
    for (int b = 0; b < NBETA; b++) {
      PROB_DATATYPE myrand = curand_uniform (&seed0);
      //PROB_DATATYPE myrand = curand (&seed0);
      gpu_lcg (&seed1);
      aaa += myrand;
    }
  }

  plattice[TperB * blockIdx.x + bidx] = (MSC_DATATYPE) aaa;


  // copy seed back
  pgpuseed[TperB * blockIdx.x + bidx] = seed0;
  pgpuseed1[TperB * blockIdx.x + bidx] = seed1;
  pcnt[TperB * blockIdx.x + bidx] = rand123_cnt;
}


__device__ void
stencil_swap (int *temp_idx_shared, float *temp_beta_shared, float *E, int mod)
{
}


__global__ void
kernel_rearrange ()
{
}

