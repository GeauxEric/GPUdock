#include "COPYING"


/*
  ALLOCATION = SHARED
  DENSE      = COMPACT
  MSCT       = 4
*/

__device__ void
stencil_unified (int *temp_idx_shared, float *temp_beta_shared, float *E, int iter, int swap_mod, int mod)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB


  curandState seed0 = pgpuseed[TperB * blockIdx.x + bidx];

  // signed 16 bit integer: -32K ~ 32K, never overflows
  // sizeof (shot) * 24 * 512 = 24 KB
  __shared__ short E_shared[NBETA_PER_WORD][TperB];
  __shared__ float __align__ (32) Eh[NBETA];


  /// temperature scratchpad
  __shared__ PROB_DATATYPE temp_prob_shared[NBETA_PER_WORD][NPROB_MAX];
  gpu_init_temp (temp_prob_shared, bidx);

  /// lattice scratchpad
  // sizeof(u_int32_t) * 16 * 16 * 16 = 16 KB
  __shared__ MSC_DATATYPE l[SZ][SZ][SZ];

  // index for read/write glocal memory
  const int xx = (SZ_HF * (threadIdx.y & 1)) + threadIdx.x;
  const int yy = (SZ_HF * (threadIdx.z & 1)) + (threadIdx.y >> 1);

  // index for reading scratchpad
  const int y = threadIdx.y;
  //const int ya = (y + SZ - 1) % SZ;
  //const int yb = (y + 1) % SZ;
  const int ya = (y + SZ - 1) & (SZ - 1);
  const int yb = (y + 1) % (SZ - 1);


  for (int word = 0; word < NWORD; word++) {
    int lattice_offset = SZ_CUBE * NWORD * blockIdx.x + SZ_CUBE * word;

    // initilize temperature scratchpad
    gpu_compute_temp (temp_prob_shared, temp_beta_shared, bidx, word);

    // import lattice scratchpad
    for (int z_offset = 0; z_offset < SZ; z_offset += (BDz0 >> 1)) {
      int zz = z_offset + (threadIdx.z >> 1);
      l[zz][yy][xx] = plattice[lattice_offset + SZ * SZ * z_offset + bidx];
    }

    // reset partial status
    for (int b = 0; b < NBETA_PER_WORD; b++)
      E_shared[b][bidx] = 0;

    __syncthreads ();



    for (int i = 0; i < iter; i++) {
      // two phases update
      for (int run = 0; run < 2; run++) {
	int mybool = (mod == 1) && (i % ITER_SWAP_KERNFUNC == 0) && (run == 0);

	int x0 = (threadIdx.z & 1) ^ (threadIdx.y & 1) ^ run;	// initial x
	int x = (threadIdx.x << 1) + x0;
	int x1 = x + 1 - x0 * 2;
	//int xa = (x + SZ - 1) % SZ;
	//int xb = (x + 1) % SZ;
	int xa = (x + SZ - 1) & (SZ - 1);
	int xb = (x + 1) & (SZ - 1);

	// data reuse among z ???
	for (int z_offset = 0; z_offset < SZ; z_offset += BDz0) {
	  int z = z_offset + threadIdx.z;
	  //int za = (z + SZ - 1) % SZ;
	  //int zb = (z + 1) % SZ;
	  int za = (z + SZ - 1) & (SZ - 1);
	  int zb = (z + 1) & (SZ - 1);

	  MSC_DATATYPE c = l[z][y][x];	// center
	  MSC_DATATYPE n0 = l[z][y][xa];	// left
	  MSC_DATATYPE n1 = l[z][y][xb];	// right
	  MSC_DATATYPE n2 = l[z][ya][x];	// up
	  MSC_DATATYPE n3 = l[z][yb][x];	// down
	  MSC_DATATYPE n4 = l[za][y][x];	// front
	  MSC_DATATYPE n5 = l[zb][y][x];	// back

	  n0 = MASK_A * ((c >> SHIFT_J0) & 1) ^ n0 ^ c;
	  n1 = MASK_A * ((c >> SHIFT_J1) & 1) ^ n1 ^ c;
	  n2 = MASK_A * ((c >> SHIFT_J2) & 1) ^ n2 ^ c;
	  n3 = MASK_A * ((c >> SHIFT_J3) & 1) ^ n3 ^ c;
	  n4 = MASK_A * ((c >> SHIFT_J4) & 1) ^ n4 ^ c;
	  n5 = MASK_A * ((c >> SHIFT_J5) & 1) ^ n5 ^ c;

	  MSC_DATATYPE flip = 0;
	  for (int offset = 0; offset < NBETA_PER_SEG; ++offset) {
	    MSC_DATATYPE e =
	      ((n0 >> offset) & MASK_S) +
	      ((n1 >> offset) & MASK_S) +
	      ((n2 >> offset) & MASK_S) +
	      ((n3 >> offset) & MASK_S) +
	      ((n4 >> offset) & MASK_S) +
	      ((n5 >> offset) & MASK_S);
	    e = (e << 1) + ((c >> offset) & MASK_S);

	    for (int seg = 0; seg < NSEG_PER_WORD; ++seg) {
	      int b = NBETA_PER_SEG * seg + offset;
	      int position = (seg << 2) + offset;
	      PROB_DATATYPE val = temp_prob_shared[b][(e >> (seg << 2)) & MASK_E];
	      PROB_DATATYPE myrand = curand (&seed0);	// range: [0,U_INT32_T_MAX]
	      flip |= ((myrand < val) << position);	// myrand < val ? 1 : 0;

	    if (mybool)
	      E_shared[b][bidx] += e ^ (myrand < val);
	    }
	  }

	  l[z][y][x] = c ^ flip;

	}			// z_offset


	if (mybool) {
	  gpu_reduction (E, E_shared, bidx, word);
	  __syncthreads ();

	  /// energy contribute by external field
	  for (int b = 0; b < NBETA_PER_WORD; b++)
	    E_shared[b][bidx] = 0;
    
	  for (int z_offset = 0; z_offset < SZ; z_offset += BDz0) {
	    int z = z_offset + threadIdx.z;
	    MSC_DATATYPE c0 = l[z][y][x];
	    MSC_DATATYPE c1 = l[z][y][x1];

	    for (int s = 0; s < NBETA_PER_SEG; s++) {
	      for (int shift = 0; shift < SHIFT_MAX; shift += MSCT) {
		int ss = shift + s;
		E_shared[ss][bidx] += ((c0 >> ss) & 1) + ((c1 >> ss) & 1);
	      }
	    }
	  }

	  gpu_reduction (Eh, E_shared, bidx, word);
	  __syncthreads ();


	  if (bidx < NBETA) {
	    E[bidx] = E[bidx] * 2 - 6 * (SZ / BDz0) * TperB;
	    Eh[bidx] = Eh[bidx] * 2 - SZ_CUBE;
	    E[bidx] = E[bidx] + Eh[bidx] * H;
	  }
	  __syncthreads ();


	  gpu_shuffle (temp_idx_shared, temp_beta_shared, E, bidx, mod);
	}


	__syncthreads ();
      }				// run



    

    }				// i

    // export lattice scratchpad
    for (int z_offset = 0; z_offset < SZ; z_offset += (BDz0 >> 1)) {
      int zz = z_offset + (threadIdx.z >> 1);
      plattice[lattice_offset + SZ * SZ * z_offset + bidx] = l[zz][yy][xx];
    }

  }				// word


  // copy seed back
  pgpuseed[TperB * blockIdx.x + bidx] = seed0;
}







__global__ void
kernel_unified (int rec, int mod)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB

  /// temperature
  // (4 * 32) * 2 = 256 B
  __shared__ int __align__ (32) temp_idx_shared[NBETA];
  __shared__ float __align__ (32) temp_beta_shared[NBETA];

  /// lattice energy
  // sizeof (float) * 32 = 128 B
  __shared__ float __align__ (32) E[NBETA];

  // load temperature
  if (bidx < NBETA) {
    temp_idx_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].idx;
    temp_beta_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].beta;
  }


  if (mod == 0)
    for (int i = 0; i < ITER_WARMUP_KERN; i += ITER_WARMUP_KERNFUNC) {
      stencil_unified (temp_idx_shared, temp_beta_shared, E, ITER_WARMUP_KERNFUNC, 0, 0);
    }
  else
    for (int i = 0; i < ITER_SWAP_KERN; i += ITER_SWAP_KERNFUNC) {
      int swap_mod = (i / ITER_SWAP_KERNFUNC) & 1;
      stencil_unified (temp_idx_shared, temp_beta_shared, E, ITER_SWAP_KERNFUNC, swap_mod, 1);
    }



  // store temperature
  if (bidx < NBETA) {
    ptemp[NBETA_MAX * blockIdx.x + bidx].idx = temp_idx_shared[bidx];
    ptemp[NBETA_MAX * blockIdx.x + bidx].beta = temp_beta_shared[bidx];
  }
  __syncthreads ();

  // store energy status
  //  if (bidx < NBETA)
  //  pst[rec].e[blockIdx.x][temp_idx_shared[bidx]] = E[bidx];
  //__syncthreads ();
}





