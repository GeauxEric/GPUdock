/*
  PATTERN = CHECKERBOARD
  DENSE   = COMPACT
  MSCT    = 4
*/



__device__ void
stencil (int iter)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB
  curandState myseed = pgpuseed[TperB * blockIdx.x + bidx];	// needed by curand

  /// temperature scratchpad
  __shared__ float temp_prob_shared[NBETA_PER_WORD][NPROB_MAX];
  gpu_init_temp (temp_prob_shared, bidx);


  /// lattice scratchpad
  // sizeof(u_int32_t) * 18 * 16 * 16 = 18 KB
  __shared__ DATATYPE l[SZ + 2][SZ][SZ];

  // index for import/export scratchpad
  const int yy = (SZ_HF * (threadIdx.z & 1)) + (threadIdx.y >> 1);
  const int xx = (SZ_HF * (threadIdx.y & 1)) + threadIdx.x;

  // index for reading scratch pad
  const int y = threadIdx.y;
  const int ya = (y + SZ - 1) % SZ;
  const int yb = (y + 1) % SZ;


  for (int word = 0; word < NWORD; word++) {

    // initilize temperature scatchpad
    gpu_compute_temp (temp_prob_shared, bidx, word);
    int lattice_offset = SZ_CUBE * NWORD * blockIdx.x + SZ_CUBE * word + bidx;


    for (int i = 0; i < iter; i++) {

      // two phases update
      for (int run = 0; run < 2; run++) {
	int x0 = (threadIdx.z & 1) ^ (threadIdx.y & 1) ^ run;	// initial x
	int x = (threadIdx.x << 1) + x0;
	//int xa = (x + SZ - 1) % SZ;
	//int xb = (x + 1) % SZ;
	int xa = (x + SZ - 1) & (SZ - 1);
	int xb = (x + 1) & (SZ - 1);


	for (int tile_offset = 0; tile_offset < SZ_CUBE;
	     tile_offset += SZ_TILE) {


	  // lattice scratchpad import
	  for (int z_offset = 0; z_offset < SZ; z_offset += (BDz0 >> 1)) {
	    int zz = z_offset + (threadIdx.z >> 1);
	    l[zz + 1][yy][xx] =
	      plattice[lattice_offset + tile_offset + (SZ * SZ * zz)];
	  }

	  // left and right bound
	  if (bidx < SZ * SZ) {
	    int tile_offset_a = (tile_offset + SZ_CUBE - SZ_TILE) % SZ_CUBE;
	    int tile_offset_b = (tile_offset + SZ_TILE) % SZ_CUBE;
	    l[0][yy][xx] =
	      plattice[lattice_offset + tile_offset_a + (SZ * SZ * (SZ - 1))];
	    l[SZ + 1][yy][xx] = plattice[lattice_offset + tile_offset_b];
	  }
	  __syncthreads ();


	  for (int z_offset = 0; z_offset < SZ; z_offset += BDz0) {
	    int z = z_offset + threadIdx.z + 1;

	    DATATYPE c = l[z][y][x];	// center
	    DATATYPE n0 = l[z][y][xa];	// left
	    DATATYPE n1 = l[z][y][xb];	// right
	    DATATYPE n2 = l[z][ya][x];	// up
	    DATATYPE n3 = l[z][yb][x];	// down
	    DATATYPE n4 = l[z - 1][y][x];	// front
	    DATATYPE n5 = l[z + 1][y][x];	// back

	    n0 = MASK_A * ((c >> SHIFT_J0) & 1) ^ n0 ^ c;
	    n1 = MASK_A * ((c >> SHIFT_J1) & 1) ^ n1 ^ c;
	    n2 = MASK_A * ((c >> SHIFT_J2) & 1) ^ n2 ^ c;
	    n3 = MASK_A * ((c >> SHIFT_J3) & 1) ^ n3 ^ c;
	    n4 = MASK_A * ((c >> SHIFT_J4) & 1) ^ n4 ^ c;
	    n5 = MASK_A * ((c >> SHIFT_J5) & 1) ^ n5 ^ c;

	    for (int s = 0; s < MSCT; s++) {
	      DATATYPE h =
		((n0 >> s) & MASK_S) +
		((n1 >> s) & MASK_S) +
		((n2 >> s) & MASK_S) +
		((n3 >> s) & MASK_S) +
		((n4 >> s) & MASK_S) + ((n5 >> s) & MASK_S);
	      h = (h << 1) + ((c >> s) & MASK_S);

	      DATATYPE flip = 0;

	      for (int shift = 0; shift < 25; shift += MSCT) {
		float val = temp_prob_shared[shift >> 2][(h >> shift) & MASK_E];
		float myrand = curand_uniform (&myseed);	// range: [0,1]
		flip |= ((myrand < val) << shift);	// myrand < val ? 1 : 0;
	      }


	      flip = flip << s;
	      l[z][y][x] = c ^ flip;
	    }

	  }			// end of "for z_offset" (compute inside a tile)

	  // lattice scratchpad export
	  for (int z_offset = 0; z_offset < SZ; z_offset += (BDz0 >> 1)) {
	    int zz = z_offset + (threadIdx.z >> 1);
	    plattice[lattice_offset + tile_offset + (SZ * SZ * zz)] =
	      l[zz + 1][y][x];
	  }
	  __syncthreads ();

	}			// end of "for tile_offset"

      }				// end of "for run" (1 iteration of QM)
      //__syncthreads ();
    }				// end of "for i" (multiple iteration of QM)
  }				// end of "for word"
}
