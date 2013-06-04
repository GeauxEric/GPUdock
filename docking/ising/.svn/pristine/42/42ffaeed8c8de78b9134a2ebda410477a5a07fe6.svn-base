/*
  PATTERN = TIMEITER
  DENSE   = SPARSE
  MSCT    = 3 4
*/


#if MSCT != 1


__device__ void
stencil_10 (int iter)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB
  curandState myseed = pgpuseed[TperB * blockIdx.x + bidx];	// needed by curand

  /// temperature scratchpad
  __shared__ float temp_prob_shared[NBETA_PER_WORD][NPROB_MAX];
  gpu_init_temp (temp_prob_shared, bidx);


  /// lattice scratchpad
  // sizeof(u_int32_t) * 18 * 16 * 16 = 18 KB
  __shared__ DATATYPE l[SZ + 2][SZ][SZ];


  // index for reading scratch pad
  const int y = threadIdx.y;
  const int ya = (y + SZ - 1) % SZ;
  const int yb = (y + 1) % SZ;
  const int x = threadIdx.x;
  const int xa = (x + SZ - 1) % SZ;
  const int xb = (x + 1) % SZ;



  for (int word = 0; word < NWORD; word++) {

    // initilize temperature scatchpad
    gpu_compute_temp (temp_prob_shared, bidx, word);
    int lattice_offset = SZ_CUBE * NWORD * blockIdx.x + SZ_CUBE * word + bidx;


    for (int i = 0; i < iter; i++) {
      int oddeven = iter & 1;


      for (int tile_offset = 0; tile_offset < SZ_CUBE; tile_offset += SZ_TILE) {

	// lattice scratchpad import
	for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
	  int z = z_offset + threadIdx.z;
	  l[z + 1][y][x] =
	    plattice[lattice_offset + tile_offset + (SZ * SZ * z)];
	}

	// left and right bound
	if (bidx < SZ * SZ) {
	  int tile_offset_a = (tile_offset + SZ_CUBE - SZ_TILE) % SZ_CUBE;
	  int tile_offset_b = (tile_offset + SZ_TILE) % SZ_CUBE;
	  l[0][y][x] =
	    plattice[lattice_offset + tile_offset_a + (SZ * SZ * (SZ - 1))];
	  l[SZ + 1][y][x] = plattice[lattice_offset + tile_offset_b];
	}
	__syncthreads ();



	// data reuse among zi ???
	for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
	  int z = z_offset + threadIdx.z + 1;

	  DATATYPE c = l[z][y][x];	// center
	  DATATYPE n0 = l[z][y][xa];	// left
	  DATATYPE n1 = l[z][y][xb];	// right
	  DATATYPE n2 = l[z][ya][x];	// up
	  DATATYPE n3 = l[z][yb][x];	// down
	  DATATYPE n4 = l[z - 1][y][x];	// front
	  DATATYPE n5 = l[z + 1][y][x];	// back

	  // for profiling purpose
	  //float energy = (e0 ^ e1 ^ e2 ^ e3 ^ e4 ^ e5) & 2;
	  //float val = 0.7;
	  //float myrand = 0.4;
	  //DATATYPE flip = n0 ^ n1 ^ n2 ^ n3 ^ n4 ^ n5;
	  //flip |= ((curand_uniform (&myseed) < val) << b);

	  ///*
	  n0 = MASK_A * ((c >> SHIFT_J0) & 1) ^ n0 ^ c;
	  n1 = MASK_A * ((c >> SHIFT_J1) & 1) ^ n1 ^ c;
	  n2 = MASK_A * ((c >> SHIFT_J2) & 1) ^ n2 ^ c;
	  n3 = MASK_A * ((c >> SHIFT_J3) & 1) ^ n3 ^ c;
	  n4 = MASK_A * ((c >> SHIFT_J4) & 1) ^ n4 ^ c;
	  n5 = MASK_A * ((c >> SHIFT_J5) & 1) ^ n5 ^ c;
	  c = (c >> oddeven) & MASK_S;
	  n0 = (n0 >> oddeven) & MASK_S;
	  n1 = (n1 >> oddeven) & MASK_S;
	  n2 = (n2 >> oddeven) & MASK_S;
	  n3 = (n3 >> oddeven) & MASK_S;
	  n4 = (n4 >> oddeven) & MASK_S;
	  n5 = (n5 >> oddeven) & MASK_S;


#if MSCT == 3
	  DATATYPE h = n0 + n1 + n2 + n3 + n4 + n5;
#endif
#if MSCT == 4
	  DATATYPE h = n0 + n1 + n2 + n3 + n4 + n5;
	  h = (h << 1) + c;
#endif

	  // process a lattice element
	  DATATYPE flip = 0;




#if MSCT == 3
	  //#pragma unroll
	  for (int shift = 0; shift < 29; shift += MSCT) {
	    int energy = h >> shift & MASK_E;	// range: [0,6]
	    int spin = c >> shift & 1;
	    float val = temp_prob_shared[shift / MSCT][(energy << 1) + spin];
	    float myrand = curand_uniform (&myseed);	// range: [0,1]
	    flip |= ((myrand < val) << shift);	// myrand < val ? 1 : 0;
	  }
#endif



#if MSCT == 4
	  //#pragma unroll
	  for (int shift = 0; shift < LENG; shift += MSCT) {
	    float val = temp_prob_shared[shift >> 2][(h >> shift) & MASK_E];
	    float myrand = curand_uniform (&myseed);	// range: [0,1]
	    flip |= ((myrand < val) << shift);	// myrand < val ? 1 : 0;
	  }
#endif


	  flip = flip << (!oddeven);
	  //*/
	  l[z][y][x] &= flip;

	}			// end of "for z_offset" (compute inside a tile)


	// lattice scratchpad export
	for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
	  int z = z_offset + threadIdx.z;
	  plattice[lattice_offset + tile_offset + (SZ * SZ * z)] =
	    l[z + 1][y][x];
	}

	__syncthreads ();

      }				// end of "for tile_offset"

    }				// end of "for i" (multiple iteration of QM)


  }				// end of "for word"
}


#endif
