/*
  PATTERN = TIMEITER
  DENSE   = COMPACT
  MSCT    = 4
*/


#if MSCT != 1

__device__ void
stencil_11 (int iter)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB
  curandState myseed = pgpuseed[TperB * blockIdx.x + bidx];	// needed by curand

  /// temperature scratchpad
  __shared__ float temp_prob_shared[NBETA_PER_WORD][NPROB_MAX];
  gpu_init_temp (temp_prob_shared, bidx);


  /// lattice scratchpad
  // sizeof(u_int32_t) * 16 * 16 * 16 = 16 KB
  __shared__ DATATYPE l[SZ][SZ][SZ];


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


    // lattice scratchpad import
    for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
      int z = z_offset + threadIdx.z;
      l[z][y][x] = plattice[lattice_offset + (SZ * SZ * z)];
    }
    __syncthreads ();


    //int x0, x, xa, xb, z, za, zb;
    //DATATYPE center, e0, e1, e2, e3, e4, e5;

    for (int i = 0; i < iter; i++) {


      // data reuse among zi ???
      for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
	int z = z_offset + threadIdx.z;
	//int za = (z + SZ - 1) % SZ;
	//int zb = (z + 1) % SZ;
	int za = (z + SZ - 1) & (SZ - 1);
	int zb = (z + 1) & (SZ - 1);

	DATATYPE c = l[z][y][x];	// center
	DATATYPE n0 = l[z][y][xa];	// left
	DATATYPE n1 = l[z][y][xb];	// right
	DATATYPE n2 = l[z][ya][x];	// up
	DATATYPE n3 = l[z][yb][x];	// down
	DATATYPE n4 = l[za][y][x];	// front
	DATATYPE n5 = l[zb][y][x];	// back

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



	// process a lattice element
	for (int j = 0; j < 3; j += 2) {
	  int s = (iter & 1) + j;
	  //int s = 0;
	  DATATYPE h =
	    ((n0 >> s) & MASK_S) +
	    ((n1 >> s) & MASK_S) +
	    ((n2 >> s) & MASK_S) +
	    ((n3 >> s) & MASK_S) +
	    ((n4 >> s) & MASK_S) + ((n5 >> s) & MASK_S);
	  h = (h << 1) + ((c >> s) & MASK_S);

	  DATATYPE flip = 0;

	  //#pragma unroll
	  for (int shift = 0; shift < 25; shift += 4) {
	    float val = temp_prob_shared[shift >> 2][(h >> shift) & MASK_E];
	    float myrand = curand_uniform (&myseed);	// range: [0,1]
	    flip |= ((myrand < val) << shift);	// myrand < val ? 1 : 0;
	  }

	  flip = flip << s;
	  l[z][y][x] &= flip;
	}

	//*/

      }				// end of "for z_offset" (compute inside a tile)

      __syncthreads ();
    }				// end of "for i" (multiple iteration of QM)

    // lattice scratchpad export
    for (int z_offset = 0; z_offset < SZ; z_offset += BDz1) {
      int z = z_offset + threadIdx.z;
      plattice[lattice_offset + (SZ * SZ * z)] = l[z][y][x];
    }
    __syncthreads ();

  }				// end of "for word"
}


#endif
