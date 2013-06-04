#include "COPYING"




// initilize seeds for curand
// CURAND_Library.pdf, pp21
__global__ void
kernel_init_seed ()
{
  const int gidx = TperB * blockIdx.x + threadIdx.x;
  curand_init (seed, gidx, 0, &pgpuseed[gidx]);

  // seed, subsequence, offset, gpuseed
  // skipahead(100000, &gpuseed[gidx]);
}





__global__ void
kernel_warmup ()
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB

  /// temperature
  // (4 * 32) * 2 = 256 B
  __shared__ float __align__ (32) temp_beta_shared[NBETA];
  if (bidx < NBETA)
    temp_beta_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].beta;


  for (int i = 0; i < ITER_WARMUP_KERN; i += ITER_WARMUP_KERNFUNC) {
    stencil (temp_beta_shared, ITER_WARMUP_KERNFUNC);
  }
}




__global__ void
kernel_swap (int rec)
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


  for (int i = 0; i < ITER_SWAP_KERN; i += ITER_SWAP_KERNFUNC) {
    int swap_mod = (i / ITER_SWAP_KERNFUNC) & 1;
    stencil_swap (temp_idx_shared, temp_beta_shared, E, swap_mod);
    stencil (temp_beta_shared, ITER_SWAP_KERNFUNC);
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






__global__ void
kernel_compute_q (int rec)
{
  const int bidx = threadIdx.x;

  // sizeof(u_int32_t) * 16 * 16 * 16 = 16 KB 
  __shared__ MSC_DATATYPE l1[SZ_CUBE];

  const int lattice_offset0 = SZ_CUBE * (blockIdx.x << 1);
  const int lattice_offset1 = lattice_offset0 + SZ_CUBE;

  for (int offset = 0; offset < SZ_CUBE; offset += TperB) {
    l1[offset + bidx] =		// xord_word 
      plattice1[lattice_offset0 + offset + bidx] ^
      plattice1[lattice_offset1 + offset + bidx];
  }

  // is double an overkill?
  if (bidx < NBETA) {
    float q0 = 0.0f;
    double qk_real = 0.0;
    double qk_imag = 0.0;
    double qk2_real = 0.0;
    double qk2_imag = 0.0;
    double angel1, angel2;

    MSC_DATATYPE xor_word;
    int xor_bit;

    for (int i = 0; i < SZ_CUBE; i++) {
      xor_word = l1[i];
      xor_bit = (xor_word >> bidx) & 0x1;
      xor_bit = 1 - (xor_bit << 1);	// parallel: +1, reverse: -1

      // 2 * pi / L * x_i
      angel1 = (double) (i % SZ) * 2 * PI / SZ;
      // 2 * pi / L * (x_i + y_i)
      angel2 = (double) (i % SZ + (i / SZ) % SZ) * 2 * PI / SZ;

      q0 += xor_bit;
      qk_real += xor_bit * cos (angel1);
      qk_imag += xor_bit * sin (angel1);
      qk2_real += xor_bit * cos (angel2);
      qk2_imag += xor_bit * sin (angel2);
    }

    // save measurements in "st"
    pst[rec].q[blockIdx.x][bidx] = q0;
    pst[rec].qk_real[blockIdx.x][bidx] = qk_real;
    pst[rec].qk_imag[blockIdx.x][bidx] = qk_imag;
    //      (float) sqrt (qk_real * qk_real + qk_imag * qk_imag);
    pst[rec].qk2_real[blockIdx.x][bidx] = qk2_real;
    pst[rec].qk2_imag[blockIdx.x][bidx] = qk2_imag;
    //      (Float) sqrt (qk2_real * qk2_real + qk2_imag * qk2_imag);
  }
}


//not used
#if 0
__global__ void
kernel_compute_q (int rec)
{
  const int bidx = threadIdx.x;

  // sizeof(u_int32_t) * 16 * 16 * 16 = 16 KB
  __shared__ MSC_DATATYPE l0[SZ_CUBE];
  __shared__ MSC_DATATYPE l1[SZ_CUBE];

  //  __shared__ float m0[NBETA];
  //__shared__ float m1[NBETA];

  const int lattice_offset0 = SZ_CUBE * (blockIdx.x << 1);
  const int lattice_offset1 = lattice_offset0 + SZ_CUBE;


  for (int offset = 0; offset < SZ_CUBE; offset += TperB) {
    l0[offset + bidx] =		// xord_word
      plattice1[lattice_offset0 + offset + bidx];
    l1[offset + bidx] = plattice1[lattice_offset1 + offset + bidx];
  }

  __syncthreads();

  if (bidx < NBETA) {
    float m0 = 0.0f;
    float m1 = 0.0f;
    MSC_DATATYPE xor_word1, xor_word0;
    int xor_bit1, xor_bit0;
    for (int i = 0; i < SZ_CUBE; i++) {
      xor_word1 = l1[i];
      xor_word0 = l0[i];
      xor_bit1 = (xor_word1 >> bidx) & 0x1;
      xor_bit0 = (xor_word0 >> bidx) & 0x1;
      xor_bit1 = (xor_bit1 << 1) - 1;	// parallel: +1, reverse: -1
      xor_bit0 = (xor_bit0 << 1) - 1;	// parallel: +1, reverse: -1

      m1 += xor_bit1;
      m0 += xor_bit0;
    }
    m0 /= (float)SZ_CUBE;
    m1 /= (float)SZ_CUBE;


    float q0 = 0.0f;
    float qk_real = 0.0f;
    float qk_imag = 0.0f;
    float qk2_real = 0.0f;
    float qk2_imag = 0.0f;

    //    MSC_DATATYPE xor_word0,xor_word1;
    //int xor_bit0,xor_bit1;
    for (int i = 0; i < SZ_CUBE; i++) {
      xor_word1 = l1[i];
      xor_word0 = l0[i];
      xor_bit0 = (xor_word0 >> bidx) & 0x1;
      xor_bit1 = (xor_word1 >> bidx) & 0x1;
      //      xor_bit = 1 - (xor_bit << 1);     // parallel: +1, reverse: -1
      xor_bit0 = (xor_bit0 << 1) - 1;
      xor_bit1 = (xor_bit1 << 1) - 1;

      float xor_bit = ((float)xor_bit0 - m0) * ((float)xor_bit1 - m1);

      q0 += xor_bit;
      qk_real += xor_bit * cos ((double) (i % SZ) * 2 * PI / SZ);	// cos(2pi/L*x_i);
      qk_imag += xor_bit * sin ((double) (i % SZ) * 2 * PI / SZ);	// sin(2pi/L*x_i)
      qk2_real += xor_bit * cos ((double) (i % SZ + (i / SZ) % SZ) * 2 * PI / SZ);	// cos(2pi/L*(x_i+y_i))
      qk2_imag += xor_bit * sin ((double) (i % SZ + (i / SZ) % SZ) * 2 * PI / SZ);	// sin(2pi/L*(x_i+y_i))
    }

    // save measurements in "st"
    pst[rec].q[blockIdx.x][bidx] = q0;
    pst[rec].qk[blockIdx.x][bidx] =
      (float) sqrt (qk_real * qk_real + qk_imag * qk_imag);
    pst[rec].qk2[blockIdx.x][bidx] =
      (float) sqrt (qk2_real * qk2_real + qk2_imag * qk2_imag);
  }
  
  
}

#endif
