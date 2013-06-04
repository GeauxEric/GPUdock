


#ifdef DEBUG_M

  __shared__ short M_shared[NBETA_TperB];
  __shared__ int M[LENG]; // reduction of M_shared
  __shared__ int M2[LENG]; // reverse shuffle to original order

...

for (int b = 0; b < NBETA; b++)
  M_shared[TperB * b + bidx] = 0;

...

M_shared[TperB * b + bidx] += spin ^ (myrand < val);

...


gpu_reduction (E_shared, E, bidx);
#if TEMP_ALG == 0
      if (bidx < LENG)
	E[bidx] = E[bidx] / 2;
#elif TEMP_ALG == 1
      // convert from [0,6] to [-6,6], e = e * 2 - 6
      // E = sum_TperB sum_run sum_ZITER (e * 2 - 6)
      //   = sum_TperB sum_run (2 * sum_ZITER e - 6 * ZITER)
      //   = sum_TperB (2 * sum_ZITER_run e - 6 * ZITER * run)
      //   = (2 * sum_ZITER_run_TperB e - 6 * ZITER * run * TperB)
      if (bidx < LENG)
	E[bidx] = E[bidx] - 3 * ZITER * 2 * TperB;
#endif
      if (bidx < LENG)
	E2[temp_idx_shared[bidx]] = E[bidx];
      __syncthreads ();
      if (bidx < LENG)
	st[iter + iter0].E[LENG * blockIdx.x + bidx] = E2[bidx];


      gpu_reduction (M_shared, M, bidx);
#if TEMP_ALG == 1
      // convert from [0,1] to [-1,1], spin = spin * 2 - 1
      if (bidx < LENG)
	M[bidx] = M[bidx] * 2 - TperB * ZITER * 2;
#endif
      if (bidx < LENG)
	M2[temp_idx_shared[bidx]] = M[bidx];
      __syncthreads ();
      if (bidx < LENG)
	st[iter + iter0].M[LENG * blockIdx.x + bidx] = M2[bidx];

...












void cpu_conclude_st (St * st, Avrg * avrg)
{
  int E_rec = 0, M_rec = 0;
  int64_t E_sum = 0, M_sum = 0;
  double M1, M2_sum = 0, M4_sum = 0; // intermediates for calc U
  
  // M
  M_sum += (int64_t) abs (st[iter_rec].M[i]);
  M1 = (double) st[iter_rec].M[i];
  M2_sum += M1 * M1;
  M4_sum += M1 * M1 * M1 * M1;
  M_rec += 1;

  avrg->M[i] = int (M_sum / M_rec);
  avrg->U[i] = float (1 - M4_sum / (3 * M2_sum * M2_sum) * M_rec);
}
#endif





void
cpu_conclude_st (St * st, Avrg * avrg)
{
  for (int lattice = 0; lattice < GD; lattice++) {
    for (int b = 0; b < NBETA; b++) {
      int i = LENG * lattice + b;

      int E_rec = 0, M_rec = 0;
      int64_t E_sum = 0, M_sum = 0;
      double M1, M2_sum = 0, M4_sum = 0;	// intermediates for calc U

      for (int iter = ITER_WARMUP; iter < ITER_TOTAL; iter++) {
	int iter_rec = iter - ITER_WARMUP;
	// E
	E_sum += st[iter_rec].E[i];
	E_rec += 1;

	// M
	M_sum += (int64_t) abs (st[iter_rec].M[i]);
	M1 = (double) st[iter_rec].M[i];
	M2_sum += M1 * M1;
	M4_sum += M1 * M1 * M1 * M1;
	M_rec += 1;
      }

      avrg->E[i] = int (E_sum / E_rec);
      avrg->M[i] = int (M_sum / M_rec);
      avrg->U[i] = float (1 - M4_sum / (3 * M2_sum * M2_sum) * M_rec);
    }
  }
}













void
cpu_save_st_debug (Temp * temp, St * st, Avrg * avrg, char *mystime)
{
  char myfile[STR_LENG];
  sprintf (myfile, "output/st_%s_debug.txt", mystime);

  printf ("write debug output to\t %s\n", myfile);
  FILE *fp = fopen (myfile, "w");
  if (fp == NULL) {
    fprintf (stderr, "failed to open %s\n", myfile);
    exit (1);
  }

//  E = -6, -4, -2, 0, 2, 4, 6
//  M = -1, 1

  for (int iter = 0; iter < ITER_REC; iter++) {
    //for (int lattice = 0; lattice < GD; lattice++) {
    for (int lattice = 0; lattice < 1; lattice++) {
      fprintf (fp, "#lattice %2d\n", lattice);
      for (int b = 0; b < NBETA; b++) {
	int i = LENG * lattice + b;
	fprintf (fp, "%5u %6d %5d", iter, st[iter].E[i], st[iter].M[i]);
	fprintf (fp, "\t%2d %.6f\n", b, temp[temp[i].idx].beta);
      }
    }
  }
  fclose (fp);
}




void
cpu_save_st (Temp * temp, St * st, Avrg * avrg, char *mystime)
{
  char myfile[STR_LENG];
  sprintf (myfile, "output/st_%s.txt", mystime);

  printf ("write output to\t\t %s\n", myfile);
  FILE *fp = fopen (myfile, "w");
  if (fp == NULL) {
    fprintf (stderr, "failed to open %s\n", myfile);
    exit (1);
  }

  fprintf (fp, "{      E,      M,        U,         beta }\n\n");

  for (int lattice = 0; lattice < GD; lattice++) {
    fprintf (fp, "#lattice %2d\n", lattice);
    for (int b = 0; b < NBETA; b++) {
      int i = LENG * lattice + b;
      fprintf (fp, "{ %6d, %6d, %8.4f, %12.8f },\n",
	       avrg->E[temp[i].idx], avrg->M[temp[i].idx],
	       avrg->U[temp[i].idx], temp[i].beta);
    }
  }

  fclose (fp);
}







  /*
     // st - status records
     St *st, *st_dev;
     size_t st_sz = sizeof (St) * ITER_REC;
     //printf ("st_sz=%d MB\n", st_sz / 1024 / 1024);
     st = (St *) malloc (st_sz);
     cudaMalloc ((void **) &st_dev, st_sz);

cudaFree (st_dev);
free (st)
   */



/*
     cudaMemcpy (st, st_dev, st_sz, cudaMemcpyDeviceToHost);
     Avrg avrg;
     cpu_conclude_st (st, &avrg);
     cpu_save_st (temp, st, &avrg, mystime);
     cpu_save_st_debug (temp, st, &avrg, mystime);
*/











__device__ void
gpu_reduction0 (short *a_shared, int *a, const int bidx)
{
  short s = 0;
  if (bidx < NBETA) {
    int bb = TperB * bidx;
    for (int i = 0; i < 256; i++)
      s += a_shared[bb + i];	// bank conflict
    a[bidx] = s;
  }
}





// load beta from shared memory, compute, save prob in shared memory
__device__ void
gpu_compute_temp_s
(float *beta, float *prob, const int bidx)
{
  if (bidx < 8) {		// only compute first 8 elements of prob[14]
    for (int b = 0; b < NBETA; b++) {
      float mybeta = beta[b];
      float energy = bidx - 6 + H - (H * 2.0f + 1.0f) * (bidx % 2);
      prob[16 * b + bidx] = __expf (2 * energy * mybeta);
    }
  }
}







// did not varify the correctness
void
cpu_compute_q_1 (DATATYPE * lattice, Q * q, int rec, Temp * temp)
{
  DATATYPE xor_word;
  DATATYPE w0, w1; // lattice element, MSCT
  DATATYPE s0, s1; // rearranged spins, in compact format
  int temp_idx0[NBETA];
  int temp_idx1[NBETA];
  int q0[NBETA];


  for (int pair = 0; pair < GD; pair += 2) {	// for every pair

    for (int b = 0; b < NBETA_PER_WORD; b++) {
      q0[b] = 0;
      temp_idx0[b] = temp[NBETA_MAX * pair + b].idx;
      temp_idx1[b] = temp[NBETA_MAX * (pair + 1) + b].idx;
    }

    for (int i = 0; i < SZ_CUBE; i++) {	// for every element
      w0 = lattice[SZ_CUBE * pair + i];
      w1 = lattice[SZ_CUBE * (pair + 1) + i];
      s0 = 0;
      s1 = 0;

      // rearrange the spins so that they matches the order of temperature
      for (int b = 0; b < NBETA; b++) {
	int shift = MSCT * b;
	s0 |= ((w0 >> shift) & 1) << temp_idx0[b];
	s1 |= ((w1 >> shift) & 1) << temp_idx1[b];
      }

      xor_word = s0 ^ s1; // parallel: 0, reverse: 1
      for (int b = 0; b < NBETA; b++)
	q0[b] += 1 - 2 * ((xor_word >> b) & 1); // parallel: +1, reverse: -1
    }

    for (int b = 0; b < NBETA_PER_WORD; b++)
      q[rec].q0[pair / 2][b] = q0[b];

  }
}











//__global__ void
__device__ void
gpu_ising_swap (int mod)
{
  const int bidx = blockDim.x * blockDim.y * threadIdx.z + blockDim.x * threadIdx.y + threadIdx.x;	// within a TB

  /*
     if (bidx == 0 && blockIdx.x == 0)
     printf ("hello from __glocal__ gpu_ising_swap ()\n");
   */


  /// temperature scratchpad
  // sizeof(float) * 32 = 128 B
  __shared__ int __align__(32) temp_idx_shared[NBETA];
  // sizeof(float) * 32 = 128 B
  __shared__ float __align__(32) temp_beta_shared[NBETA];
  if (bidx < NBETA) {
    temp_idx_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].idx;
    temp_beta_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].beta;
  }



  /// E scratchpads
  /*
     Should "short" datatype degrade performance?
     involve less threads
   */

  // signed 16 bit integer: -32K ~ 32K, never overflows
  // sizeof(shot) * 6 * 256 = 3 KB
  __shared__ short E_shared[NBETA_PER_WORD][TperB];
  // sizeof(int) * 32 = 128 B
  __shared__ int __align__(32) E[NBETA];	// reduction of E_shared


  /// lattice scratchpad
  // sizeof(u_int32_t) * 16 * 16 * 16 = 16 KB
  __shared__ DATATYPE l[SZ][SZ][SZ];


  // index for global/shared memory transactions
  // temporority fix the code for BDx=8, BDy=16, BDz=2
  const int yy = (SZ_HF * threadIdx.z) + (threadIdx.y >> 1);
  const int xx = (SZ_HF * (threadIdx.y & 1)) + threadIdx.x;

  // index for reading scratch pad
  const int y = threadIdx.y;
  const int ya = (y + SZ - 1) % SZ;
  const int yb = (y + 1) % SZ;




  for (int word = 0; word < NWORD; word++) {

    // lattice scratchpad import
    int lattice_offset = SZ_CUBE * NWORD * blockIdx.x + SZ_CUBE * word + bidx;
    for (int zz = 0; zz < SZ; zz++)
      l[zz][yy][xx] = plattice[lattice_offset + (SZ * SZ * zz)];
    
    // reset partial status
    for (int b = 0; b < NBETA_PER_WORD; b++)
      E_shared[b][bidx] = 0;
    
    __syncthreads ();
    


    // updating a bipartite lattice consists two phases
    for (int run = 0; run < 2; run++) {
      int x0 = (threadIdx.z & 1) ^ (threadIdx.y & 1) ^ run;	// initial x
      int x = (threadIdx.x << 1) + x0;
      int xa = (x + SZ - 1) % SZ;
      int xb = (x + 1) % SZ;

      for (int zi = 0; zi < ZITER; zi++) {
	int z = ZITER * threadIdx.z + zi;

	DATATYPE center = l[z][y][x];	// center
	DATATYPE e0 = l[z][y][xa];	// left
	DATATYPE e1 = l[z][y][xb];	// right
	DATATYPE e2 = l[z][ya][x];	// up
	DATATYPE e3 = l[z][yb][x];	// down
	DATATYPE e4 = l[(z + SZ - 1) % SZ][y][x];	// front
	DATATYPE e5 = l[(z + 1) % SZ][y][x];	// back

	e0 = MASK_S * ((center >> SHIFT_J0) & 1) ^ e0 ^ center;
	e1 = MASK_S * ((center >> SHIFT_J1) & 1) ^ e1 ^ center;
	e2 = MASK_S * ((center >> SHIFT_J2) & 1) ^ e2 ^ center;
	e3 = MASK_S * ((center >> SHIFT_J3) & 1) ^ e3 ^ center;
	e4 = MASK_S * ((center >> SHIFT_J4) & 1) ^ e4 ^ center;
	e5 = MASK_S * ((center >> SHIFT_J5) & 1) ^ e5 ^ center;

	/*
	// sum J S S = S sum J S

	// process a lattice element
	for (int b = 0; b < NBETA_PER_WORD; b++) {
	  int shift =  MSCT * b;
	  int energy =		// range: [0, 6]
	    (e0 >> shift & 1) +
	    (e1 >> shift & 1) +
	    (e2 >> shift & 1) +
	    (e3 >> shift & 1) +
	    (e4 >> shift & 1) +
	    (e5 >> shift & 1);
	  // energy = (energy * 2 - 6) / 2;
	  energy = energy - 3;
	  E_shared[b][bidx] += energy;
	}
	*/



#if MSCT != 1
	  e0 = e0 + e1 + e2 + e3 + e4 + e5;
#endif

	  //#pragma unroll
	  for (int b = 0; b < NBETA_PER_WORD; b++) {
#if MSCT == 1
	    int energy =		// range: [0,6]
	      (e0 >> b & MASK_E) +
	      (e1 >> b & MASK_E) +
	      (e2 >> b & MASK_E) +
	      (e3 >> b & MASK_E) +
	      (e4 >> b & MASK_E) +
	      (e5 >> b & MASK_E);
#else
	    int shift = MSCT * b;
	    int energy = (e0 >> shift) & MASK_E;	// range: [0,6]
#endif

	    energy = energy - 3;
	    E_shared[b][bidx] += energy;
	  }



      }				// end of "for zi"
      __syncthreads ();
    }				// end of "for run"
    gpu_reduction (E_shared, E, bidx, word);
    __syncthreads ();
  }				// end of "for word"



  gpu_shuffle (temp_idx_shared, temp_beta_shared, E, bidx, mod);

  // copy from shared memory to global memory
  // betas
  if (bidx < NBETA) {
    temp_idx_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].idx;
    temp_beta_shared[bidx] = ptemp[NBETA_MAX * blockIdx.x + bidx].beta;
  }

}













void
host_compute_q (DATATYPE * lattice1, St * st)
{
  DATATYPE xor_word, xor_bit;

  for (int pair = 0; pair < GD; pair += 2) {	// for every pair

    int tmp[NBETA];
    for (int i = 0; i < NBETA; i++)
      tmp[i] = 0;

    for (int i = 0; i < SZ_CUBE; i++) {	// for every element
      xor_word =
	lattice1[SZ_CUBE * pair + i] ^ lattice1[SZ_CUBE * (pair + 1) + i];
      for (int b = 0; b < NBETA; b++) {
	xor_bit = (xor_word >> b) & 0x1;
	tmp[b] += 1 - (xor_bit << 1);	// parallel: +1, reverse: -1
      }
    }

    for (int b = 0; b < NBETA; b++)
      st->q[pair / 2][b] = tmp[b];
  }

}




void
host_save_q (St * st, char *mydir, int rec)
{
  //printf ("write outputs to\t %s/\n", mydir);

  //for (int rec = 0; rec < Q_REC; rec++) {

  char myfile[STR_LENG];
  sprintf (myfile, "%s/q_%04d.txt", mydir, rec);
  FILE *fp = fopen (myfile, "a");
  if (fp == NULL) {
    fprintf (stderr, "failed to open %s\n", myfile);
    exit (1);
  }

  for (int pair = 0; pair < GD; pair += 2) {	// for every pair
    for (int b = 0; b < NBETA; b++)
      fprintf (fp, "%6d", st->q[pair / 2][b]);
    fprintf (fp, "\n");
  }

  fclose (fp);
}








#ifdef DEBUG_PRINT_PROB
    if (blockIdx.x == 0 && bidx == 0)
      for (int b = 0; b < NBETA_PER_WORD; b++) {
	printf ("prob: b=%02d %f, %d \t",
		NBETA_PER_WORD * word + b,
		ptemp[NBETA_MAX * blockIdx.x + NBETA_PER_WORD * word + b].beta,
		ptemp[NBETA_MAX * blockIdx.x + NBETA_PER_WORD * word + b].idx);
	for (int i = 0; i < NPROBZ_MAX; i++) {
	  printf ("%f ", temp_prob_shared[b][i]);
	}
	printf ("\n");
      }
#endif

