/*
#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "gpu.cuh"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
*/





__device__ void
InitAcs_d (const int bidx)
{
  if (blockIdx.x == 0) {
    for (int i = bidx; i < MAXREP; i += TperB) {
      acs_dc[i] = 0;
    }
  }
}

__device__ void
InitLigRecord_d (const int bidx, const int myreplica, const int rep_begin)
{
  for (int s2s3 = 0; s2s3 < steps_per_dump_dc; ++s2s3) {
    LigRecordSingleStep *myrecord = &ligrecord_dc[myreplica - rep_begin].step[s2s3];
    myrecord->replica.idx_rep = 0;
    myrecord->replica.idx_prt = 0;
    myrecord->replica.idx_tmp = 0;
    myrecord->replica.idx_lig = 0;

    for (int i = 0; i < MAXWEI; ++i)
      myrecord->energy.e[i] = 0.0f;

    for (int i = 0; i < 6; ++i)
      myrecord->movematrix[i] = 0.0f;

    myrecord->step = 0;
  }

}


/*
__forceinline__
__device__ void
BackupLigCoord_d (const int bidx, Ligand *mylig)
{

  const LigCoord *src = &mylig->coord_old;
  LigCoord *dst = &mylig->coord_bkup;

  for (int atom = bidx; atom < lna_dc; atom += TperB) {
    dst->x[atom] = src->x[atom];
    dst->y[atom] = src->y[atom];
    dst->z[atom] = src->z[atom];
  }
  if (bidx < 3)
    dst->center[bidx] = src->center[bidx];

}
*/



__device__ void
ComputeMoveMatrix_d (const int bidx, const int myreplica, Ligand *mylig)
{
  if (bidx < 6)
    ligmovematrix_dc[myreplica].a[bidx] = mylig->movematrix_old[bidx];

}






__device__ void
RecordLigand_d (const int bidx, const int s1, const int s2s3, const int myreplica, const int rep_begin, const Ligand * mylig)
{
  /*
  if (bidx == 0) // && myreplica == 0)
    printf ("rep %d, iter %d, rep_begin %d, n_rep %d, idx_rep %d\n",
	    myreplica, s2, rep_begin, n_rep_dc, replica_dc[myreplica].idx_rep);

  if (myreplica == 0) {
    PrintEnergy2_d (bidx, mylig, myreplica, s1 + s2s3, 2);
  }
  */

  if (bidx == 0) {
    LigRecordSingleStep *myrecord = &ligrecord_dc[myreplica - rep_begin].step[s2s3];
    myrecord->replica = replica_dc[myreplica];
    myrecord->energy = mylig->energy_old;
    for (int i = 0; i < 6; ++i)
      myrecord->movematrix[i] = mylig->movematrix_old[i];
    myrecord->step = s1 + s2s3;
  }

}






__forceinline__
__device__ float
MyRand_d ()
{

  float randdd; 

  if (is_random_dc == 0) {
     randdd = 20.0f;
    //randdd = 0.0f;
  }
  else {
    const int gidx =
      blockDim.x * blockDim.y * blockIdx.x +
      blockDim.x * threadIdx.y +
      threadIdx.x;
    curandState myseed = curandstate_dc[gidx];
    randdd = curand_uniform (&myseed);
    curandstate_dc[gidx] = myseed;
  }

  // printf("%f\n", randdd);

  return randdd;
}




/*
__forceinline__
__device__ int
Mininal_int_d (const int a, const int b)
{
 return a < b ? a : b;
}
*/





__forceinline__
__device__ void
SumReduction1D_d (const int bidx, float *a)
{
  __syncthreads ();

  for (int stride = TperB / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride)
      a[bidx] += a[stride + bidx];
    __syncthreads ();
  }
}


__forceinline__
__device__ void
SumReduction1D_5_d (const int bidx, float *a, float *b, float *c, float *d, float *e)
{
  __syncthreads ();

  for (int stride = TperB / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride) {
      a[bidx] += a[stride + bidx];
      b[bidx] += b[stride + bidx];
      c[bidx] += c[stride + bidx];
      d[bidx] += d[stride + bidx];
      e[bidx] += e[stride + bidx];
    }
    __syncthreads ();
  }
}


__forceinline__
__device__ void
SumReduction2D_d (float a[BDy][BDx])
{
  __syncthreads ();

  for (int stride = BDx / 2; stride >= 1; stride >>= 1) {
    if (threadIdx.x < stride) {
      a[threadIdx.y][threadIdx.x] += a[threadIdx.y][stride + threadIdx.x];
    }
    __syncthreads ();
  }
}


__forceinline__
__device__ void
SumReduction2D_2_d (float a[BDy][BDx], float b[BDy][BDx])
{
  __syncthreads ();

  for (int stride = BDx / 2; stride >= 1; stride >>= 1) {
    if (threadIdx.x < stride) {
      a[threadIdx.y][threadIdx.x] += a[threadIdx.y][stride + threadIdx.x];
      b[threadIdx.y][threadIdx.x] += b[threadIdx.y][stride + threadIdx.x];
    }
    __syncthreads ();
  }
}



