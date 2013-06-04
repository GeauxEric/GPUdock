#ifndef TABLE_CUH
#define TABLE_CUH


#include <stdlib.h>
#include <stdio.h>
//#include <stdint.h>
#include <math.h>
#include <time.h>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
//#include <cutil.h>
#include "../mycudacall.cuh"




// table size
#define SZ (20 * 1024 * 1024)

// blocksPerGrid
#define GD 16
// threadsPerBlock
#define BD 64

#define NT (BD * GD)


// number of betas combined in a integer
#define NBETA 14
#define NBETA_ALL (NBETA * GD)

#define BETA_BEGIN 1.5
#define BETA_END 1.8



#define BITMASK_S 0xFC11111111111111
#define BITMASK_J 0xFC00000000000000



#endif
