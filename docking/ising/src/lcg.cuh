#ifndef LCG_CUH
#define LCG_CUH

#include "lcg.h"


__device__ void gpu_lcg (LCG_DATATYPE *seed);
__device__ float gpu_lcg_f (LCG_DATATYPE *seed);


#endif
