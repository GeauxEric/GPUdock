#include "COPYING"


#ifndef SIM_H
#define SIM_H


//#include <stdint.h>




#define PI 3.141592654f

// beta ~ 1/temperature
//#define BETA_LOW 0.05f
//#define BETA_HIGH 0.4f

#define BETA_LOW 0.2f
#define BETA_HIGH 0.25f

// external field
#define H 0.0f
//#define H 1.0f

// randomly initialize J
//#define RANDJ
#define RANDS


// using the same random number for every spins integrated in a world
//#define SHARERAND







// iteration parameters


// total monte-carlo sweeps: 2,000,000
// status samples:           absolutely no more than 10,000, tipically 1,000

// for one realization, each status sample consumes memory:
//   (NBETA_MAX * sizeof (int)) * 2 = 256B
// assume assign 32 realiztions on a GPU, memory for saving status:
//   256 * 32 * 10,000 = 64MB

// simulation time estimzation for 16 realizations, 16^3 cubic lattice
// NBETA * (16 ^ 3) * 16 * (2 * 10^6) * (50PS/spin) = 170 seconds


/*
#define ITER_WARMUP          40000
#define ITER_WARMUP_KERN     10000
#define ITER_WARMUP_KERNFUNC 200
#define ITER_SWAP            0
#define ITER_SWAP_KERN       0
#define ITER_SWAP_KERNFUNC   0
*/

/*
#define ITER_WARMUP          4000
#define ITER_WARMUP_KERN     1000
#define ITER_WARMUP_KERNFUNC 200
#define ITER_SWAP            4000
#define ITER_SWAP_KERN       1000
#define ITER_SWAP_KERNFUNC   10
*/


/*
#define ITER_WARMUP          8000
#define ITER_WARMUP_KERN     2000
#define ITER_WARMUP_KERNFUNC 200
#define ITER_SWAP            8000
#define ITER_SWAP_KERN       1000
#define ITER_SWAP_KERNFUNC   10
*/

///*
#define ITER_WARMUP          0
#define ITER_WARMUP_KERN     0
#define ITER_WARMUP_KERNFUNC 0
#define ITER_SWAP            40000
#define ITER_SWAP_KERN       1000
#define ITER_SWAP_KERNFUNC   10
//*/

/*
#define ITER_WARMUP          4000
#define ITER_WARMUP_KERN     400
#define ITER_WARMUP_KERNFUNC 100
#define ITER_SWAP            0
#define ITER_SWAP_KERN       0
#define ITER_SWAP_KERNFUNC   0
*/

/*
#define ITER_WARMUP          0
#define ITER_WARMUP_KERN     0
#define ITER_WARMUP_KERNFUNC 0
#define ITER_SWAP            1
#define ITER_SWAP_KERN       1
#define ITER_SWAP_KERNFUNC   1
*/









// lattice size
// SZ must be even
// SZz must divides SZ

#define SZ 16
#define SZz SZ

#define SZ_HF (SZ / 2)
#define SZ_CUBE (SZ * SZ * SZz)
#define SZ_CUBE_HF (SZ_CUBE / 2)
#define SZ_TILE (SZ * SZ * SZ)


// SM per GPU
// should implement GD = func (prop.multiProcessorCount);

// GD - blocksPerGrid, must be even
// BD - threadsPerBlock
// when modifing "GD", should also update "GD_HF",

#define GD 64
#define GD_HF 32

#define TperB 256
// checkerboard 3D block
#define BDx0 SZ_HF
#define BDy0 SZ
#define BDz0 2








//#define DEBUG0
//#define DEBUG1
//#define DEBUG2
//#define DEBUG3

//#define PRINT_E
//#define PRINT_PROB
//#define PRINT_SPEED
//#define NO_OUTPUT





// calculate probablities
// __expf (2 * energy * temp_beta_shared[b]);
// temp_prob_shared[16 * b + (energy << 1) + spin];
#define DIRECT_COMPUTE 0
#define TABLE_LOOKUP 1
//typedef float PROB_DATATYPE;
typedef u_int32_t PROB_DATATYPE;





// storage allocation for 2 sublattices
#define SHARED 0
#define SEPARATED 1
#define INTEGRATED 2


// bit packaging format
#define SPARSE 0
#define COMPACT 1





#define PROB_GEN TABLE_LOOKUP

#define ALLOCATION SHARED
//#define ALLOCATION SEPARATED
//#define ALLOCATION INTEGRATED

//#define DENSE SPARSE
#define DENSE COMPACT

// Multispin Coding, MSCT = 1, 3, 4
#define MSCT 4

// MSC encoding alternatives, MSC_FORMAT = 0, 1
// !!!!!!!!  1 IS NOT WORKING !!!!!!!!
#define MSC_FORMAT 0






// Multispin Coding for 32 bit unsighed integer
#include "bit32.h"

// string length for file names, etc.
#define STR_LENG 64


 /*
   probability = expf (2 * beta * (energy - H * spin))
   prob []

   index energy spin   (energy - H * spin)
   00    -6       -1    -6+H
   01    -6       +1    -6-H
   
   02    -4       -1    -4+H
   03    -4       +1    -4-H
   
   04    -2       -1    -2+H
   05    -2       +1    -2-H
   
   06     0       -1     0+H
   07     0       +1     0-H
   
   08     2       -1     2+H
   09     2       +1     2-H
   
   10     4       -1     4+H
   11     4       +1     4-H
   
   12     6       -1     6+H
   13     6       +1     6-H
 */
#define NPROB 14
#define NPROB_MAX 16





// s structure should better allign on 4 byte boundary

typedef struct
{
  int E[NBETA][GD];
  int M[NBETA][GD];
  float U[NBETA][GD];		// U = 1 - M4 / (3 * M2 * M2)
} Avrg;


typedef struct
{
  float Q0[NBETA][TperB];
  float Qk_real[NBETA][TperB];
  float Qk_imag[NBETA][TperB];
  float Qk2_real[NBETA][TperB];
  float Qk2_imag[NBETA][TperB];                                               
} Qk;


typedef struct
{
  //  float e[GD][NBETA_MAX];
  float q[GD_HF][NBETA_MAX];
  float qk_real[GD_HF][NBETA_MAX];
  float qk_imag[GD_HF][NBETA_MAX];
  float qk2_real[GD_HF][NBETA_MAX];
  float qk2_imag[GD_HF][NBETA_MAX];
} St;


typedef struct
{
  int idx;
  float beta;
} Temp;


typedef struct
{
  double start; // start time stamp
  double span; // accumulated time span
  double avrg; // time span average
  int n; // number of records
} Timing;




// host_func.cc
void host_timing_init (Timing * t, int n);
void host_timing (Timing * t, int idx, int mode);
double host_time_now ();
void host_init_J (MSC_DATATYPE * l);
void host_init_S (MSC_DATATYPE * l);
void host_init_lattice (MSC_DATATYPE * l);
void host_save_st (St * st, char *mydir, int node, int device);
void host_init_temp (Temp * temp, float beta_low, float beta_high);
void host_report_speed_title ();
void host_report_speed (double start, double stop, int iter, char *event);
void host_usage (char *bin);
void host_summary (float beta_low, float beta_high, char *mydir);
void host_makedir (char *mydir);

//host_launcher.cu
void host_launcher (float beta_low, float beta_high, char* mydir, int node, int device);

//host_kernel.cu
void host_kernel_warmup (MSC_DATATYPE * lattice, Temp * temp);
void host_kernel_swap (MSC_DATATYPE * lattice, Temp * temp, St * st, int rec);



#endif /* SIM_H */

