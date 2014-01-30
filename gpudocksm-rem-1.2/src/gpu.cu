#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include "dock.h"
#include "toggle.h"
#include "include_cuda/cutil_inline.h"
#include "gpu.cuh"
#include "size_gpu.cuh"


// _dc stands for gpu constant memory

// array pointers
__constant__ Protein *prt_dc;
__constant__ Psp *psp_dc;
__constant__ Kde *kde_dc;
__constant__ Mcs *mcs_dc;
__constant__ EnePara *enepara_dc;
__constant__ Temp *temp_dc;

__constant__ Ligand *lig_dc;
__constant__ Replica *replica_dc;
__constant__ float *etotal_dc;
__constant__ LigMoveVector *ligmovevector_dc;
__constant__ LigRecord *ligrecord_dc;
__constant__ TmpEnergy *tmpenergy_dc;
__constant__ int *acs_mc_dc;
__constant__ int *acs_temp_exchg_dc;
__constant__ ConfusionMatrix *ref_matrix_dc;

// PRNG seeds
__constant__ int seed_dc;
__constant__ curandState *curandstate_dc;



// monte carlo parameters
__constant__ int steps_per_exchange_dc;
__constant__ int steps_per_dump_dc;
__constant__ int steps_total_dc;
__constant__ int is_random_dc;
__constant__ float * move_scale_dc;


__constant__ float enepara_lj0_dc;
__constant__ float enepara_lj1_dc;
__constant__ float enepara_el0_dc;
__constant__ float enepara_el1_dc;
__constant__ float enepara_a1_dc;
__constant__ float enepara_b1_dc;
__constant__ float enepara_kde2_dc;
__constant__ float enepara_kde3_dc;



// residue numbers (of per replica)
__constant__ int lna_dc;
__constant__ int pnp_dc;
__constant__ int pnk_dc;
__constant__ int n_pos_dc;

// replica numbers
__constant__ int n_lig_dc;
__constant__ int n_prt_dc;
__constant__ int n_tmp_dc;
__constant__ int n_rep_dc;


#include "initcurand_d.cu"
#include "exchangereplicas_d.cu"
#include "montecarlo_d.cu"
#include "move_d.cu"
#include "calcenergy_d.cu"
#include "accept_d.cu"
#include "util_d.cu"
#include "calcmcc_d.cu"

