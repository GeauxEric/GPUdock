
// mode 0: 0 swap 1 , 2 swap 3 , 4 swap 5 , ...
// mode 1: 0 , 1 swap 2 , 3 swap 4 , ...
// mode 3: random
// mode 4: not exchange

__global__ void
ExchangeReplicas_d (const int mode_l, const int mode_t)
{

  if (blockIdx.x == 0) {
    const int bidx = blockDim.x * threadIdx.y + threadIdx.x;      // within a TB

    // should implement a function to check if MAX_XXX is in proper size
    // eg: assert ((MAXREP > n_rep) && (not exaust the shared memory))

    
    __shared__ int idx_tmp[MAXREP];
    __shared__ int idx_lig[MAXREP];

    for (int r = bidx; r < n_rep_dc; r += TperB) {
      idx_tmp[r] = replica_dc[r].idx_tmp;
      idx_lig[r] = replica_dc[r].idx_lig;
    }
    
    __syncthreads ();

    // exchange temperature
    if (mode_t < 2) {
      int bidx_max = n_rep_dc / 2 + n_rep_dc % 2 - mode_t;
      int left_t = (bidx << 1) + mode_t;
      int right_t = left_t + 1;
      if (bidx < bidx_max) {
	int t = idx_tmp[left_t];
	idx_tmp[left_t] = idx_tmp[right_t];
	idx_tmp[right_t] = t;
      }
    }



    // exchange ligand
    if (mode_l < 2) {
      int bidx_max = n_rep_dc / 2 + n_rep_dc % 2 - mode_t;
      int left_l = bidx << 1 + mode_l;
      int right_l = left_l + 1;
      if (bidx < bidx_max) {
	int l = idx_lig[left_l];
	idx_lig[left_l] = idx_tmp[right_l];
	idx_lig[right_l] = l;
      }
    }
 

    __syncthreads ();

    for (int r = bidx; r < n_rep_dc; r += TperB) {
      replica_dc[r].idx_tmp = idx_tmp[r];
      replica_dc[r].idx_lig = idx_lig[r];
    }


  }


}



