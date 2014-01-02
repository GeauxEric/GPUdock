
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

    
    //__shared__ int idx_tmp[MAXREP];
    //__shared__ int idx_lig[MAXREP];
    

    //__shared__ int common_temp_replica[n_prt_dc][n_lig_dc][n_tmp_dc];
    //__shared__ int common_temp_replica[n_tmp_dc];
    __shared__ int common_temp_replica[MAXTMP];
    __shared__ int common_temp_order[MAXTMP];


    if (bidx == 0) {

      for (int p = 0; p < n_prt_dc; ++p) {
	for (int l = 0; l < n_lig_dc; ++l) {

	  // copy index from the replica structure
	  for (int t = 0; t < n_tmp_dc; ++t) {
	     const int linear_addr =
	       n_tmp_dc * n_lig_dc * p + n_lig_dc * t + l;
	     common_temp_replica[t] = replica_dc[linear_addr].idx_tmp;
	     common_temp_order[common_temp_replica[t]] = t;
	  }

	  /*
	  // bubble sort temperature in increment order
	  for (int i = 0; i < n_tmp_dc - 1; ++i) {
	    for (int j = i + 1; j < n_tmp_dc; ++j) {
	      const int left = common_temp_replica[i];
	      const int right = common_temp_replica[j];
	      if (left > right) {
		common_temp_replica[i] = right;
		common_temp_replica[j] = left;
	      }
            }
	  }
	  */


          // exchange temperature, assume even number of temperatures
	  for (int t = mode_t; t < n_tmp_dc - 1; t += 2) {
	    int left_t = common_temp_order[t];
	    int right_t = common_temp_order[t + 1];
	    const int tt = common_temp_replica[left_t];
	    common_temp_replica[left_t] = common_temp_replica[right_t];
	    common_temp_replica[right_t] = tt;
	  }

 
	  // copy index to the replica structure
	  for (int t = 0; t < n_tmp_dc; ++t) {
	    const int linear_addr =
	      n_tmp_dc * n_lig_dc * p + n_lig_dc * t + l;
	    replica_dc[linear_addr].idx_tmp = common_temp_replica[t];
	  }

        }
      }


    }





/*

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

*/

  }


}





