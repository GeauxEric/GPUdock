

  int i = 0;
  cudaSetDevice (i);
  CUDAKERNEL (MonteCarlo_Init_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

  for (int s1 = 0; s1 < mcpara->steps_total; s1 += mcpara->steps_per_dump) {

    double t0 = host_time_now ();
    printf ("\t%d / %d \r", s1, mcpara->steps_total);
    fflush (stdout);

    for (int s2 = 0; s2 < mcpara->steps_per_dump; s2 += mcpara->steps_per_exchange) {
      CUDAKERNEL (MonteCarlo_d, dim_grid, dim_block, rep_begin[i], rep_end[i], s1, s2);
# if IS_EXCHANGE == 1
      const int mode_l = 4; // ligand exchange mode
      const int mode_t = (s2 / mcpara->steps_per_exchange) % 2; // change exchange mode alternatively
      CUDAKERNEL (ExchangeReplicas_d, dim_grid, dim_block, mode_l, mode_t);
# endif
    }

    // accumulate for compute time
    mclog->t0 += host_time_now () - t0;

#if IS_OUTPUT == 1
    // copy ligand record from GPU to CPU memory
    CUDAMEMCPY (&ligrecord[rep_begin[i]], ligrecord_d[i], ligrecord_sz_per_gpu[i], cudaMemcpyDeviceToHost);

    // dump ligand record from CPU memory to disk
    char myoutputfile[MAXSTRINGLENG];
    sprintf(myoutputfile, "%s/%s_%04d.h5", mcpara->outputdir, mcpara->outputfile, s1 / mcpara->steps_total);
    DumpLigRecord (ligrecord, n_rep, myoutputfile);
#endif
    
    // accumulate for wall time (compute time plus I/O time)
    mclog->t1 += host_time_now () - t0;
  }

