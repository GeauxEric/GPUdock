
#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )


inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err )
    {
      fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
	       file, line, cudaGetErrorString( err ) );
      exit( -1 );
    }
#endif

  return;
}


inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err )
    {
      fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
	       file, line, cudaGetErrorString( err ) );
      exit( -1 );
    }

  // More careful checking. However, this will affect performance.
  // Comment away if needed.
  err = cudaDeviceSynchronize();
  if( cudaSuccess != err )
    {
      fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
	       file, line, cudaGetErrorString( err ) );
      exit( -1 );
    }
#endif

  return;
}						\




/*
    cudaError_t error = cudaGetLastError ();
    if (error != cudaSuccess) {
      printf ("CUDA error: %s\n", cudaGetErrorString (error));
      exit (-1);
    }




# define CUDA_CALL ( x ) do { if (( x ) != cudaSuccess ) { \
      printf (" Error at % s :% d \ n " , __FILE__ , __LINE__ ) ; \
      return EXIT_FAILURE ;}} while (0)

# define CURAND_CALL ( x ) do { if (( x ) != CURAND_STATUS_SUCCESS ) { \
      printf (" Error at % s :% d \ n " , __FILE__ , __LINE__ ) ; \
      return EXIT_FAILURE ;}} while (0)

*/

