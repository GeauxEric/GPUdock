#ifndef _CUTIL_YEAH_
#define _CUTIL_YEAH_

#include <cutil.h>

#define CUDAVERSION 4


/*

1.
CheckError (cudaMemcpyToSymbol (...));
CheckError2 (cudaMemcpyToSymbol (...), "Error message");

2.
kernel <<< ... >>> (...);
cutilCheckMsgAndSync ("kernel failed\n");


// call this function before the first kernal invocation
cudaGetLastError ();

*/


inline void CheckError2(int err, char* pMsg)
{
  if (err != 0) {
    printf(" Error on %s\n" , pMsg);
    cutilCheckMsg(pMsg);
    exit(-1);
  }
}


inline void CheckError(int err)
{
  if (err != 0) {
    printf("Error!\n");
    exit(-1);
  }
}


#define CUDAKERNEL(funcname, dim_grid, dim_block, ...) \
  funcname <<< dim_grid, dim_block >>> (__VA_ARGS__); \
  cutilCheckMsgAndSync ("#funcname kernel failed\n")

#define CUDAKERNELSTREAMSYNC(funcname, dim_grid, dim_block, n, stream, ...) \
  funcname <<< dim_grid, dim_block, n, stream >>> (__VA_ARGS__);	\
  cutilCheckMsgAndSync ("#funcname kernel failed\n")

#define CUDAKERNELSTREAM(funcname, dim_grid, dim_block, n, stream, ...) \
  funcname <<< dim_grid, dim_block, n, stream >>> (__VA_ARGS__);	\
  cutilCheckMsg ("#funcname kernel failed\n")




#define CUDAMEMCPY(dst, src, sz, direction) \
  CheckError2 (cudaMemcpy (dst, src, sz, direction), #dst)



#if CUDAVERSION == 4
#define CUDAMEMCPYTOSYMBOL(dst, src, type) \
  CheckError2 (cudaMemcpyToSymbol (#dst, src, sizeof (type), 0, cudaMemcpyHostToDevice), #dst)
#elif CUDAVERSION ==5
#define CUDAMEMCPYTOSYMBOL(dst, src, type) \
  CheckError2 (cudaMemcpyToSymbol (dst, src, sizeof (type), 0, cudaMemcpyHostToDevice), #dst)
#else
  printf ("cuda version not supportted by this utility box\n");
#endif




#endif

