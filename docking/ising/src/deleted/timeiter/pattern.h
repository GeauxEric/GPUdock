

// checkerboard decomposition
// bipartite lattice, two phases update
//     pro
//       no storage overhead

// time iteration
// odd-even instance
// integrate two lattice of new/old instance using reserved bits of MSCT
//     pro
//       "adjacent access pattern", better shared memory BW utilizations
//       no wrap around variables, no "x0", straight forward index generation
//       if update multiple integers in every iteration, more data reuse
//     con
//     	 extra operations to distinguish instance
//       iteration must be even

#define CHECKERBOARD 0
#define TIMEITER 1




#define PATTERN CHECKERBOARD
//#define PATTERN TIMEITER

