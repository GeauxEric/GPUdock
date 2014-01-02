#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "util.h"
#include "size.h"
#include "hdf5io.h"

#define N_REP 240


int
main (int argc, char **argv)
{
  if (argc != 2) {
    fprintf (stderr, "usage: %s <input file>\n", argv[0]);
  }

  LigRecord *ligrecord;
  const size_t ligrecord_sz = sizeof (LigRecord) * N_REP;
  ligrecord = (LigRecord *) malloc (ligrecord_sz);
  
  ReadLigRecord (ligrecord, N_REP, argv[1]);

  const int myreplica = 0;
  const int repp_begin = 0;
  const int repp_end = 22;
  const int iter_begin = 0;
  //const int iter_end = STEPS_PER_DUMP - 1;
  const int iter_end = minimal_int (STEPS_PER_DUMP, 20) - 1;
  const int arg = 2;

  ComplexSize complexsize;
  complexsize.n_prt = 3;
  complexsize.n_tmp = MAXTMP;
  complexsize.n_lig = 20;
  complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;







  PrintLigRecord (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);
  //PrintRepRecord (ligrecord, STEPS_PER_DUMP, repp_begin, repp_end, iter_begin, iter_end, arg);
  PrintRepRecord2 (ligrecord, complexsize, STEPS_PER_DUMP, 0, 0, iter_begin, iter_end, arg);
  //PrintMoveRecord (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);

  free (ligrecord);
  return 0;
}
