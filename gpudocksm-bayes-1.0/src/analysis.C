#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "dock.h"
#include "util.h"
#include "size.h"
#include "hdf5io.h"

using namespace std;


int
main (int argc, char **argv)
{
  if (argc < 2) {
    fprintf (stderr, "usage: %s  <input file>\n", argv[0]);
    printf ("-nl <number of show lines>\n");
    printf ("-l <ligand conf number>\n");
    printf ("-p <protein conf number>\n");
  }

  // default settings
  int lig_conf = 0;
  int prt_conf = 0;
  int show_energy = 0;
  int show_rep = 0;
  int myreplica = 0;

  ComplexSize complexsize;

  complexsize.n_prt = 11;
  complexsize.n_tmp = MAXTMP;
  complexsize.n_lig = 31;
  complexsize.lna = 0; // unused, the value does not matter
  complexsize.pnp = 0; // unused, the value does not matter
  complexsize.pnk = 0; // unused, the value does not matter
  complexsize.pos = 0; // unused, the value does not matter

  // ./analysis -num_temp 20 -rep 589 -nl 4500 -l 1 -p 2 -e XXX.h5

  for ( int i = 0; i < argc; i++ ) {
    if ( !strcmp(argv[i],"-l")  && i < argc ) 
      lig_conf = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-e")  && i < argc ) 
      show_energy = 1;
    if ( !strcmp(argv[i],"-r")  && i < argc ) 
      show_rep = 1;
    if ( !strcmp(argv[i],"-p")  && i < argc ) 
      prt_conf = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-np")  && i < argc ) 
      complexsize.n_prt = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-nl")  && i < argc ) 
      complexsize.n_lig = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-nt")  && i < argc ) 
      complexsize.n_tmp = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-rep")  && i < argc ) 
      myreplica = atoi(argv[i+1]);
  }


  complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;

  	

  /*
  int lig_conf = 0;

  int prt_conf = 0;

  int replica_same_temp[MAXTMP];

  for (int t = 0; t < complexsize.n_tmp; t++) {
	int flatten_addr =
	       complexsize.n_tmp * complexsize.n_lig * prt_conf + complexsize.n_lig * t + lig_conf;
	// cout << "flatten addr :" << flatten_addr << endl;
	replica_same_temp[t] = flatten_addr;
	// cout << flatten_addr << endl;
  }
  */


  LigRecord *ligrecord;
  size_t ligrecord_sz = sizeof (LigRecord) * complexsize.n_rep;
  ligrecord = (LigRecord *) malloc (ligrecord_sz);
  ReadLigRecord (ligrecord, complexsize.n_rep, argv[argc-1]);

  int iter_begin = 0;
  int iter_end = STEPS_PER_DUMP - 1;
  int arg = 2;

  PrintTrack (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);

  if (show_energy == 1)
    PrintLigRecord (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);

  if (show_rep == 1)
    PrintRepRecord2 (ligrecord, complexsize, STEPS_PER_DUMP, 
		     lig_conf, prt_conf, 
		     iter_begin, iter_end, arg);

  free (ligrecord);
  return 0;
}

