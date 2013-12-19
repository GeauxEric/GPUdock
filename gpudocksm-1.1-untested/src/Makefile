#==============================================================================================
#     __________________ ____ ___              .___             __      _________   _____   
#    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
#   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
#   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
#    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
#           \/                                  \/           \/     \/        \/         \/ 
#
#   GPU-accelerated hybrid-resolution ligand docking using Replica Exchange Monte Carlo
#
#==============================================================================================


EXE = gpudocksm gpudocksm_lhm

OBJ = gpudocksm.o coords.o data.o complex.o walltime.o rmsd.o remc.o montecarlo.o random.o

SH = sh

CPP = g++

FC = gfortran

CFLAGS = -O2 -Wall -I../lib -I.

FFLAGS = -O2 -Wall -static

LDFLAGS = -lm -L../lib

default: $(EXE)

all: $(EXE)

gpudocksm: $(OBJ)
	$(CPP) -o $@ $(OBJ) $(LDFLAGS)
	@mv gpudocksm ../bin/

gpudocksm_lhm:
	$(SH) gpudocksm_lhm.shar
	@chmod +x gpudocksm_lhm
	@mv gpudocksm_lhm ../bin/

gpudocksm.o: gpudocksm.C
	$(CPP) $(CFLAGS) -c -o gpudocksm.o gpudocksm.C

coords.o: coords.C
	$(CPP) $(CFLAGS) -c -o coords.o coords.C

data.o: data.C
	$(CPP) $(CFLAGS) -c -o data.o data.C

complex.o: complex.C
	$(CPP) $(CFLAGS) -c -o complex.o complex.C

remc.o: remc.C
	$(CPP) $(CFLAGS) -c -o remc.o remc.C

montecarlo.o: montecarlo.C
	$(CPP) $(CFLAGS) -c -o montecarlo.o montecarlo.C

rmsd.o: rmsd.f
	$(FC) $(FFLAGS) -c -o rmsd.o rmsd.f

walltime.o: walltime.C
	$(CPP) $(CFLAGS) -c -o walltime.o walltime.C

random.o: random.C
	$(CPP) $(CFLAGS) -c -o random.o random.C

clean:
	@(rm -f ${EXE} ../bin/gpudocksm ../bin/gpudocksm_lhm ${OBJ})

