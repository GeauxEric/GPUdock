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


EXE = decoy_gen

OBJ = decoy_gen.o coords.o data.o complex.o walltime.o rmsd.o remc.o montecarlo.o

CPP = g++

FC = gfortran

#CFLAGS = -O2 -Wall -I/project/michal/apps/gsl/include/ -I.
CFLAGS = -O2 -Wall -I/usr/local/packages/gsl/1.15/Intel-13.0.0/include/ -I.

FFLAGS = -O2 -Wall -static

#LDFLAGS = -lm -lgsl -lgslcblas -L/project/michal/apps/gsl/lib/
LDFLAGS = -lm -lgsl -lgslcblas -L/usr/local/packages/gsl/1.15/Intel-13.0.0/lib/

default: $(EXE)

all: $(EXE)

$(EXE): $(OBJ)
	$(CPP) -o $@ $(OBJ) $(LDFLAGS)
	@cp $(EXE) ../test/
	@echo "cp $(EXE) ../test/"
	@cp $(EXE) ../bin/
	@echo "cp $(EXE) ../bin/"
	@cp $(EXE) ~/bin/
	@echo "cp $(EXE) ~/bin/"
	@cp $(EXE) ~/work/bin/
	@echo "cp $(EXE) ~/work/bin/"


%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.C
	$(CPP) $(CFLAGS) -c $<

# decoy_gen_gapless.o: decoy_gen_gapless.C
# 	$(CPP) $(CFLAGS) -c -o decoy_gen_gapless.o decoy_gen_gapless.C
# 
# coords.o: coords.C
# 	$(CPP) $(CFLAGS) -c -o coords.o coords.C
# 
# data.o: data.C
# 	$(CPP) $(CFLAGS) -c -o data.o data.C
# 
# complex.o: complex.C
# 	$(CPP) $(CFLAGS) -c -o complex.o complex.C
# 
# remc.o: remc.C
# 	$(CPP) $(CFLAGS) -c -o remc.o remc.C
# 
# montecarlo.o: montecarlo.C
# 	$(CPP) $(CFLAGS) -c -o montecarlo.o montecarlo.C
# 
# rmsd.o: rmsd.f
# 	$(FC) $(FFLAGS) -c -o rmsd.o rmsd.f
# 
# walltime.o: walltime.C
# 	$(CPP) $(CFLAGS) -c -o walltime.o walltime.C

clean:
	@(rm -f ${EXE} ../bin/decoy_gen_gapless ${OBJ})

