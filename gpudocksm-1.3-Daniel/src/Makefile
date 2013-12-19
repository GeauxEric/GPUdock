
EXE = hpcmcc

OBJ = hpcmcc.o coords.o data.o complex.o mpidecoys.o rmsd.o

CPP = g++

FC = gfortran

CFLAGS = -O2 -Wall -I/project/michal/apps/gsl/include/ -I.

FFLAGS = -O2 -Wall -static

LDFLAGS = -lm -lgsl -lgslcblas -L/project/michal/apps/gsl/lib/

default: $(EXE)

all: $(EXE)

hpcmcc: $(OBJ)
	$(CPP) -o $@ $(OBJ) $(LDFLAGS)
	@mv hpcmcc ../bin/

hpcmcc.o: hpcmcc.C
	$(CPP) $(CFLAGS) -c -o hpcmcc.o hpcmcc.C

coords.o: coords.C
	$(CPP) $(CFLAGS) -c -o coords.o coords.C

data.o: data.C
	$(CPP) $(CFLAGS) -c -o data.o data.C

complex.o: complex.C
	$(CPP) $(CFLAGS) -c -o complex.o complex.C

rmsd.o: rmsd.f
	$(FC) $(FFLAGS) -c -o rmsd.o rmsd.f

mpidecoys.o: mpidecoys.C
	$(CPP) $(CFLAGS) -c -o mpidecoys.o mpidecoys.C

clean:
	@(rm -f ${EXE} ../bin/hpcmcc ${OBJ})

