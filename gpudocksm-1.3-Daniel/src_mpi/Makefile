EXE = mpi_mcc

OBJ = mpi_mcc.o coords.o data.o complex.o rmsd.o mpidecoys.o
SRC_MPI = mpidecoys.C
SRC_CPU = mpi_mcc.C coords.C data.C complex.C rmsd.f 

#CPP = g++
#FC = gfortran
#MPICC = /usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpic++

CPP = /usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpicxx
FC =  /usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpif90
MPICC = /usr/local/packages/mpich/3.0.2/Intel-13.0.0/bin/mpicc


CFLAGS = -O2 -Wall -I/usr/local/packages/gsl/1.15/Intel-13.0.0/include -I.
FFLAGS = -O2 -Wall
LDFLAGS = -L/usr/local/packages/gsl/1.15/Intel-13.0.0/lib -lgsl -lgslcblas -lm


default: $(EXE)
all: $(EXE)

$(EXE): $(OBJ)
## $(MPICC) $(LDFLAGS) -o $@ $(OBJ)
		$(MPICC) -o $@ $(OBJ) $(LDFLAGS)
		@cp $(EXE) ../bin/
		@cp $(EXE) ../test/
		@cp $(EXE) ~/bin/
		@cp $(EXE) /work/jaydy/bin/
		@echo "cp $(EXE) ../bin/"
		@echo "cp $(EXE) ../test/"
		@echo "cp $(EXE) ~/bin/"
		@echo "cp $(EXE) /work/jaydy/bin/"

# highest priority
mpidecoys.o: mpidecoys.C
		$(MPICC) -c $<

%.o: %.f
		$(FC) $(FFLAGS) -c $<

%.o: %.C
		$(CPP) $(CFLAGS) -c $<

clean:
		@(rm -f ${EXE} ../bin/hpcmcc ${OBJ} )
