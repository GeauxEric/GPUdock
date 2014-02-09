EXE := bayes_dock
OBJ_CPU := dock.o load.o data.o rmsd.o util.o timing.o quicksort.o hdf5io.o 
OBJ_GPU := run.o

SH := sh
CPP := g++
CXX := nvcc
FC := gfortran

HEADPATH := -I. -Iinclude_cuda -Ih5
#LIBPATH := -L/usr/local/cuda/4.1.28/lib64
#LIBPATH := -L/usr/local/cuda-4.2/lib64
#LIBPATH := -L/usr/local/cuda-4.2.9/lib64
#LIBPATH := -L/usr/local/cuda-5.5/lib64
#LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
LIBPATH := -L/usr/local/packages/cuda/4.2.9/cuda/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
LINKFLAGS := -lm -lcudart -lhdf5 #-lcudadevrt

ARCHFLAG := -gencode arch=compute_20,code=sm_20
#ARCHFLAG := -gencode=arch=compute_20,code=\"sm_20,compute_20\"
#ARCHFLAG := -gencode arch=compute_30,code=sm_30
#ARCHFLAG := -gencode arch=compute_35,code=\"sm_35,compute_35\"
#ARCHFLAG := -arch=sm_35

OPTFLAGS := -O2
#OPTFLAGS := -O3


CPPFLAGS := -Wall $(HEADPATH) $(OPTFLAGS)
CXXFLAGS := $(HEADPATH) $(OPTFLAGS) $(ARCHFLAG) # --ptxas-options=-v #-ccbin=g++  #-keep #-rdc=true
FFLAGS := -Wall -static $(OPTFLAGS)




default: $(EXE) analysis

$(EXE): $(OBJ_CPU) $(OBJ_GPU) 
	$(CPP) $(LINKFLAGS) $(LIBPATH) $(OBJ_CPU) $(OBJ_GPU) -o $@
	@cp $(EXE) ../bin/
	@cp $(EXE) ~/bin/
	@cp $(EXE) ../test/

hdf5io.o: hdf5io.C
	h5c++ -c $<

%.o: %.C
	$(CPP) $(CPPFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.cu
	$(CXX) $(CXXFLAGS) -c $<


count: count.o load.o data.o rmsd.o
	$(CPP) $(LINKFLAGS) -O0 -g $^ -o $@
	@cp count ../test/
	@echo "copy count to ../test/"

analysis: analysis.o util.o hdf5io.o
	$(CPP) $(LINKFLAGS) $(LIBPATH) $^ -o $@
	@cp analysis ../bin/
#	@cp analysis ~/bin/
	@cp analysis ../test/

clean:
	@(rm -f ${EXE} ${OBJ_CPU} ${OBJ_GPU} analysis count *i *cubin* *.ptx *cudafe* *.fatbin* *.hash count count.o ../bin/dock)

