EXE := dock
OBJ_CPU := dock.o load.o data.o rmsd.o util.o hdf5io.o
OBJ_GPU := run.o


SH := sh
CPP_HOST := g++
CPP_DEV := nvcc
FC := gfortran


HOSTFLAGS :=
DEVFLAGS := -gencode arch=compute_20,code=sm_20
#DEVFLAGS := -gencode=arch=compute_20,code=\"sm_20,compute_20\"
#DEVFLAGS := -gencode arch=compute_35,code=\"sm_35,compute_35\"
#DEVFLAGS := -arch=sm_35
#-ccbin=g++ # --ptxas-options=-v #-keep #-rdc=true


OPTFLAGS := -O3

HEADPATH := -I./include -I./include/nvidia_gpucomputingsdk_4.2.9_c_common_inc

#LIBPATH := -L/usr/local/packages/cuda/4.2.9/cuda/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
#LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.9/Intel-13.0.0-openmpi-1.6.2/lib

LINKFLAGS := -lcudart -lhdf5 -lm


HOSTFLAGS += -Wall $(OPTFLAGS) $(HEADPATH)
DEVFLAGS += $(OPTFLAGS) $(HEADPATH)
FFLAGS := $(OPTFLAGS)







default: $(EXE)

$(EXE): $(OBJ_CPU) $(OBJ_GPU) 
	$(CPP_HOST) $(HOSTFLAGS) $(LIBPATH) $(OBJ_CPU) $(OBJ_GPU) -o $@ $(LINKFLAGS)

hdf5io.o: hdf5io.C
	h5c++ -c $<

%.o: %.C
	$(CPP_HOST) $(HOSTFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.cu
	$(CPP_DEV) $(DEVFLAGS) -c $<

install:
	@cp $(EXE) ../work/

clean:
	@(rm -f ${EXE} ${OBJ_CPU} ${OBJ_GPU} count *i *cubin* *.ptx *cudafe* *.fatbin* *.hash count count.o ../bin/dock)


