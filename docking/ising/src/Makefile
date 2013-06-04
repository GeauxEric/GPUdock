.SUFFIXES: .cu .cuh .h

#INC = -I.
LIBDIR = -L/usr/local/cuda-4.2/lib64
#LIBDIR = -L/usr/local/cuda/4.1.28/lib64
#LIBDIR = -L/usr/local/packages/cuda/4.0.17/cuda/lib64
LIB = -lcudart #-lcutil_x86_64
#-lm

CC = gcc
CPP = g++
CXX = nvcc
MPICC = mpicc

# -v
ARCHFLAG = -gencode arch=compute_20,code=sm_20
OPTFLAG0 = -O0 -g -G
OPTFLAG1 = -O0
OPTFLAG2 = -O2
OPTFLAG3 = -O3
ASMFLAG = -ptx
PROFFLAG = --ptxas-options=-v -keep
#WFLAGS = -Wall # not supported by nvcc

CFLAGS = -std=c99 -pthread
CPPFLAGS = -pthread
CXXFLAGS = $(INC) $(LIB) $(LIBDIR) $(ARCHFLAG)

srcdir = src
gpusrc = gpu.cu threads.cu
cpusrc = host_main.cpp host_func.cpp host_launcher.cu lcg.cpp
exec = ising
#exec = mpi_ising postprocess_ea postprocess_ea2 postprocess_q postprocess_qq postprocess_r12

default: $(exec) 


#.cu: Makefile
#	$(CXX) $(CFLAGS) -o $@ $*.cu
#


%.o: %.c
	$(CC) $(OPTFLAG3) $(CFLAGS) -c $<

%.o: %.cpp
	$(CPP) $(OPTFLAG3) $(CPPFLAGS) -c $<

%.o: %.cu
	$(CXX) $(OPTFLAG2) $(CXXFLAGS) -c $<

mpiprocess.o: mpiprocess.cpp
	$(MPICC) -c $<

mpi_ising: host_main.o mpiprocess.o threads.o host_func.o host_launcher.o lcg.o gpu.o host_kernel.o
	$(MPICC) $(OPTFLAG2) $(LIB) $(LIBDIR) -o $@ $^

ising: host_main.o threads.o host_func.o host_launcher.o lcg.o gpu.o host_kernel.o
	$(CXX) $(OPTFLAG2) $(CXXFLAGS) -o $@ $^

postprocess_%: postprocess_%.c
	$(CPP) -O2 $(CPPFLAGS) -o $@ $<



prof: $(gpusrc) $(cpusrc)
	$(CXX) $(OPTFLAG1) $(CXXFLAGS) $(PROFFLAG) $^

profclean: $(gpusrc) $(cpusrc)
	$(CXX) $(OPTFLAG1) $(CXXFLAGS) $(PROFFLAG) -clean $^

g_ising: $(gpusrc) $(cpusrc)
	$(CXX) $(OPTFLAG0) $(CXXFLAGS) -o $@ $^

#ptx: $(src)
#	$(CXX) $(OPTFLAG1) $(CXXFLAGS) $(ASMFLAG) $^



install:
	cp $(exec) /work/sfeng/ising/l16_test/

clean:
	rm -r *.o $(exec)

cleanoutput:
	rm -r output_*

