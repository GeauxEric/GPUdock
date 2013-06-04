# reference:
# LSU EE 7700-2 (Sp 2011), GPU Microarchitecture
# gp/hw/hw2

INC = -I/nfs/snow/yfang11/app/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc
LIB = -L/nfs/snow/yfang11/app/NVIDIA_GPU_Computing_SDK/C/lib -lcutil_x86_64

CC = gcc
CXX = nvcc

.SUFFIXES: .cu .cuh


# -v
#DEBUGFLAG = -g -G
DEBUGFLAG =
ARCHFLAG = -gencode arch=compute_20,code=sm_20
#OPTFLAG = -O0 -g -G
OPTFLAG = -O3

FLAGS = $(INC) $(ARCHFLAG) $(OPTFLAG) $(DEBUGFLAG)
CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)

default: table


#.cu: Makefile
#	$(CXX) $(CFLAGS) -o $@ $*.cu

ising: ising.cu ising.cuh
	$(CXX) $(CFLAGS) -o $@ $^

table: table.cu table.cuh
	$(CXX) $(CFLAGS) -o $@ $^


clean:
	rm ising table
