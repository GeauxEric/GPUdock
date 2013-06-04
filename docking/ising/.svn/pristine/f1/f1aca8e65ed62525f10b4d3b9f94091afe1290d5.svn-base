.SUFFIXES: .cu .cuh .h

INC = -I.
LIB = -lcutil_x86_64

CC = gcc
CPP = g++
CXX = nvcc

# -v
DEBUGFLAG = 
ARCHFLAG = -gencode arch=compute_20,code=sm_20
OPTFLAG0 = -O0 -g -G
OPTFLAG1 = -O0
OPTFLAG2 = -O2
ASMFLAG = -ptx
PROFFLAG = --ptxas-options=-v -keep
#WFLAGS=-Wall # not supported by nvcc
LIBS= -lm

FLAGS = $(INC) $(ARCHFLAG) $(DEBUGFLAG) $(LIBS) #-v
CFLAGS =
CPPFLAGS =
CXXFLAGS = $(FLAGS)

srcdir = src
gpusrc = gpu.cu
cpusrc = host_main.cc host_func.cc host_launcher.cu

default: postprocess_e postprocess_q


postprocess_e: postprocess_e.c
	$(CPP) -g $(CPPFLAGS) -o $@ $^

postprocess_q: postprocess_q.c
	$(CPP) -g $(CPPFLAGS) -o $@ $^

postprocess_qq: postprocess_qq.c
	$(CPP) -g $(CPPFLAGS) -o $@ $^
