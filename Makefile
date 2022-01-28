#=======================================================================
# Makefile for Xcompact3D
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
#   -DDEBG        - debuggin xcompact3d.f90
# generate a Git version string
GIT_VERSION := $(shell git describe --tag --long --always)

DEFS = -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

LCL = local# local,lad,sdu,archer
IVER = 17# 15,16,17,18
CMP = gcc# intel,gcc
FFT = generic# generic,fftw3,mkl

#######CMP settings###########
ifeq ($(CMP),intel)
FC = mpiifort
#FFLAGS = -fpp -O3 -xHost -heap-arrays -shared-intel -mcmodel=large -safe-cray-ptr -g -traceback
FFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr -I$(MPI_ROOT)/lib
##debuggin test: -check all -check bounds -chintel eck uninit -gen-interfaces -warn interfaces
else ifeq ($(CMP),gcc)
FC = mpif90
#FFLAGS = -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -x f95-cpp-input
FFLAGS = -cpp -O3 -funroll-loops -floop-optimize -g -Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none
#-ffpe-trap=invalid,zero
else ifeq ($(CMP),nagfor)
FC = mpinagfor
FFLAGS = -fpp
else ifeq ($(CMP),cray)
FC = ftn
FFLAGS = -eF -g -O3 -N 1023
endif


MODDIR = ./mod
DECOMPDIR = ./decomp2d
SRCDIR = ./src

### List of files for the main code
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/parameters.f90 #$(SRCDIR)/*.f90
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/module_param.f90 $(SRCDIR)/variables.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/derive.f90 $(SRCDIR)/schemes.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/mom.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/xcompact3d.f90

#######OPTIONS settings###########
OPT = -I$(SRCDIR) -I$(DECOMP_INCDIR) $(FFLAGS)
LINKOPT = $(FFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

DECOMP_LIB = 2decomp_fft
DECOMP_LIBDIR = $(DECOMPDIR)/lib
DECOMP_INCDIR = $(DECOMPDIR)/include
DECOMP.A = $(DECOMP_LIBDIR)/lib$(DECOMP_LIB).a
include $(DECOMPDIR)/src/Makefile.inc

all: xcompact3d

xcompact3d : $(DECOMP.A) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJ) -L$(DECOMP_LIBDIR) -l$(DECOMP_LIB) $(LIBFFT)

$(DECOMP.A):
	make -C $(DECOMPDIR) F90=$(FC) OPTIONS="$(FFLAGS) $(DEFS) $(DEFS2)" lib

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${SRCDIR}
	#mv *.mod ${SRCDIR}

## This %.o : %.f90 doesn't appear to be called...
%.o : %.f90
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) $(INC) -c $<

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean


clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod xcompact3d post

.PHONY: clean-decomp
clean-decomp:
	make -C $(DECOMPDIR) clean

.PHONY: cleanall
cleanall: clean clean-decomp
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
