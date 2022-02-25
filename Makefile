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
CMP = nvhpc# intel,gcc,nagfor,cray,nvhpc
FFT = generic# fftw3,fftw3_f03,generic,mkl

BUILD ?= # debug can be used with gcc
FCFLAGS ?= # user can set default compiler flags
LDFLAGS ?= # user can set default linker flags
FFLAGS = $(FCFLAGS)
LFLAGS = $(LDFLAGS)

#######CMP settings###########
ifeq ($(CMP),intel)
  FC = mpiifort
  FFLAGS += -fpp -O3 -mavx2 -march=core-avx2 -mtune=core-avx2
  FFLAGS += -fopenmp
  LFLAGS += -fopenmp
else ifeq ($(CMP),gcc)
  FC = mpif90
  FFLAGS += -cpp
  ifeq "$(shell expr `gfortran -dumpversion | cut -f1 -d.` \>= 10)" "1"
    FFLAGS += -fallow-argument-mismatch
  endif
  ifeq ($(BUILD),debug)
    DEFS += -DDEBUG
    FFLAGS += -g3 -Og
    FFLAGS += -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none
  else
    FFLAGS += -O3 -march=native
    FFLAGS += -fopenmp -ftree-parallelize-loops=12
    LFLAGS += -fopenmp
  endif
else ifeq ($(CMP),nagfor)
  FC = mpinagfor
  FFLAGS += -fpp
else ifeq ($(CMP),cray)
  FC = ftn
  FFLAGS += -eF -g -O3 -N 1023
else ifeq ($(CMP),nvhpc)
  FC = mpif90
  FFLAGS += -cpp -O3 -march=native
  FFLAGS += -Minfo=accel -stdpar -acc -target=multicore
#  FFLAGS = -cpp -Mfree -Kieee -Minfo=accel -g -acc -target=gpu -fast -O3 -Minstrument
  LFLAGS += -acc -lnvhpcwrapnvtx
endif

DECOMPDIR = ./decomp2d
SRCDIR = ./src

### List of files for the main code
SRCDECOMP = $(DECOMPDIR)/decomp_2d.f90 $(DECOMPDIR)/glassman.f90 $(DECOMPDIR)/fft_$(FFT).f90 
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/x3d_precision.f90 $(SRCDIR)/module_param.f90 $(SRCDIR)/tools.f90 $(SRCDIR)/x3d_transpose.f90 $(SRCDIR)/var.f90 $(SRCDIR)/thomas.f90 $(SRCDIR)/x3d_operator_x_data.f90 $(SRCDIR)/x3d_operator_y_data.f90 $(SRCDIR)/x3d_operator_z_data.f90 $(SRCDIR)/x3d_operator_1d.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/x3d_derive.f90 $(SRCDIR)/x3d_staggered.f90 $(SRCDIR)/x3d_filters.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/bc_tgv2d.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/x3d_tools.f90 $(SRCDIR)/xcompact3d.f90

#######FFT settings##########
ifeq ($(FFT),fftw3)
  #FFTW3_PATH=/usr
  #FFTW3_PATH=/usr/lib64
  FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH=/usr                                #ubuntu # apt install libfftw3-dev
  #FFTW3_PATH=/usr/lib64                         #fedora # dnf install fftw fftw-devel
  #FFTW3_PATH=/usr/local/Cellar/fftw/3.3.7_1     #macOS  # brew install fftw
  INC=-I$(FFTW3_PATH)/include
  LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  INC=
  LIBFFT=
else ifeq ($(FFT),mkl)
  SRCDECOMP := $(DECOMPDIR)/mkl_dfti.f90 $(SRCDECOMP)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC=-I$(MKLROOT)/include
endif

#######OPTIONS settings###########
OPT = -I$(SRCDIR) -I$(DECOMPDIR)
LINKOPT = $(LFLAGS)
#-----------------------------------------------------------------------
# Normally no need to change anything below

all: xcompact3d

xcompact3d : $(OBJDECOMP) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(LIBFFT)

$(OBJDECOMP):$(DECOMPDIR)%.o : $(DECOMPDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $<
	mv $(@F) ${DECOMPDIR}

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $<
	mv $(@F) ${SRCDIR}

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean


clean:
	rm -f $(DECOMPDIR)/*.o $(DECOMPDIR)/*.mod
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod xcompact3d post

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
