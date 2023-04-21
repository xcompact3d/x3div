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
PARAMOD = mpi # multicore,gpu


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
    #FFLAGS += -fopenmp -ftree-parallelize-loops=12
    #LFLAGS += -fopenmp
  endif
else ifeq ($(CMP),nagfor)
  FC = mpinagfor
  FFLAGS += -fpp
else ifeq ($(CMP),cray)
  FC = ftn
  FFLAGS += -eF -g -O3 -N 1023
  FFLAGS += -h omp -h thread_do_concurrent
  LFLAGS += -h omp -h thread_do_concurrent
else ifeq ($(CMP),nvhpc)
  FC = mpif90
  FFLAGS += -cpp
  ifeq ($(PARAMOD),multicore)
     FFLAGS += -O3 -Minfo=accel -stdpar -acc -target=multicore
     LFLAGS += -acc -lnvhpcwrapnvtx
  else ifeq ($(PARAMOD),gpu)
     CCXY=80
     MANAGED=no
     ifeq ($(MANAGED),yes)
       GPUOPT=-gpu=cc${CCXY},managed,lineinfo
     else
       GPUOPT=-gpu=cc${CCXY},lineinfo
     endif
     NCCL=no
     #FFLAGS += -D_GPU
     #ifeq ($(NCCL),yes)
     #  FFLAGS += -D_NCCL
     #endif
     FFLAGS += -Mfree -Kieee -Minfo=accel,stdpar ${GPUOPT} -acc -target=gpu -traceback -O3 -DUSE_CUDA -cuda
     ifeq ($(NCCL),yes)
       #FFLAGS += -cudalib=cufft,nccl
       LFLAGS += -cudalib=cufft,nccl
     else
       #FFLAGS += -cudalib=cufft
       LFLAGS += -cudalib=cufft
     endif
     #FFLAGS += -D_GPU -Mfree -Kieee -Minfo=accel,stdpar -stdpar=gpu -gpu=cc80,managed,lineinfo -acc -target=gpu -traceback -O3 -DUSE_CUDA -cuda -cudalib=cufft
     #FFLAGS += -Mfree -Kieee -Minfo=accel,stdpar -stdpar=gpu -gpu=cc80,managed,lineinfo -acc -target=gpu -traceback -O3 -DUSE_CUDA -cuda -cudalib=cufft
     LFLAGS += -lnvhpcwrapnvtx
  else
    FFLAGS += -cpp -O3 -march=native
  endif
  #FFLAGS += -cpp -O3 -Minfo=accel -stdpar -acc -target=multicore
  #FFLAGS = -cpp -Mfree -Kieee -Minfo=accel -g -acc -target=gpu -fast -O3 -Minstrument
endif

SRCDIR = ./src

### List of files for the main code
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/x3d_precision.f90 $(SRCDIR)/module_param.f90 $(SRCDIR)/time_integrators.f90 $(SRCDIR)/x3d_transpose.f90 $(SRCDIR)/var.f90 $(SRCDIR)/thomas.f90 $(SRCDIR)/x3d_operator_x_data.f90 $(SRCDIR)/x3d_operator_y_data.f90 $(SRCDIR)/x3d_operator_z_data.f90 $(SRCDIR)/x3d_operator_1d.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/x3d_derive.f90 $(SRCDIR)/x3d_staggered.f90 $(SRCDIR)/x3d_filters.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/mom.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/x3d_tools.f90 $(SRCDIR)/xcompact3d.f90

########Decomp2d settings##########
DECOMP_ROOT = /Users/rfj82982/GIT_projects/2decomp-fft
DECOMP_BUILD_DIR = $(DECOMP_ROOT)/build
DECOMP_INSTALL_DIR ?= $(DECOMP_ROOT)/build/opt
# Use default unless set by user

INC += -I$(DECOMP_INSTALL_DIR)/include

# Users build/link targets
LIBS = -L$(DECOMP_INSTALL_DIR)/lib64 -L$(DECOMP_INSTALL_DIR)/lib -ldecomp2d

#######FFT settings##########
ifeq ($(FFT),fftw3)
  FFTW3_PATH=/opt/local/lib
  LIBFFT=-L$(FFTW3_PATH) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH=/usr                                #ubuntu # apt install libfftw3-dev
  LIBFFT=-L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  LIBFFT=
else ifeq ($(FFT),mkl)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  #INC+=-I$(MKLROOT)/include
else ifeq ($(FFT),cufft)
  #CUFFT_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.1/math_libs                                
  #INC+=-I${NVHPC}/Linux_x86_64/${EBVERSIONNVHPC}/compilers/include
endif

#######OPTIONS settings###########
#OPT = -I$(SRCDIR) -I$(DECOMPDIR) -I$(DECOMP_INSTALL_DIR)/include 
OPT = -I$(SRCDIR) -I$(DECOMP_INSTALL_DIR)/include 
#LINKOPT = $(FFLAGS) -lnvhpcwrapnvtx
LINKOPT = $(FFLAGS)
LIBS += $(LIBFFT) $(LFLAGS) 
#-----------------------------------------------------------------------
# Normally no need to change anything below

all: xcompact3d

xcompact3d : $(OBJDECOMP) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OPT) $(OBJ) $(LIBS)

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(DEFS) $(INC) -c $<
	mv $(@F) ${SRCDIR}

.PHONY: post
post:
	$(FC) $(FFLAGS) $(DEFS) $(DEFS2) post.f90 -c
	$(FC) $(FFLAGS) -o $@ $(PSRC:.f90=.o)

.PHONY: clean


clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod xcompact3d post

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
