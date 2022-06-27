#
GIT_VERSION := $(shell git describe --tag --long --always)

DEFS = -DDOUBLE_PREC -DVERSION=\"$(GIT_VERSION)\"

IVER = 17# 15,16,17,18
CMP = nvhpc# intel,gcc,nvhpc
FFT = generic# generic,fftw3,mkl,cufft

#######CMP settings###########
ifeq ($(CMP),intel)
FC = mpiifort
FFLAGS = -fpp -O3 -mavx2 -march=core-avx2 -mtune=core-avx2
FFLAGS += -fopenmp
else ifeq ($(CMP),gcc)
FC = mpif90
FFLAGS = -cpp -O3 -march=native
NTHREADS = 1 # Specify NTHREADS on commandline to override
FFLAGS += -fopenmp -ftree-parallelize-loops=$(NTHREADS)
else ifeq ($(CMP),nagfor)
FC = mpinagfor
FFLAGS = -fpp
else ifeq ($(CMP),cray)
FC = ftn
FFLAGS = -eF -g -O3 -N 1023
else ifeq ($(CMP),nvhpc)
FC = mpif90
#FFLAGS += -Minfo=accel -stdpar -acc -target=multicore
FFLAGS = -cpp -D_GPU -D_NCCL -Mfree -Kieee -Minfo=accel,stdpar -stdpar=gpu -gpu=cc80,managed,lineinfo -acc -target=gpu -traceback -O3 -DUSE_CUDA -cuda -cudalib=cufft,nccl
#FFLAGS = -cpp -D_GPU -Mfree -Kieee -Minfo=accel,ftn,inline,loop,vect,opt,stdpar -stdpar=gpu -gpu=cc80,managed,lineinfo -acc -target=gpu -traceback -O3 -DUSE_CUDA -cuda -cudalib=cufft  
endif


MODDIR = ./mod
DECOMPDIR = ./decomp2d
SRCDIR = ./src

### List of files for the main code
SRCDECOMP = $(DECOMPDIR)/decomp_2d.f90 $(DECOMPDIR)/glassman.f90 $(DECOMPDIR)/fft_$(FFT).f90 
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
OBJ = $(SRC:%.f90=%.o)
SRC = $(SRCDIR)/x3d_precision.f90 $(SRCDIR)/module_param.f90 $(SRCDIR)/time_integrators.f90 $(SRCDIR)/tools.f90 $(SRCDIR)/x3d_transpose.f90 $(SRCDIR)/var.f90 $(SRCDIR)/thomas.f90 $(SRCDIR)/x3d_operator_x_data.f90 $(SRCDIR)/x3d_operator_y_data.f90 $(SRCDIR)/x3d_operator_z_data.f90 $(SRCDIR)/x3d_operator_1d.f90 $(SRCDIR)/poisson.f90 $(SRCDIR)/x3d_derive.f90 $(SRCDIR)/x3d_staggered.f90 $(SRCDIR)/x3d_filters.f90 $(SRCDIR)/bc_tgv.f90 $(SRCDIR)/navier.f90 $(SRCDIR)/parameters.f90 $(SRCDIR)/bc_tgv2d.f90 $(SRCDIR)/case.f90 $(SRCDIR)/transeq.f90 $(SRCDIR)/x3d_tools.f90 $(SRCDIR)/xcompact3d.f90

#######FFT settings##########
ifeq ($(FFT),fftw3)
  #FFTW3_PATH=/usr
  #FFTW3_PATH=/usr/lib64
  FFTW3_PATH=/scratch21/eb/gpu/software/FFTW/3.3.10-gompi-2021b
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
  LIBFFT=#-lnvhpcwrapnvtx
else ifeq ($(FFT),mkl)
  SRCDECOMP := $(DECOMPDIR)/mkl_dfti.f90 $(SRCDECOMP)
  LIBFFT=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC=-I$(MKLROOT)/include
else ifeq ($(FFT),cufft)
  #CUFFT_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.1/math_libs                                
  INC=-I${NVHPC}/Linux_x86_64/${EBVERSIONNVHPC}/compilers/include
  #LIBFFT=-L$(CUFFT_PATH)/lib64 -Mcudalib=cufft 
endif

INC=-I${NVHPC}/Linux_x86_64/${EBVERSIONNVHPC}/

#######OPTIONS settings###########
OPT = -I$(SRCDIR) -I$(DECOMPDIR) 
LINKOPT = $(FFLAGS) -lnvhpcwrapnvtx
#-----------------------------------------------------------------------
# Normally no need to change anything below

all: xcompact3d

xcompact3d : $(OBJDECOMP) $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(LIBFFT)

$(OBJDECOMP):$(DECOMPDIR)%.o : $(DECOMPDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(DEFS2) $(INC) -c $<
	mv $(@F) ${DECOMPDIR}
	#mv *.mod ${DECOMPDIR}


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
	rm -f $(DECOMPDIR)/*.o $(DECOMPDIR)/*.mod
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod xcompact3d post

.PHONY: cleanall
cleanall: clean
	rm -f *~ \#*\# out/* data/* stats/* planes/* *.xdmf *.log *.out nodefile core sauve*
