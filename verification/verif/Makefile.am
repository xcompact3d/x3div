if HAVE_PFUNIT

USEMPI=YES
include @PFUNITDIR@/include/PFUNIT.mk

AM_FCFLAGS = -I$(top_builddir)/decomp2d -I$(top_builddir)/src -I$(top_builddir)/verification/lib $(PFUNIT_EXTRA_FFLAGS)
AM_LDFLAGS = -L$(top_builddir)/decomp2d -L$(top_builddir)/src -L$(top_builddir)/verification/lib $(PFUNIT_LIBRARIES) -fopenmp

pffiles=\
initx00.pf \
final.pf \
testx00.pf

f90files = $(pffiles:%.pf=%.F90)

$(f90files): $(blddir)%.F90 : $(blddir)%.pf
	@PFUNITDIR@/bin/funitproc $< $@

bin_PROGRAMS = test_x_00
test_x_00_SOURCES = x3d_driver.F90 initx00.F90 testx00.F90 final.F90
test_x_00_FCFLAGS = ${AM_FCFLAGS} -Dx3dinit=initx00_suite -Dx3dtest=testx00_suite -Dx3dfin=final_suite
test_x_00_LDADD = -ldecomp2d -lx3d -lx3dtest
if HAVE_FFTW3
test_x_00_LDADD += @FFTW3LDFLAGS@
endif
if HAVE_MKL
test_x_00_LDADD += @MKLLDFLAGS@
endif
test_x_00_DEPENDENCIES = $(f90files)

endif
