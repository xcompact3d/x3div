SUBDIRS = decomp2d src verification .

bin_PROGRAMS = xcompact3d
xcompact3d_SOURCES = src/xcompact3d.f90
xcompact3d_FCFLAGS = -I$(top_builddir)/decomp2d -I$(top_builddir)/src
xcompact3d_LDFLAGS = -L$(top_builddir)/decomp2d -L$(top_builddir)/src
xcompact3d_LDADD = -ldecomp2d -lx3d
if HAVE_FFTW3
xcompact3d_LDADD += @FFTW3LDFLAGS@
endif
if HAVE_MKL
xcompact3d_LDADD += @MKLLDFLAGS@
endif

include_HEADERS = \
$(top_builddir)/decomp2d/decomp_2d.mod \
$(top_builddir)/decomp2d/glassman.mod \
$(top_builddir)/decomp2d/decomp_2d_fft.mod \
$(top_builddir)/src/bc_tgv2d.mod \
$(top_builddir)/src/case.mod \
$(top_builddir)/src/decomp_2d_poisson.mod \
$(top_builddir)/src/navier.mod \
$(top_builddir)/src/param.mod \
$(top_builddir)/src/parameters.mod \
$(top_builddir)/src/thomas.mod \
$(top_builddir)/src/transeq.mod \
$(top_builddir)/src/variables.mod \
$(top_builddir)/src/var.mod \
$(top_builddir)/src/x3d_derive.mod \
$(top_builddir)/src/x3d_filters.mod \
$(top_builddir)/src/x3d_operator_1d.mod \
$(top_builddir)/src/x3d_operator_x_data.mod \
$(top_builddir)/src/x3d_operator_y_data.mod \
$(top_builddir)/src/x3d_operator_z_data.mod \
$(top_builddir)/src/x3d_precision.mod \
$(top_builddir)/src/x3d_staggered.mod \
$(top_builddir)/src/x3d_transpose.mod
