# Decomp2d library
lib_LTLIBRARIES = libdecomp2d.la
libdecomp2d_la_SOURCES = \
decomp_2d.f90 \
glassman.f90 \
fft_generic.f90
libdecomp2d_la_LDFLAGS = -no-undefined

# Allow parallel make
glassman.lo: decomp_2d.lo
fft_generic.lo: glassman.lo

clean-local:
	-rm -f *.mod
