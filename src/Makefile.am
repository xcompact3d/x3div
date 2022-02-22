AM_FCFLAGS = -I$(top_builddir)/decomp2d

lib_LTLIBRARIES = libx3d.la
libx3d_la_SOURCES = \
x3d_precision.f90 \
module_param.f90 \
x3d_transpose.f90 \
var.f90 \
thomas.f90 \
x3d_operator_x_data.f90 \
x3d_operator_y_data.f90 \
x3d_operator_z_data.f90 \
x3d_operator_1d.f90 \
poisson.f90 \
x3d_derive.f90 \
x3d_staggered.f90 \
x3d_filters.f90 \
navier.f90 \
parameters.f90 \
bc_tgv2d.f90 \
case.f90 \
transeq.f90 \
x3d_tools.f90
libx3d_la_LDFLAGS = -no-undefined

# Allow parallel make
x3d_transpose.lo: module_param.lo
var.lo: module_param.lo
thomas.lo: module_param.lo
x3d_operator_1d.f90: module_param.lo x3d_operator_x_data.lo x3d_operator_y_data.lo x3d_operator_z_data.lo
poisson.lo: x3d_operator_x_data.lo x3d_operator_y_data.lo x3d_operator_z_data.lo
x3d_derive.lo: x3d_operator_1d.lo
x3d_staggered.lo: x3d_operator_1d.lo
x3d_filters.lo: x3d_operator_1d.lo
navier.lo: poisson.lo x3d_staggered.lo
parameters.lo: module_param.lo
bc_tgv2d.lo: x3d_operator_x_data.lo x3d_operator_y_data.lo
case.lo: bc_tgv2d.lo
transeq.lo: x3d_derive.lo
x3d_tools.lo: x3d_derive.lo

clean-local:
	-rm -f *.mod
