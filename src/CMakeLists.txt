set(SRCFILES
  xcompact3d.f90
  bc_tgv2d.f90
  bc_tgv.f90
  case.f90
  module_param.f90
  navier.f90
  parameters.f90
  poisson.f90
  thomas.f90
  time_integrators.f90
  tools.f90
  transeq.f90
  var.f90
  x3d_derive.f90
  x3d_filters.f90
  x3d_operator_1d.f90
  x3d_operator_x_data.f90
  x3d_operator_y_data.f90
  x3d_operator_z_data.f90
  x3d_staggered.f90
  x3d_tools.f90
  x3d_transpose.f90
  x3d_precision.f90
  )

add_executable(xcompact3d ${SRCFILES})

target_include_directories(xcompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_link_libraries(xcompact3d PRIVATE decomp2d)
target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)

include(cmake/gfortran_flags.cmake)
if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10.3)
  list(APPEND gfortran_flags -fallow-argument-mismatch)
  # See https://github.com/xcompact3d/x3div/pull/20#discussion_r814948614
endif()
target_compile_options(xcompact3d PRIVATE ${gfortran_flags})
