file(GLOB files_decomp decomp_2d.f90 glassman.f90)
file(GLOB files_fft fft_generic.f90)
set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)

include(cmake/gfortran_flags.cmake)

if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10.3)
  list(APPEND gfortran_flags -fallow-argument-mismatch)
  # See https://github.com/xcompact3d/x3div/pull/20#discussion_r814948614
endif()
target_compile_options(decomp2d PRIVATE ${gfortran_flags})
