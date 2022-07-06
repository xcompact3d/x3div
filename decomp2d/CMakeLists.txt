file(GLOB files_decomp decomp_2d.f90 glassman.f90)


if(FFTW_FOUND)
  message (STATUS "Compiling using FFTW3")
  file(GLOB files_fft fft_fftw3.f90)
else(FFTW_FOUND)
  message (STATUS "Compiling using Generic FFT")
  file(GLOB files_fft fft_generic.f90)
endif(FFTW_FOUND)


set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
if(FFTW_FOUND)
  target_include_directories(decomp2d PRIVATE ${FFTW_INCLUDE_DIRS})
endif(FFTW_FOUND)

include(cmake/gfortran_flags.cmake)

if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10.3)
  list(APPEND gfortran_flags -fallow-argument-mismatch)
  # See https://github.com/xcompact3d/x3div/pull/20#discussion_r814948614
endif()
target_compile_options(decomp2d PRIVATE ${gfortran_flags})
