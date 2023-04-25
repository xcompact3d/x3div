cmake_minimum_required(VERSION 3.12)                        
cmake_policy(SET CMP0074 NEW)

project(xcompact3d LANGUAGES Fortran)


set(BUILD_TARGET "mpi" CACHE STRING "Target for acceleration (mpi (default) or gpu)")
set_property(CACHE BUILD_TARGET PROPERTY STRINGS mpi gpu)

if (BUILD_TARGET MATCHES "gpu")
  option(ENABLE_OPENACC "Allow user to activate/deactivate OpenACC support" ON)
  option(ENABLE_CUDA "Allow user to activate/deactivate CUDA support" ON)
  option(ENABLE_MANAGED "Allow user to activate/deactivate automatic memory managment from NVHPC" OFF)
  option(ENABLE_NCCL "Allow user to activate/deactivate Collective Comunication NCCL" OFF)

  if (ENABLE_CUDA)
    message(STATUS "Before enable CUDA")
    enable_language(CUDA)
    message(STATUS "After enable CUDA")
  endif()

endif(BUILD_TARGET MATCHES "gpu")

#set (ENABLE_PROFILER "OFF" CACHE STRING "Activate/select the profiler")
#set_property(CACHE ENABLE_PROFILER PROPERTY STRINGS OFF caliper)

set(AUTHOR "Stefano Rolfo;Charles Moulinec;Paul Bartholomew")
set(AUTHOR_DETAILS "stefano.rolfo@stfc.ac.uk;charles.moulinec@stfc.ac.uk;p.bartholomew@epcc.ed.ac.uk")
set(DESCRIPTION "Building x3div using cmake")

message(STATUS "building ${PROJECT_NAME}")

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/opt" CACHE PATH "..." FORCE)
endif()

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Dev Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/compilers" "${CMAKE_SOURCE_DIR}/cmake/fft")
find_package(DECOMP2D REQUIRED)

if (BUILD_TARGET MATCHES "gpu")
  include(X3D_GPU)
endif (BUILD_TARGET MATCHES "gpu")

include(X3D_MPI)

include(X3D_Compilers)
#include(D2D_Profilers)

## FFT options
#if (BUILD_TARGET MATCHES "gpu")
#  set(FFT_Choice "cufft" CACHE STRING "FFT for XCompact3d project (with GPU cufft is the default)")
#else ()
#  set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
#endif (BUILD_TARGET MATCHES "gpu")
#set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw fftw_f03 mkl cufft)
#include(fft)

# Add the sources
add_subdirectory(src)


# Add a prettify target
#add_custom_target(format sh ${CMAKE_SOURCE_DIR}/scripts/format.sh
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})


