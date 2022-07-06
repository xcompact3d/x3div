cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(xcompact3d LANGUAGES Fortran)

find_package(MPI REQUIRED)
if (MPI_Fortran_COMPILER)
  message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
else (MPI_Fortran_COMPILER)
  message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)

# Add defualt definitions such as Version and Double precision
execute_process(
  COMMAND git describe --tag --long --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DVERSION=\"${GIT_VERSION}\"")
option(DOUBLE_PRECISION "Build Xcompact with double precision" ON)
if (DOUBLE_PRECISION)
  add_definitions("-DDOUBLE_PREC")
endif()
# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE MATCHES DEBUG)
  message(STATUS "Build DEBUG ${CMAKE_BUILD_TYPE}")
  add_definitions("-DDEBUG")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

# Create a static library for the fft
add_subdirectory(decomp2d)

# Create the Xcompact3d executable
add_subdirectory(src)

include(GNUInstallDirs)
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
install(TARGETS xcompact3d
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

