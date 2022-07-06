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

# FFT options
set(FFT_Choice "generic" CACHE STRING "FFT for XCompact3d project (generic is the default)")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw)

# Look for fftw is required
if(${FFT_Choice} MATCHES "fftw")
  configure_file(cmake/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
      else()
      message("CMake step for findFFTW completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
  endif()
  
  set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
  
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
  find_package(FFTW)
  message(STATUS "FFTW_FOUND     : ${FFTW_FOUND}")
  message(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
  message(STATUS "FFTW_INCLUDE   : ${FFTW_INCLUDE_DIRS}")

  #add_definitions("-lfftw3 -lfftw3f")

endif(${FFT_Choice} MATCHES "fftw")

# Create a static library for decomposition and fft
add_subdirectory(decomp2d)

# Create the Xcompact3d executable
add_subdirectory(src)

include(GNUInstallDirs)
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
install(TARGETS xcompact3d
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

