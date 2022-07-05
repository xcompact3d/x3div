cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(xcompact3d LANGUAGES Fortran)

find_package(MPI REQUIRED)
if (MPI_Fortran_COMPILER)
  message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
else (MPI_Fortran_COMPILER)
  message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)

# Create a static library for the fft
add_subdirectory(decomp2d)

# Create the Xcompact3d executable
add_subdirectory(src)

include(GNUInstallDirs)
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
install(TARGETS xcompact3d
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
