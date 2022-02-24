cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(xcompact3d LANGUAGES Fortran)

find_package(MPI REQUIRED)
if (MPI_Fortran_COMPILER)
  message(STATUS "MPI_Fortran_COMPILER found: ${MPI_Fortran_COMPILER}")
else (MPI_Fortran_COMPILER)
  message(SEND_ERROR "This application cannot compile without MPI")
endif(MPI_Fortran_COMPILER)

# Warning if Include are not found => can be fixed with more recent cmake version
if (MPI_FOUND)
  message(STATUS "MPI FOUND: ${MPI_FOUND}")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  message(STATUS "MPI INCL ALSO FOUND: ${MPI_INCLUDE_PATH}")
else (MPI_FOUND)
  message(STATUS "NO MPI include have been found. The executable won't be targeted with MPI include")
  message(STATUS "Code will compile but performaces can be compromised")
  message(STATUS "Using a CMake vers > 3.10 should solve the problem")
  message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (MPI_FOUND)


set(CMAKE_Fortran_FLAGS "-cpp -funroll-loops -floop-optimize -g
-Warray-bounds -fcray-pointer -fbacktrace -ffree-line-length-none")

# Create a static library for the fft
add_subdirectory(decomp2d)

# Create the Xcompact3d executable
add_subdirectory(src)
