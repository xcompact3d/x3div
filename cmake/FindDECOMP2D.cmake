# - Find the 2decomp-fft library

if( NOT D2D_ROOT AND DEFINED ENV{D2D_DIR} )
  set( D2D_ROOT $ENV{D2D_DIR} )
else()
  message(STATUS "2decomp-fft PATH not available we'll try to download and install")
  configure_file(${CMAKE_SOURCE_DIR}/cmake/decomp2d/downloadBuild2decomp.cmake.in decomp2d-build/CMakeLists.txt)
  #message("Second CMAKE_GENERATOR ${CMAKE_GENERATOR}") 
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build )
  if(result)
      message(FATAL_ERROR "CMake step for 2decomp-fft failed: ${result}")
  else()
      message("CMake step for 2decomp-fft completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
         RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build )
  if(result)
      message(FATAL_ERROR "Build step for 2decomp-fft failed: ${result}")
  endif()
  set(D2D_ROOT ${CMAKE_CURRENT_BINARY_DIR}/decomp2d-build/downloadBuild2decomp-prefix/src/downloadBuild2decomp-build/opt)
endif()
message(STATUS "D2D_ROOT ${D2D_ROOT}")

# find libs
find_library(
    DECOMP2D_FFT_LIB
    NAMES "decomp2d" libdecomp2d
    PATHS ${D2D_ROOT}
    PATH_SUFFIXES "lib" "lib64"
    NO_DEFAULT_PATH
    REQUIRED
    )
set(2decomp_INCLUDE_DIR "${D2D_ROOT}/include")
message(STATUS "D2D found at: ${D2D_ROOT}")


