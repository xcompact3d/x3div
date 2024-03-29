# GPU CMakeLists
message(STATUS "Check GPU")

if (ENABLE_OPENACC)
  include(FindOpenACC)
  if(OpenACC_Fortran_FOUND)
    message(STATUS "OpenACC for Fotran Compiler Found, version ${OpenACC_Fortran_VERSION_MAJOR}.${OpenACC_Fortran_VERSION_MINOR}")
  else()
    message(ERROR_CRITICAL "No OpenACC support detected")
  endif()
endif()

if (ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  if(${CMAKE_VERSION} VERSION_LESS_EQUAL "3.13.4")
    cuda_select_nvcc_arch_flags(ARCH_FLAGS "Auto") # optional argument for arch to add
    message(STATUS "ARCH_FLAGS = ${ARCH_FLAGS}")
    string(REPLACE "-gencode;" "--generate-code=" ARCH_FLAGS "${ARCH_FLAGS}")
    string(APPEND CMAKE_CUDA_FLAGS "${ARCH_FLAGS}")
    message(STATUS "ARCH_FLAGS WITH CUDA = ${ARCH_FLAGS}")
  else()
    include(FindCUDA/select_compute_arch)
    CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
    string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
    string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
    string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_3}")
    SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
    set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
    message(STATUS "CUDA_ARCHITECTURES ${CUDA_ARCH_LIST}")
    list(GET CUDA_ARCH_LIST 0 CUDA_ARCH_FIRST)
    message(STATUS "CUDA_ARCH_FIRST ${CUDA_ARCH_FIRST}")
  endif()
  message(STATUS "CUDA_LIBRARIES ${CUDA_LIBRARIES}")
endif()

