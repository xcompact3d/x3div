#Compilers Flags for NVIDIA

set(X3D_FFLAGS "-cpp -Mfree -Kieee")
set(X3D_FFLAGS_RELEASE "-O3")
set(X3D_FFLAGS_DEBUG   "-O0 -g -traceback -Mbounds -Mchkptr -Ktrap=fp")
set(X3D_FFLAGS_DEV     "${X3D_FFLAGS_DEBUG}")
if (BUILD_TARGET MATCHES "mpi")
  set(X3D_FFLAGS_RELEASE "${X3D_FFLAGS_RELEASE} -fast -march=native")
endif (BUILD_TARGET MATCHES "mpi")

if (BUILD_TARGET MATCHES "gpu")
  set(X3D_FFLAGS "${X3D_FFLAGS} -Minfo=accel -target=gpu")
  add_definitions("-D_GPU")
  if (ENABLE_OPENACC)
    set(X3D_FFLAGS "${X3D_FFLAGS} -acc")
  endif()
  if (ENABLE_CUDA)
    add_definitions("-DUSE_CUDA")
    set(X3D_FFLAGS "${X3D_FFLAGS} -cuda")
    # Add Compute Capabilities and memory managemnt 
    if (ENABLE_MANAGED)
      set(X3D_FFLAGS "${X3D_FFLAGS} -gpu=cc${CUDA_ARCH_FIRST},managed,lineinfo")
    else (ENABLE_MANAGED)
      set(X3D_FFLAGS "${X3D_FFLAGS} -gpu=cc${CUDA_ARCH_FIRST},lineinfo")
    endif(ENABLE_MANAGED)
    # Add NCCL cuFFT
    if (ENABLE_NCCL)
      add_definitions("-D_NCCL")
      set(X3D_FFLAGS "${X3D_FFLAGS} -cudalib=nccl,cufft")
    else(ENABLE_NCCL)
      set(X3D_FFLAGS "${X3D_FFLAGS} -cudalib=cufft")
    endif(ENABLE_NCCL)
  endif(ENABLE_CUDA)
  # Add profiler
  #if (ENABLE_PROFILER)
  #  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lnvhpcwrapnvtx")
  #endif(ENABLE_PROFILER)
endif (BUILD_TARGET MATCHES "gpu")
