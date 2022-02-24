set(SRCFILES
  xcompact3d.f90
  case.f90
  derive.f90
  module_param.f90
  mom.f90
  navier.f90
  parameters.f90
  poisson.f90
  schemes.f90
  transeq.f90
  variables.f90
  )

add_executable(xcompact3d ${SRCFILES})

target_include_directories(xcompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_link_libraries(xcompact3d PRIVATE decomp2d)
target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)

install(TARGETS xcompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
