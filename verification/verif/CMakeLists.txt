include_directories(${PROJECT_BINARY_DIR}/decomp2d
                    ${PROJECT_BINARY_DIR}/src
                    ${PROJECT_BINARY_DIR}/verification/lib)

add_pfunit_ctest (verif_x_00
  TEST_SOURCES initx00.pf testx00.pf final.pf
  LINK_LIBRARIES decomp2d x3d x3dverif MPI::MPI_Fortran PFUNIT::pfunit)


install(TARGETS verif_x_00
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
