include_directories(${2decomp_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src)
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

set(X3DFILES 
            case.f90
            module_param.f90
            mom.f90
            navier.f90
            parameters.f90
            poisson.f90
            thomas.f90
            time_integrators.f90
	    tools.f90
            transeq.f90
            var.f90
	    x3d_derive.f90
	    x3d_filters.f90
	    x3d_operator_1d.f90
	    x3d_operator_x_data.f90
	    x3d_operator_y_data.f90
	    x3d_operator_z_data.f90
	    x3d_staggered.f90
            x3d_tools.f90
	    x3d_transpose.f90
            x3d_precision.f90)

add_library(x3d STATIC ${X3DFILES})
add_executable(xcompact3d xcompact3d.f90)

target_link_libraries(x3d PRIVATE decomp2d)
target_link_libraries(xcompact3d PRIVATE decomp2d x3d)

if (MPI_FOUND)
  target_link_libraries(x3d PRIVATE MPI::MPI_Fortran)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS xcompact3d x3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


