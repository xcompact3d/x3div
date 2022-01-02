include_directories(${CMAKE_SOURCE_DIR}/decomp2d
	            ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d 
                    x3d_tools.f90
                    x3d_precision.f90
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
                    xcompact3d.f90)
	    #tools.f90
#add_executable(xcompact3d ${files_xcompact})
target_include_directories(xcompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_link_libraries(xcompact3d PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS xcompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
