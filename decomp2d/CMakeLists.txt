file(GLOB files_decomp decomp_2d.f90 glassman.f90)
file(GLOB files_fft fft_generic.f90)
set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)

include(cmake/gfortran_flags.cmake)
target_compile_options(decomp2d PRIVATE ${gfortran_flags})

install(TARGETS decomp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
