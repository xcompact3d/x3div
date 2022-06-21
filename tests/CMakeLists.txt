find_package(PFUNIT REQUIRED)
enable_testing()

set(test_src test_deriv.pf)
set(other_src
  ../src/x3d_precision.f90
  )
add_pfunit_ctest(tests
  TEST_SOURCES ${test_src}
  OTHER_SOURCES ${other_src}
  LINK_LIBRARIES decomp2d # module param uses "mytype" from decomp2d
  )
target_include_directories(tests PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_compile_options(tests PRIVATE -cpp)

