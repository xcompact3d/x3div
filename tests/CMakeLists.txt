find_package(PFUNIT REQUIRED)
enable_testing()

set(test_src test_dummy.pf)
add_pfunit_ctest(tests TEST_SOURCES ${test_src})

