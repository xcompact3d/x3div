!---------------------------------------------------------------------------
! This file is external to the pfunit library and is used to represent
! information that cannot be determined at run-time.   In particular,
! Fortran lacks reflection/introspection and therefore cannot automatically
! collect user-defined test suites.
!
!
! Typical usage is for the user to compile this file while leaving it in the
! pFUnit installation directory.    Use "-D<_TEST_SUITES>" to specify a file
! that contains the list of test_suites to be builtlinked and executed.
!
!
! For serial runs, the user links with the FUnit library, while for parallel
! runs the user links with FUnit _and_ pFUnit.
!---------------------------------------------------------------------------

module loader
   use FUnit, only: TestSuite
   implicit none

contains

   function load_tests() result(suite)

      type (TestSuite) :: suite

#define ADDULE_TEST_SUITE(m,s) ! do nothing
#define ADD_TEST_SUITE(s) type (TestSuite), external :: s
#ifdef x3dinit
ADD_TEST_SUITE(x3dinit)
#endif
#ifdef x3dtest
ADD_TEST_SUITE(x3dtest)
#endif
#ifdef x3dfin
ADD_TEST_SUITE(x3dfin)
#endif
#undef ADD_TEST_SUITE
#undef ADDULE_TEST_SUITE

      suite = TestSuite()

#define ADD_TEST_SUITE(s) call suite%addTest(s())
#define ADDULE_TEST_SUITE(m,s) call suite%addTest(s())
#ifdef x3dinit
ADD_TEST_SUITE(x3dinit)
#endif
#ifdef x3dtest
ADD_TEST_SUITE(x3dtest)
#endif
#ifdef x3dfin
ADD_TEST_SUITE(x3dfin)
#endif
#undef ADD_TEST_SUITE
#undef ADDULE_TEST_SUITE

   end function load_tests

end module loader

program main
   use FUnit, only : stub
   use loader
#ifdef PFUNIT_EXTRA_USE
      ! Use external code for whatever suite-wide fixture is in use.
      use PFUNIT_EXTRA_USE
#endif
   implicit none

   procedure(), pointer :: extra_initialize
   procedure(), pointer :: extra_finalize

#ifdef PFUNIT_EXTRA_INITIALIZE
#  ifndef PFUNIT_EXTRA_USE
   external :: PFUNIT_EXTRA_INITIALIZE
#  endif
   extra_initialize => PFUNIT_EXTRA_INITIALIZE
#else
   extra_initialize => stub
#endif

#ifdef PFUNIT_EXTRA_FINALIZE
#  ifndef PFUNIT_EXTRA_USE
   external :: PFUNIT_EXTRA_FINALIZE
#  endif
   extra_finalize => PFUNIT_EXTRA_FINALIZE
#else
   extra_finalize => stub
#endif

   call funit_main(load_tests, extra_initialize, extra_finalize)

end program main
