!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!This file is adapted from https://github.com/Goddard-Fortran-Ecosystem/pFUnit/

!Copyright 2005 United States Government as represented by the
!Administrator of The National Aeronautics and Space
!Administration. All Rights Reserved.

module loader
   use FUnit, only: TestSuite
   implicit none

contains

   function load_tests() result(suite)

#define ADDULE_TEST_SUITE(m,s) use m, only: s
#define ADD_TEST_SUITE(s) ! do nothing
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
