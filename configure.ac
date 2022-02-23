# Name / version / website
AC_INIT([x3div], [alpha], [https://github.com/xcompact3d/x3div/])
AC_LANG(Fortran)
AC_CONFIG_MACRO_DIR([./m4])

AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])

# Save the flags provided by the user
user_FCFLAGS=$FCFLAGS
#user_LDFLAGS=$LDFLAGS

# Check for a MPI Fortran compiler if FC is not provided
AX_PROG_FC_MPI([],[],[])
AC_MSG_NOTICE([Compiler detected : $FC])

# Single precision build, off by default
AC_ARG_ENABLE([single], AS_HELP_STRING([--enable-single], [single precision build (default is double)]))

# Debug mode, off by default
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [enable debugging (default is no)]))

# FFT, fftw3
AC_ARG_WITH(fftw3,
            [AS_HELP_STRING([--with-fftw3=PATH],[prefix directory for fftw3])],
            [if test "x$withval" = "x"; then
		with_fftw3=check
	    fi],
	    [with_fftw3=check])
AC_ARG_WITH(fftw3-include,
	    [AS_HELP_STRING([--with-fftw3-include=PATH], [directory for fftw3 include files])],
	    [if test "x$with_fftw3" = "xcheck"; then
		with_fftw3=yes
	     fi
	     FFTW3_FCFLAGS="-I$with_fftw3_include"],
	    [if test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xyes" \
		  -a "x$with_fftw3" != "xcheck"; then
		FFTW3_FCFLAGS="-I$with_fftw3/include"
	     fi])
AC_ARG_WITH(fftw3-lib,
	    [AS_HELP_STRING([--with-fftw3-lib=PATH], [directory for fftw3 library])],
	    [if test "x$with_fftw3" = "xcheck"; then
		with_fftw3=yes
	     fi
	     FFTW3_LDFLAGS="-L$with_fftw3_lib"],
	    [if test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xyes" \
		  -a "x$with_fftw3" != "xcheck"; then
	     FFTW3_LDFLAGS="-L$with_fftw3/lib"
	     fi])
# FFT, mkl
AC_ARG_WITH(mkl,
            [AS_HELP_STRING([--with-mkl=PATH],[prefix directory for mkl])],
            [if test "x$withval" = "x"; then
		with_mkl=check
	    fi],
	    [with_mkl=check])
AC_ARG_WITH(mkl-include,
	    [AS_HELP_STRING([--with-mkl-include=PATH], [directory for mkl include files])],
	    [if test "x$with_mkl" = "xcheck"; then
		with_mkl=yes
	     fi
	     MKL_FCFLAGS="-I$with_mkl_include"],
	    [if test "x$with_mkl" != "xno" -a "x$with_mkl" != "xyes" \
		  -a "x$with_mkl" != "xcheck"; then
		MKL_FCFLAGS="-I$with_mkl/include"
	     fi])
AC_ARG_WITH(mkl-lib,
	    [AS_HELP_STRING([--with-mkl-lib=PATH], [directory for mkl library])],
	    [if test "x$with_mkl" = "xcheck"; then
		with_mkl=yes
	     fi
	     MKL_LDFLAGS="-L$with_mkl_lib"],
	    [if test "x$with_mkl" != "xno" -a "x$with_mkl" != "xyes" \
		  -a "x$with_mkl" != "xcheck"; then
	     MKL_LDFLAGS="-L$with_mkl/lib"
	     fi])
# Default FFT engine is generic
AC_SUBST(FFTENGINE, generic)
if test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xcheck"; then
	AC_MSG_NOTICE([FFT engine fftw3 was provided])
	AS_IF([test "x$enable_single" = "xyes"],
	      [FFTW3_LDFLAGS="$FFTW3_LDFLAGS -lfftw3f"],
              [FFTW3_LDFLAGS="$FFTW3_LDFLAGS -lfftw3"])
	AC_SUBST(FFTENGINE, fftw3)
	AC_SUBST(FFTW3FCFLAGS, ${FFTW3_FCFLAGS})
	AC_SUBST(FFTW3LDFLAGS, ${FFTW3_LDFLAGS})
fi
AM_CONDITIONAL(HAVE_FFTW3, [test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xcheck"])
if test "x$with_mkl" != "xno" -a "x$with_mkl" != "xcheck"; then
	AC_MSG_NOTICE([FFT engine mkl was provided])
	AC_SUBST(FFTENGINE, mkl)
	AC_SUBST(MKLFCFLAGS, ${MKL_FCFLAGS})
        AC_SUBST(MKLLDFLAGS, ${MKL_LDFLAGS})
fi
AM_CONDITIONAL(HAVE_MKL, [test "x$with_mkl" != "xno" -a "x$with_mkl" != "xcheck"])

#
# Check the Fortran compiler
#
Fcompiler=none
#
# 1. gfortran ?
#
echo `$FC --version 2>&1 | head -1 | grep 'GNU Fortran'` > /dev/null
if test "$?" = "0" ; then
	Fcompiler=ok
	fc_version="`$FC -dumpversion | cut -f1 -d.`"
	AC_MSG_NOTICE([GNU Fortran version $fc_version])
	if test "x$enable_single" = "xyes"; then
		FCFLAGS="-cpp"
	else
                FCFLAGS="-cpp -DDOUBLE_PREC"
        fi
	AS_IF([test "x$enable_debug" = "xyes"],
	      [FCFLAGS="$FCFLAGS -Og -g3 -ffpe-trap=invalid,zero -fcheck=bounds -fimplicit-none -DDEBUG"],
	      [FCFLAGS="$FCFLAGS -O3 -march=native"])
fi
#
# 2. ifort ?
#
if test "x$Fcompiler" != "xok"; then
echo `$FC --version 2>&1 | head -1 | grep 'ifort'` > /dev/null
if test "" = ""; then
	Fcompiler=ok
	fc_version=""
fi
fi
#
# 3. nvfortran ?
#

#
# 4+ TODO
#

# Combine default flags with the user ones
FCFLAGS="$FCFLAGS $user_FCFLAGS"

# Maintener mode if a .git folder is found
AS_IF([test -d ${srcdir}/.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-static])

AC_CONFIG_FILES([Makefile src/Makefile decomp2d/Makefile])

AC_OUTPUT
