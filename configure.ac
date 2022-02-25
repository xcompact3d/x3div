# Name / version / website
AC_INIT([x3div], [alpha], [https://github.com/xcompact3d/x3div/])
AC_LANG(Fortran)
AC_CONFIG_MACRO_DIR([./m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])

#
# Make sure flags are initialized.
#
AM_FCFLAGS="${AM_FCFLAGS}"
AM_LDFLAGS="${AM_LDFLAGS}"
FCFLAGS="${FCFLAGS}"
LDFLAGS="${LDFLAGS}"

#
# Save the flags provided by the user
#
user_FCFLAGS="$FCFLAGS"
#user_LDFLAGS="$LDFLAGS"

#
# Check for a MPI Fortran compiler if FC is not provided
#
AX_PROG_FC_MPI([],[],[])
AC_MSG_NOTICE([Compiler detected : $FC])

#
# Single precision build, off by default
#
AC_ARG_ENABLE([single], AS_HELP_STRING([--enable-single], [single precision build [default is double]]))

#
# Single precision snapshots, off by default
#
AC_ARG_ENABLE([single-snapshot], AS_HELP_STRING([--enable-single-snapshot], [single precision snapshot for double precision build [default is double]]))

#
# Debug mode, off by default
#
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [enable debugging [default is no]]))

#
# FFT engine fftw3 ?
#
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

#
# FFT engine intel mkl ?
#
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
	     MKL_LDFLAGS="-L$with_mkl/lib/intel64"
	     fi])

#
# FFT engine is generic if fftw3 or mkl is not provided
#
if test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xcheck"; then
	AC_MSG_NOTICE([FFT engine fftw3 was provided])
	AC_SUBST(FFTENGINE, fftw3)
	AC_SUBST(FFTW3FCFLAGS, ${FFTW3_FCFLAGS})
	AC_SUBST(FFTW3LDFLAGS, ${FFTW3_LDFLAGS})
	if test "x$enable_single" = "xyes"; then
		AC_SUBST(FFTW3LDADD, -lfftw3f)
	else
		AC_SUBST(FFTW3LDADD, -lfftw3)
	fi
else
if test "x$with_mkl" != "xno" -a "x$with_mkl" != "xcheck"; then
	AC_MSG_NOTICE([FFT engine mkl was provided])
	AC_SUBST(FFTENGINE, mkl)
	AC_SUBST(MKLFCFLAGS, ${MKL_FCFLAGS})
        AC_SUBST(MKLLDFLAGS, ${MKL_LDFLAGS})
else
	AC_MSG_NOTICE([FFT engine is generic])
	AC_SUBST(FFTENGINE, generic)
fi
fi
AM_CONDITIONAL(HAVE_FFTW3, [test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xcheck"])
AM_CONDITIONAL(HAVE_MKL, [test "x$with_mkl" != "xno" -a "x$with_mkl" != "xcheck"])

#
# pfunit unit tests
#
AC_ARG_WITH(pfunit,
            [AS_HELP_STRING([--with-pfunit=PATH],[prefix directory for pfunit])],
            [if test "x$withval" = "x"; then
                with_pfunit=no
            fi],
            [with_pfunit=no])
if test "x$with_pfunit" != "xno"; then
        AC_MSG_NOTICE([Unit tests activated])
        AC_SUBST(PFUNITDIR, ${with_pfunit})
fi
AM_CONDITIONAL(HAVE_PFUNIT, [test "x$with_pfunit" != "xno"])

#
# Check the Fortran compiler
#
Fcompiler=none
#
# 1. gfortran ?
#
if test "x$Fcompiler" != "xok"; then
echo `$FC --version 2>&1 | head -1` | grep 'GNU Fortran' > /dev/null
if test "$?" = "0" ; then
	Fcompiler=ok
	fc_version="`$FC -dumpversion | cut -f1 -d.`"
	AC_MSG_NOTICE([GNU Fortran version $fc_version])
	if test "x$enable_single" = "xyes"; then
		FCFLAGS="-cpp"
	else
                FCFLAGS="-cpp -DDOUBLE_PREC"
		if test "x$enable_single_snapshot" = "xyes"; then
			FCFLAGS="$FCFLAGS -DSAVE_SINGLE"
		fi
        fi
	AS_IF([test "x$enable_debug" = "xyes"],
	      [FCFLAGS="$FCFLAGS -Og -g3 -ffpe-trap=invalid,zero -fcheck=all -fimplicit-none -DDEBUG"],
	      [FCFLAGS="$FCFLAGS -O3 -DNDEBUG"])
	OMPFLAG="-fopenmp"
fi
fi
#
# 2. ifort classic ?
#
if test "x$Fcompiler" != "xok"; then
echo `$FC --version 2>&1 | head -1` | grep 'ifort' > /dev/null
if test "$?" = "0"; then
	Fcompiler=ok
	fc_version="`$FC -dumpversion | cut -f1 -d.`"
	AC_MSG_NOTICE([Intel ifort version $fc_version])
        if test "x$enable_single" = "xyes"; then
                FCFLAGS="-fpp"
        else
                FCFLAGS="-fpp -DDOUBLE_PREC"
                if test "x$enable_single_snapshot" = "xyes"; then
                        FCFLAGS="$FCFLAGS -DSAVE_SINGLE"
                fi
        fi
	AS_IF([test "x$enable_debug" = "xyes"],
	      [FCFLAGS="$CFLAGS -g -O0 -traceback -check all -check nopointer -fpe0 -ftrapuv -DDEBUG"],
	      [FCFLAGS="$FCFLAGS -O3 -DNDEBUG"])
	OMPFLAG="-qopenmp"
fi
fi
#
# 3. nvfortran ?
#

#
# 4. TODO (Cray, IBM XL, ...)
#

#
# Provide the git version
# or the one defined in AC_INIT
#
if test -d ${srcdir}/.git; then
	X3DVERSION="`git $(srcdir) describe --tag --long --always`"
else
	X3DVERSION=${PACKAGE_VERSION}
fi
AC_SUBST(X3DVERSION)

#
# Combine default flags with the user ones
#
FCFLAGS="$FCFLAGS $user_FCFLAGS"
#LDFLAGS="$LDFLAGS $user_LDFLAGS"

#
# Propagate specific flag
#
AC_SUBST(OMPFLAG)

#
# Maintainer mode if a .git folder is found
# ("rebuild rules" for out of date Makefiles, configure, ...)
#
AS_IF([test -d ${srcdir}/.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

#
# Create libtool, needed for shared/static libraries.
#
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([dlopen])

AC_CONFIG_FILES([\
	Makefile \
	src/Makefile \
	decomp2d/Makefile \
	verification/Makefile \
	verification/lib/Makefile \
	verification/verif/Makefile])

AC_OUTPUT


#
#
# Listing at the end of configure
#
#
echo
echo
echo "Configuration options:"
echo
#
if test "x$enable_single" = "xyes"; then
	echo "   Precision: single"
else
	echo "   Precision: double"
if test "x$enable_single_snapshot" = "xyes"; then
	echo "   Snapshop precision: single"
else
	echo "   Snapshot precision: double"
fi
fi
#
if test "x$enable_debug" = "xyes"; then
	echo "   Build: debug"
else
	echo "   Build; production"
fi
#
if test "x$with_fftw3" != "xno" -a "x$with_fftw3" != "xcheck"; then
	echo "   FFT engine: FFTW3"
else
if test "x$with_mkl" != "xno" -a "x$with_mkl" != "xcheck"; then
	echo "   FFT engine: MKL"
else
	echo "   FFT engine: generic"
fi
fi
#
if test "x$with_pfunit" != "xno"; then
	echo "   unit tests activated (pFunit)"
else
	echo "   no unit test"
fi
echo
echo "Configuration is complete."
echo "Compile and install the code using:"
echo " make -j"
if test "x$with_pfunit" != "xno"; then
	echo " make -j check"
fi
echo " make -j install"
echo
echo "Please note that it is best to run configure outside of the source folder"
echo
echo "In case of issue, check the section \"Output variables\" in the file config.log"
echo
