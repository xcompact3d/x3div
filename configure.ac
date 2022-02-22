# Name / version / website
AC_INIT([x3div], [alpha], [https://github.com/xcompact3d/x3div/])
AC_CONFIG_MACRO_DIR([./m4])

AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])

# Check for a MPI library
AX_PROG_FC_MPI([],[mpifound=yes],[mpifound=no])
if test "x$mpifound" = "xyes"; then
	AC_MSG_NOTICE([MPI detected : $FC])
fi
if test "x$mpifound" = "xno"; then
	AC_MSG_FAILURE([Fortran MPI not detected])
fi

# Debug mode, off by default
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [enable debugging (default is no)]))

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
	FCFLAGS="-cpp"
	AS_IF([test "x$enable_debug" = "xyes"],
	      [FCFLAGS="-Og -g3 -ffpe-trap=invalid,zero -fcheck=bounds -fimplicit-none -DDEBUG $FCFLAGS"],
	      [FCFLAGS="-O3 -march=native $FCFLAGS"])
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

# Maintener mode if a .git folder is found
AS_IF([test -d ${srcdir}/.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-static])

AC_CONFIG_FILES([Makefile src/Makefile decomp2d/Makefile])

AC_OUTPUT
