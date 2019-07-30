AC_INIT([global], [1], [gmarcais@umd.edu])
AC_CANONICAL_HOST
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests color-tests])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([.])
AC_CONFIG_HEADERS([config.h])

# Change default compilation flags
AS_IF([test -z "$CXXFLAGS"], [cxxflags_given=no], [cxxflags_given=yes])
AC_LANG(C++)
AC_PROG_CC
AC_PROG_CXX

#
LT_PREREQ
LT_INIT

#
#AC_PROG_LIBTOOL

# Checks for libraries.
CXXFLAGS="-std=c++0x $CXXFLAGS"
AS_IF([test x$cxxflags_given = xno], [CXXFLAGS="-DNDEBUG $CXXFLAGS"])
EXTRA_CXXFLAGS="-std=c++0x"
AC_SUBST([EXTRA_CXXFLAGS])
ALL_CXXFLAGS="$EXTRA_CXXFLAG"
AC_SUBST([ALL_CXXFLAGS])

AC_CHECK_LIB([pthread], [pthread_create])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define to 1 if you have the 'clock_gettime' function])])

# Flags for Jellyfish
AC_SUBST([JELLYFISH2_0_CFLAGS], ['-I$(includedir)/jellyfish-1'])
AC_SUBST([JELLYFISH2_0_LIBS], ['-L$(libdir) -ljellyfish-2.0'])
AC_SUBST([JELLYFISH], ['$(bindir)/jellyfish'])
AC_CHECK_TYPE([__int128],
              [AC_DEFINE([HAVE_INT128], [1], [Define if type __int128 is supported])])

# Check for openmp
AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--disable-openmp], [Disable OpenMP, even if available])])
AS_IF([test "x$enable_openmp" != "xno"],
      # OpenMP enabled. Check if available
      [AX_OPENMP([AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])]
      [EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS $OPENMP_CXXFLAGS"],
      # OpenMP disabled. Ignore pragma warnings
      [EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS -Wno-unknown-pragmas"])

# Check for Boost, and install if asked to
#AS_IF([test x$BOOST_ROOT = xinstall],
#      AS_IF([test ! -f include/boost/version.hpp],
#            [AS_BOX([Installing Boost 1.59.0 from the Internet into include/boost])]
#            [AS_MKDIR_P([include])]
#            [curl -L http://downloads.sourceforge.net/project/boost/boost/1.59.0/boost_1_59_0.tar.bz2 | tar -jx -C include --strip-components=1 boost_1_59_0/boost]
#            [AS_IF([test x$? != x0], [AC_MSG_FAILURE([failed to install Boost from the Internet])])])
#      [BOOST_ROOT=`pwd`/include])
#
#BOOST_REQUIRE([1.46.0],
#              [AS_BOX([Could not find a working installation of Boost. Set BOOST_ROOT to the path where the Boost headers are installed or set BOOST_ROOT=install to have it downloaded from the Internet and installed locally. For example: BOOST_ROOT=install ./install.sh])]
#              [AC_MSG_FAILURE([cannot find Boost])])
#BOOST_FIND_HEADER([boost/icl/interval_set.hpp])

# For ufasta, disable regex unconditionally
AM_CONDITIONAL([HAVE_BOOST_REGEX], [false])

# Flags for samtools
AC_SUBST([SAMTOOLS_CFLAGS], ['-I$(includedir)/samtools'])
AC_SUBST([SAMTOOLS_LIBS], ['-L$(libdir) -lbam'])

# On MacOS X, use _NSGetExecutablePath to find path to own executable
AC_MSG_CHECKING([for _NSGetExecutablePath])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mach-o/dyld.h>]],
                                   [[_NSGetExecutablePath(0, 0);]])],
                  [AC_MSG_RESULT([yes])]
                  [AC_DEFINE([HAVE_NSGETEXECUTABLEPATH], [1], [Used to find executable path on MacOS X])],
                  [AC_MSG_RESULT([no])])

# For MUMmer, check for std::thread::sleep_for
AC_SUBST([EXTRA_CXXFLAGS])
AC_MSG_CHECKING([if std::thread::sleep_for is supported])
try_nano_sleep=
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <thread>]],
                                   [[std::this_thread::sleep_for(std::chrono::seconds(1))]])],
                  [AC_MSG_RESULT([yes])], [try_nano_sleep=yes])
AS_IF([test x$try_nano_sleep = xyes],
      [save_cxxflags="$CXXFLAGS"]
      [CXXFLAGS="-D_GLIBCXX_USE_NANOSLEEP $CXXFLAGS"]
      [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <thread>]],
                                         [[std::this_thread::sleep_for(std::chrono::seconds(1))]])],
                         [AC_MSG_RESULT([yes])]
                         [EXTRA_CXXFLAGS="-D_GLIBCXX_USE_NANOSLEEP $EXTRA_CXXFLAGS"],
                         [AC_MSG_RESULT([no])]
                         [AC_MSG_FAILURE([Compiler does not support std::this_thread::sleep_for])])]
      [CXXFLAGS="$save_cxxflags"])

# Checks for library functions.
AC_FUNC_MMAP
AC_CHECK_FUNCS([mremap])

# Check for zlib
AC_CHECK_HEADER([zlib.h], [], [AC_MSG_FAILURE([Cannot find zlib.h header. Install zlib development files (e.g. zlib1g-dev)])])
AC_CHECK_LIB([z], [zlibVersion], [], [AC_MSG_FAILURE([Cannot find zlib library. Install zlib libray])])

# Check bzip2 library
AC_CHECK_HEADER([bzlib.h], [], [AC_MSG_FAILURE([Cannot find bzilb.h header. Install bzip2 development files (e.g. bzip2-devel)])])
AC_CHECK_LIB([bz2], [BZ2_bzlibVersion], [], [AC_MSG_FAILURE([Cannot find bzip2 library. Install bzip2 library])])


## Use valgrind to check memory allocation with mmap
#AC_ARG_ENABLE([valgrind],
#              [AS_HELP_STRING([--enable-valgrind], [Instrument mmap memory allocation with valgrind])])
#AS_IF([test "x$enable_valgrind" = xyes],
#      [AC_DEFINE([HAVE_VALGRIND], [1], [Define is using Valgrind])]
#      [PKG_CHECK_MODULES([VALGRIND], [valgrind >= 1.8.0])])

# Check for builtin functions
AC_MSG_CHECKING([for __builtin_prefetch])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[int x = 0;]],
                                [[__builtin_prefetch(&x);]])],
               [AC_MSG_RESULT([yes])]
               [AC_DEFINE([HAVE_BUILTIN_PREFETCH], [1], [Defined if __builtin_prefetch is supported])],
               [AC_MSG_RESULT([no])])

# For Quorum
AM_CONDITIONAL([RELATIVE_PATHS], [true])
AM_CONDITIONAL([HAVE_RELATIVE_JF_PATH], [false])

# For Samtools, don't compile ncurses stuff
AC_SUBST([CURSES_LIB], [0])

# For CA
AC_MSG_CHECKING([for architecture])
OSTYPE=`uname`
OSVERSION=`uname -r`
MACHINETYPE=`uname -m`
AS_CASE([$MACHINETYPE],
        ["x86_64"], [MACHINETYPE=amd64],
        ["Power Macintosh"], [MACHINETYPE=ppc])
AS_IF([test x$OSTYPE = xSunOS],
      [MACHINETYPE=`uname -p`]
      [AS_IF([text x$MACHINETYPE = xsparc],
             [MACHINETYPE=sparc`/usr/bin/isainfo -b`])])
AC_SUBST([OSTYPE])
AC_SUBST([OSVERSION])
AC_SUBST([MACHINETYPE])
AC_SUBST([LOCAL_NAME], [$OSTYPE-$MACHINETYPE])
AC_MSG_RESULT([$LOCAL_NAME])

AM_CONDITIONAL([SITE_LOCAL], [true])
AM_CONDITIONAL([SITE_JCVI], [false])
AM_CONDITIONAL([SITE_BRI], [false])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])

# Check for swig and enable Perl bindings
# --enable-swig
AC_ARG_ENABLE([swig],
              [AC_HELP_STRING([--enable-swig], [enable development of swig binding])])
AS_IF([test x$enable_swig = xyes],
      [AX_PKG_SWIG([3.0.0], [], [AC_MSG_FAILURE([SWIG version 3 is required])])])
AS_IF([test -n "$SWIG"],
      [SWIG_ENABLE_CXX])
AM_CONDITIONAL([HAVE_SWIG], [test -n "$SWIG"])
AC_SUBST([MAYBE_SWIG], [swig])
# Conf Perl
PERL_EXT_LIB='$(libdir)/perl'
AX_PERL_EXT
AM_CONDITIONAL([PERL_BINDING], [true])
# Disable others
AM_CONDITIONAL([PYTHON_BINDING], [false])
AM_CONDITIONAL([RUBY_BINDING], [false])

# Output configuration files
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([jellyfish/jellyfish-2.0.pc])
m4_map_args_w([ufasta quorum SuperReadsR jellyfish],
              [AC_CONFIG_FILES(], [/Makefile)])
AC_OUTPUT
