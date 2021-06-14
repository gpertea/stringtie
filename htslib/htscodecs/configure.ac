dnl Process this file with autoconf to produce a configure script.
AC_INIT(htscodecs, 1.0)

# Some functions benefit from -O3 optimisation, so if the user didn't
# explicitly set any compiler flags, we'll plump for O3.
_=${CFLAGS="-g -O3"}

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE

# For my own benefit, here's a summary of how versioning should work.
#
# Library versioning, not to be confused with the release number:
# Update VERS_CURRENT when an interface changes in an incompatible manner.
#     eg. change of struct, removing a function, changing a prototype.
#     Also increment this when we add new functions in a backwards compatible
#     manner, but update VERS_AGE too (see below).
#
# Update VERS_REVISION whenever VERS_CURRENT hasn't changed, but the source
#     code has. (Otherwise reset revision to 0).
#
# Update VERS_AGE if new interfaces have been added, causing VERS_CURRENT to
#        be incremented too. The idea is that if the previous values are
#        VERS_CURRENT=2, VERS_AGE=0 then after adding a new function we
#        get VERS_CURRENT=3, VERS_AGE=1 (ie ABI compat with vers 2 or 3).
#
# Also see:
# http://www.gnu.org/software/libtool/manual/libtool.html#Updating-version-info

# Example: We have htscodecs 1.11.6 implementing ABI vers 0, rev 6, age 0.
# Our library has an soname of libhtscodecs.so.0 and files as follows:
#       (SONAME      libhtscodecs.so.0)
# 	libhtscodecs.so -> libhtscodecs.so.0.0.6
#       libhtscodecs.so.0 -> libhtscodecs.so.0.0.6
#       libhtscodecs.so.0.0.6
#
# We then release htscodecs 1.12.0 with modified structs, so it's incompatible.
# We bump to ABI vers 1, rev 0, age 0 - giving:
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.0.0
#       libhtscodecs.so.1 -> libhtscodecs.so.1.0.0
#       libhtscodecs.so.1.0.0
#
# If in 1.12.1 to 1.12.5 we then perform several successive minor tweaks (bug
# fix, optimisation, etc) that do not change the ABI at all, we may
# end up with ABI vers 1, rev 5, age 0:
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.0.5
#       libhtscodecs.so.1 -> libhtscodecs.so.1.0.5
#       libhtscodecs.so.1.0.5
#
# And assume in 1.13.0 we add an extra public function. Code linked
# against 1.12.x still works (aka library version 1) meaning our
# soname doesn't change, but we do update ABI version: vers 2, rev 0, age 1.
#       (SONAME      libhtscodecs.so.1)
# 	libhtscodecs.so -> libhtscodecs.so.1.1.0
#       libhtscodecs.so.1 -> libhtscodecs.so.1.1.0
#       libhtscodecs.so.1.1.0

VERS_CURRENT=2
VERS_REVISION=1
VERS_AGE=0
AC_SUBST(VERS_CURRENT)
AC_SUBST(VERS_REVISION)
AC_SUBST(VERS_AGE)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

# # by Marcelo Magallon <mmagallo@efis.ucr.ac.cr>
# # Turn around -rpath problem with libtool 1.0c
# # This define should be improbable enough to not conflict with anything.
# echo host=$host
# case ${host} in
#   *-linux-gnu | *-apple-darwin*)
#     AC_MSG_RESULT([Fixing libtool for -rpath problems.])
#     sed < libtool > libtool-2 \
#     's/^hardcode_libdir_flag_spec.*$'/'hardcode_libdir_flag_spec=""/'
#     mv libtool-2 libtool
#     chmod 755 libtool
#   ;;
# esac

dnl Turn on C compiler warnings.
VL_PROG_CC_WARNINGS

dnl Check for libraries
ZLIB_CHECK_CONFIG(,,[have_zlib=yes], [have_zlib=no])

AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [
	LIBS="-lbz2 $LIBS"
	AC_DEFINE([HAVE_LIBBZ2],1,[Define to 1 if you have the libbz2 library.])])

dnl AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [
dnl 	LIBS="-llzma $LIBS"
dnl 	AC_DEFINE([HAVE_LIBLZMA],1,[Define to 1 if you have the liblzma library.])])
dnl 
dnl # Experimental
dnl AC_CHECK_LIB([bsc], [bsc_compress], [
dnl 	LIBS="-lbsc $LIBS"
dnl 	AC_DEFINE([HAVE_LIBBSC],1,[Define to 1 if you have the libbsc library.])])

AX_LIBDEFLATE

dnl Checks for header files.
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(fcntl.h limits.h unistd.h malloc.h)
AC_CHECK_HEADER(zlib.h)

dnl Checks for library functions.
AC_SEARCH_LIBS([pthread_join], [pthread])

AC_CONFIG_FILES([Makefile htscodecs/Makefile tests/Makefile])
AC_OUTPUT
