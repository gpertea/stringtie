dnl Process this file with autoconf to produce a configure script.
AC_INIT(htscodecs, 1.4.0)

# Some functions benefit from -O3 optimisation, so if the user didn't
# explicitly set any compiler flags, we'll plump for O3.
_=${CFLAGS="-g -O3"}

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AM_EXTRA_RECURSIVE_TARGETS([fuzz])

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

VERS_CURRENT=3
VERS_REVISION=1
VERS_AGE=1
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
AC_ARG_ENABLE([bz2],
  [AS_HELP_STRING([--disable-bz2],
                  [omit support for BZ2 within Arith streams])],
  [], [enable_bz2=yes])

if test "$enable_bz2" != no; then
  bz2_devel=ok
  AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [
	LIBS="-lbz2 $LIBS"
	AC_DEFINE([HAVE_LIBBZ2],1,[Define to 1 if you have the libbz2 library.])], [bz2_devel=missing])
  if test "$bz2_devel" != "ok"; then
    AC_MSG_ERROR([libbzip2 development files not found.

This is needed for some options of the arithmetic coder.  It can be disabled
with configure --disable-bz2, but be aware some data streams will become
unable to be decompressed and arith tests will fail.])
  fi
fi

dnl AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [
dnl 	LIBS="-llzma $LIBS"
dnl 	AC_DEFINE([HAVE_LIBLZMA],1,[Define to 1 if you have the liblzma library.])])
dnl 
dnl # Experimental
dnl AC_CHECK_LIB([bsc], [bsc_compress], [
dnl 	LIBS="-lbsc $LIBS"
dnl 	AC_DEFINE([HAVE_LIBBSC],1,[Define to 1 if you have the libbsc library.])])

dnl Count parts needed to build rANS_static32x16pr_sse4.c
sse4_prerequisites=""

dnl Check if we can use our SSSE3 implementations of rANS 32x16 codec.
AX_CHECK_COMPILE_FLAG([-mssse3], [
        MSSSE3=-mssse3
	sse4_prerequisites="o$sse4_prerequisites"
	AC_SUBST([MSSSE3])
        AC_DEFINE([HAVE_SSSE3],1,[Defined to 1 if the compiler can issue SSSE3 instructions.])
], [], [], [AC_LANG_PROGRAM([[
	  #include "x86intrin.h"
	]],[[
	  __m128i a = _mm_set_epi32(1, 2, 3, 4), b = _mm_set_epi32(4, 3, 2, 1);
	  __m128i c = _mm_shuffle_epi8(a, b);
	  return *((char *) &c);
	]])
])

dnl Check if we can use popcnt instructions
AX_CHECK_COMPILE_FLAG([-mpopcnt], [
        MPOPCNT=-mpopcnt
	sse4_prerequisites="o$sse4_prerequisites"
	AC_SUBST([MPOPCNT])
        AC_DEFINE([HAVE_POPCNT],1,[Defined to 1 if the compiler can issue popcnt instructions.])
], [], [], [AC_LANG_PROGRAM([[
	  #include "x86intrin.h"
	]],[[
	  unsigned int i = _mm_popcnt_u32(1);
	  return i != 1;
	]])
])

dnl Check if we can use our SSE4.1 too.  This *may* always imply SSSE3?
dnl It may be easier just to target an old era of cpu than -mssse3 -msse4.1
dnl -mpopcnt.  Eg -march=nehalem.  I don't know how wide spread that is.
AX_CHECK_COMPILE_FLAG([-msse4.1], [
        MSSE4_1=-msse4.1
	sse4_prerequisites="o$sse4_prerequisites"
	AC_SUBST([MSSE4_1])
        AC_DEFINE([HAVE_SSE4_1],1,[Defined to 1 if the compiler can issue SSE4.1 instructions.])
], [], [], [AC_LANG_PROGRAM([[
	  #include "x86intrin.h"
	]],[[
	  __m128i a = _mm_set_epi32(1, 2, 3, 4), b = _mm_set_epi32(4, 3, 2, 1);
	  __m128i c = _mm_max_epu32(a, b);
	  return *((char *) &c);
	]])
])
AM_CONDITIONAL([RANS_32x16_SSE4],[test "x$sse4_prerequisites" = "xooo"])

dnl Check if we can use our AVX2 implementations.
AX_CHECK_COMPILE_FLAG([-mavx2], [
        MAVX2=-mavx2
	AC_SUBST([MAVX2])
        AC_DEFINE([HAVE_AVX2],1,[Defined to 1 if the compiler can issue AVX2 instructions.])
], [], [], [AC_LANG_PROGRAM([[
	  #include "x86intrin.h"
	]],[[
	  __m256i a = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
	  __m256i b = _mm256_add_epi32(a, a);
	  long long c = _mm256_extract_epi64(b, 0);
	  return (int) c;
	]])
])
AM_CONDITIONAL([RANS_32x16_AVX2],[test "x$MAVX2" != "x"])

dnl Check also if we have AVX512.  If so this overrides AVX2
AX_CHECK_COMPILE_FLAG([-mavx512f], [
        MAVX512=-mavx512f
	AC_SUBST([MAVX512])
        AC_DEFINE([HAVE_AVX512],1,[Defined to 1 if the compiler can issue AVX512F instructions.])
], [], [], [AC_LANG_PROGRAM([[
	  #include "x86intrin.h"
	]],[[
	  __m512i a = _mm512_set1_epi32(1);
	  __m512i b = _mm512_add_epi32(a, a);
	  return *((char *) &b);
	]])
])
AM_CONDITIONAL([RANS_32x16_AVX512],[test "x$MAVX512" != "x"])

dnl Detect ARM Neon availability
AC_CACHE_CHECK([whether C compiler supports ARM Neon], [htscodecs_cv_have_neon], [
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[
      #include "arm_neon.h"
    ]], [[
      int32x4_t a = vdupq_n_s32(1);
      int32x4_t b = vaddq_s32(a, a);
      return *((char *) &b);
    ]])], [htscodecs_cv_have_neon=yes], [htscodecs_cv_have_neon=no])])
AM_CONDITIONAL([RANS_32x16_NEON],[test "$htscodecs_cv_have_neon" = yes])

AC_SUBST([HTSCODECS_SIMD_SRC])

AX_LIBDEFLATE

dnl Checks for header files.
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(fcntl.h limits.h unistd.h malloc.h)
AC_CHECK_HEADER(zlib.h)

dnl Checks for library functions.
AX_SEARCH_LIBS_REV([pthread_join], [thr pthread])

AC_CONFIG_FILES([Makefile htscodecs/Makefile tests/Makefile])
AC_OUTPUT
