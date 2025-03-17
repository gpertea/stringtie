# Copyright (c) 2003, 2005-2007, 2009-2010, 2013, 2019, 2022-2023
# Genome Research Ltd.
# Author(s): James Bonfield
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
#    Institute nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific
#    prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
# LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
lib_LTLIBRARIES = libhtscodecs.la

libhtscodecs_base_src = \
	pack.c \
	pack.h \
	rle.c \
	rle.h \
	fqzcomp_qual.c \
	fqzcomp_qual.h \
	rANS_static.c \
	rANS_static.h \
	rANS_byte.h \
	rANS_static4x16pr.c \
	rANS_static4x16.h \
	rANS_word.h \
	rANS_static32x16pr.c \
	rANS_static32x16pr.h \
	rANS_static16_int.h \
	permute.h \
	tokenise_name3.c \
	tokenise_name3.h \
	pooled_alloc.h \
	arith_dynamic.c \
	arith_dynamic.h \
	c_range_coder.h \
	c_simple_model.h \
	varint.h \
	htscodecs.c \
	htscodecs.h \
	htscodecs_endian.h \
	utils.c \
	utils.h

libhtscodecs_la_SOURCES = $(libhtscodecs_base_src)
libhtscodecs_la_LIBADD =
noinst_LTLIBRARIES =

# SIMD optional extras
if RANS_32x16_SSE4
noinst_LTLIBRARIES += librANS_static32x16pr_sse4.la
librANS_static32x16pr_sse4_la_SOURCES = rANS_static32x16pr_sse4.c
librANS_static32x16pr_sse4_la_CFLAGS = @MSSE4_1@ @MSSSE3@ @MPOPCNT@
libhtscodecs_la_LIBADD += librANS_static32x16pr_sse4.la
endif
if RANS_32x16_AVX2
noinst_LTLIBRARIES += librANS_static32x16pr_avx2.la
librANS_static32x16pr_avx2_la_SOURCES = rANS_static32x16pr_avx2.c
librANS_static32x16pr_avx2_la_CFLAGS = @MAVX2@
libhtscodecs_la_LIBADD += librANS_static32x16pr_avx2.la
endif
if RANS_32x16_AVX512
noinst_LTLIBRARIES += librANS_static32x16pr_avx512.la
librANS_static32x16pr_avx512_la_SOURCES = rANS_static32x16pr_avx512.c
librANS_static32x16pr_avx512_la_CFLAGS = @MAVX512@
libhtscodecs_la_LIBADD += librANS_static32x16pr_avx512.la
endif
if RANS_32x16_NEON
libhtscodecs_la_SOURCES += rANS_static32x16pr_neon.c
endif

libhtscodecs_la_LDFLAGS = -version-info @VERS_CURRENT@:@VERS_REVISION@:@VERS_AGE@ 
libhtscodecs_la_LIBADD += -lm

# Fuzz testing version of the library.  This is build using -fsanitize=fuzzer
# and defines FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION to prevent things like
# large allocations that would otherwise upset the fuzzer.
fuzz: $(EXTRA_LIBRARIES)

# Note that we build several libraries here, so we can get automake to
# use the right options for the various parts.
# See https://www.gnu.org/software/automake/manual/html_node/Per_002dObject-Flags.html
EXTRA_LIBRARIES = libcodecsfuzz.a libcodecsfuzz_sse4.a libcodecsfuzz_avx2.a libcodecsfuzz_avx512.a
libcodecsfuzz_a_SOURCES = $(libhtscodecs_base_src)
libcodecsfuzz_a_CFLAGS = -fsanitize=fuzzer -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
libcodecsfuzz_a-htscodecs.$(OBJEXT): version.h
libcodecsfuzz_sse4_a_SOURCES = rANS_static32x16pr_sse4.c
libcodecsfuzz_sse4_a_CFLAGS = $(libcodecsfuzz_a_CFLAGS) @MSSE4_1@ @MSSSE3@ @MPOPCNT@
libcodecsfuzz_avx2_a_SOURCES = rANS_static32x16pr_avx2.c
libcodecsfuzz_avx2_a_CFLAGS = $(libcodecsfuzz_a_CFLAGS) @MAVX2@
libcodecsfuzz_avx512_a_SOURCES = rANS_static32x16pr_avx512.c
libcodecsfuzz_avx512_a_CFLAGS = $(libcodecsfuzz_a_CFLAGS) @MAVX512@

version.h: force
	@ if `git describe 2>/dev/null >/dev/null`; then \
	    echo '#define HTSCODECS_VERSION_TEXT "'`git describe --match 'v[0-9]\.[0-9]*' --dirty|sed 's/^v//'`'"' > _version.h; \
	else \
	    echo '#define HTSCODECS_VERSION_TEXT "@PACKAGE_VERSION@"' > _version.h; \
	fi;
	@ if ( [ ! -e version.h ] || ! cmp -s version.h _version.h); then \
	    echo "Update version.h: `cat _version.h`"; \
	    mv _version.h version.h; \
	else \
	    rm _version.h; \
	fi

force:

.PHONY: force

# Manually added as the auto-dependency generation won't pick up a file that
# may not yet exist.
BUILT_SOURCES = version.h

distclean-local:
	rm -f version.h
