Htscodecs
=========

See the NEWS file for a list of updates and version details.

[![Build Status](https://api.cirrus-ci.com/github/jkbonfield/htscodecs.svg?branch=master)](https://cirrus-ci.com/github/jkbonfield/htscodecs)

This repository implements the custom CRAM codecs used for "EXTERNAL"
block types.  These consist of two variants of the rANS codec (8-bit
and 16-bit renormalisation, with run-length encoding and bit-packing
also supported in the latter), a dynamic arithmetic coder, and custom
codecs for name/ID compression and quality score compression derived
from fqzcomp.


They come with small command line test tools to act as both
compression exploration programs and as part of the test harness.


Building
--------

If building from git, you'll need to recreate the configure script
using autoconf.  "autoreconf -i" should work if you have the
appropriate tools.

From then on, it follows the normal "./configure; make" or
"mkdir build; cd build; ../configure; make" rule.

The library can be used as a git sub-module or as a completely
separate entity.  If you are attempting to make use of these codecs
within your own library, such as we do within Staden io_lib, it may be
useful to configure this with `--disable-shared --with-pic'.


Testing
-------

There is a "make check" rule.  If you're using a modern clang you can
also cd to the tests directory and do "make fuzz" to build some fuzz
testing targets, but you'll likely need to modify Makefile.am first as
this has some hard-coded local paths.

We also provide test data and some command line tools to demonstrate
usage of the compression codecs.  These are in the tests directory
also. Example usage:

    ./fqzcomp_qual -s 1 < dat/q40+dir > /tmp/q40.comp
    ./fqzcomp_qual  -d < /tmp/q40.comp > /tmp/q40.uncomp
    awk '{print $1}' dat/q40+dir | md5sum;    # f91473032dd6942e72abec0868f17161
    awk '{print $1}' /tmp/q40.uncomp | md5sum;# f91473032dd6942e72abec0868f17161

The fqzcomp test format is one quality values per line, with an
optional additional parameter (0 or 1) to indicate READ1 or READ2 flag
status.

There is a larger set of test data in the htscodecs-corpus repository
(https://github.com/jkbonfield/htscodecs-corpus).  If this is cloned
into the tests subdirectory of htscodecs then the htscodecs "make
check" will also use that larger data set for testing purposes.


API
---

Many functions just take an input buffer and size and return an output
buffer, setting *out_size with the decoded size.  NULL is returned for
error.  This buffer is malloced and is expected to be freed by the
caller.  These are the *`compress` and *`uncompress` functions.

A second variant sometimes exists where the output buffer is
optionally allocated by the caller (it may be NULL in which case it
has the same operation as above).  If specified, `*out_size` must also
be set to the allocated size of `out`.  These are the `compress_to`
and `uncompress_to` functions.

The compress size sometimes needs additional options.  For the rANS
and arithmetic coder this is the "order".  Values of 0 and 1 are
simple order-0 and order-1 entropy encoder, but this is a bit field
and the more advanced codecs have additional options to pass in order
(so it should really be renamed to flags).  See below.  Fqzcomp
requires more input data - also see below.  In all cases, sufficient
information is stored in the compressed byte stream such that the
decompression will work without needing these input paramaters.

Finally the various `compress_bound` functions give the size of buffer
needed to be allocated when compressing a block of data.


### Static rANS 4x8 (introduced in CRAM v3.0)

```
#include "htscodecs/rANS_static.h"

unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
                             unsigned int *out_size, int order);
unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
                               unsigned int *out_size);
```

This is the earlier rANS entropy encoder using 4 rANS states and 8-bit
renormalisation, with Order-0 and Order-1 entropy models.

No (un)compress_to functions exist for this older codec.


### Static rANS 4x16 and 32x16 with bit-pack/RLE (CRAM v3.1):

```
#include "htscodecs/rANS_static4x16.h"

#define RANS_ORDER_X32    0x04  // 32-way unrolling instead of 4-way
#define RANS_ORDER_STRIPE 0x08  // N streams for every Nth byte (N==order>>8)
#define RANS_ORDER_NOSZ   0x10  // Don't store the original size
#define RANS_ORDER_CAT    0x20  // Nop; for tiny data segments
#define RANS_ORDER_RLE    0x40  // Run length encoding
#define RANS_ORDER_PACK   0x80  // Pack 2,4,8 or infinite symbols into a byte.

unsigned int rans_compress_bound_4x16(unsigned int size, int order);
unsigned char *rans_compress_to_4x16(unsigned char *in,  unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size,
                                     int order);
unsigned char *rans_compress_4x16(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size, int order);
unsigned char *rans_uncompress_to_4x16(unsigned char *in,  unsigned int in_size,
                                       unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_4x16(unsigned char *in, unsigned int in_size,
                                    unsigned int *out_size);
```

This is a faster version with 16-bit renormalisation and optional
transforms (RLE, small alphabet bit-packing, and interleaving of N
streams for e.g. 32-bit integer compression).  Additionally the
`order` field may include bit `RANS_ORDER_X32` in which case a 32-way
unrolled version will be used instead, with automatic CPU detection
and dispatching to an appropriate SIMD implementation if available.

### Adaptive arithmetic coding (CRAM v3.1):

```
#include "htscodecs/arith_dynamic.h"

unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
                              unsigned int *out_size, int order);

unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
                                unsigned int *out_size);

unsigned char *arith_compress_to(unsigned char *in,  unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size,
                                 int order);

unsigned char *arith_uncompress_to(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int *out_sz);

unsigned int arith_compress_bound(unsigned int size, int order);
```

These reuse the same `RANS_ORDER` bit fields and abilities above with
the exception of X32 as there is currently no unrolling of this code.

### Name tokeniser (CRAM v3.1):

```
#include "htscodecs/tokenise_name3.h"

uint8_t *encode_names(char *blk, int len, int level, int use_arith,
                      int *out_len, int *last_start_p);

uint8_t *decode_names(uint8_t *in, uint32_t sz, uint32_t *out_len);
```

This differs to the general purpose entropy encoders as it takes a
specific type of data.  The names should be newline or nul separated
for `encode_names`.  `decode_names` will alway return nul terminated
names, so you may need to swap these to newlines if you do round-trip
tests.

The compression level controls how hard it tries to find the optimum
compression method per internal token column.  By default it'll use
the rANS 4x16 codec, but with non-zero `use_arith` it'll use the
adaptive arithmetic coder instead.

If non-NULL, last_start_p can be used to point to a partial name if an
arbitrary block of names were supplied that don't end of a whole read
name. (Is this useful?  Probably not.)


### FQZComp Qual (CRAM v3.1):


```
#include "htscodecs/fqzcomp_qual.h"

#define FQZ_FREVERSE 16
#define FQZ_FREAD2 128

typedef struct {
    int num_records;
    uint32_t *len;    // of size num_records
    uint32_t *flags;  // of size num_records
} fqz_slice;

char *fqz_compress(int vers, fqz_slice *s, char *in, size_t uncomp_size,
                   size_t *comp_size, int strat, fqz_gparams *gp);
char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size,
                     int *lengths, int nlengths);
```

This is derived from the quality compression in fqzcomp.  The input
buffer is a concatenated block of quality strings, without any
separator.  In order to achieve maximum compression it needs to know
where these separators are, so they must be passed in via the
`fqz_slice` struct.

The summation of length fields should match the input uncomp_size
field.  Note the len fields may not actually be the length of the
original sequences as some CRAM features may additional quality values
(eg the "B" feature).

It can also be beneficial to supply per-record flags so fqzcomp can
determine whether orientation (complement strand) helps and whether
the READ1 vs READ2 quality distributions differ.  These are just
sub-fields from BAM FLAG.

The fqz_gparams will normally be passed in as NULL and the encoder
will automatically select parameters.  If you wish to fine tune the
compression methods, see the fqz_params and fqz_gparams structures in
the header file.  You may also find the fqz_qual_stats() utility
function helpful for gathering statistics on your quality values.

For decompression, the lengths array is optional and may be specified
as NULL.  If passed in, it must be of size nlengths and it will be
filled out with the decoded length of each quality string.  Note
regardless of whether lengths is NULL or not, the buffer returned will
be concatenated values so there is no way to tell where one record
finishes and the next starts.  (CRAM itself knows this via other means.)
