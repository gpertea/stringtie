ACLOCAL_AMFLAGS=-I m4

AM_CPPFLAGS = -I$(srcdir)/include -I$(srcdir) $(JELLYFISH2_0_CFLAGS)
AM_CXXFLAGS = -Wall -g -O2 -std=c++0x
AM_LDFLAGS = $(STATIC_FLAGS)
LDADD = $(JELLYFISH2_0_LIBS)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	@echo yaggo $< $@
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
YAGGO_SOURCES = src/error_correct_reads_cmdline.hpp	\
                src/create_database_cmdline.hpp		\
                src/merge_mate_pairs_cmdline.hpp	\
                src/split_mate_pairs_cmdline.hpp

BUILT_SOURCES = $(YAGGO_SOURCES)
noinst_HEADERS = $(YAGGO_SOURCES)
CLEANFILES =
EXTRA_DIST =

bin_PROGRAMS = quorum_error_correct_reads quorum_create_database	\
               merge_mate_pairs split_mate_pairs

quorum_error_correct_reads_SOURCES = src/error_correct_reads.cc	\
                                     src/err_log.cc

quorum_create_database_SOURCES = src/create_database.cc

merge_mate_pairs_SOURCES = src/merge_mate_pairs.cc

split_mate_pairs_SOURCES = src/split_mate_pairs.cc

noinst_HEADERS += src/error_correct_reads.hpp				\
                  src/error_correct_reads.hpp src/verbose_log.hpp	\
                  src/kmer.hpp src/mer_database.hpp src/err_log.hpp

noinst_HEADERS += include/gzip_stream.hpp include/misc.hpp	\
                  include/jflib/locks_pthread.hpp		\
                  include/jflib/pool.hpp			\
                  include/jflib/multiplexed_io.hpp


###############################
# Adapter database generation #
###############################
dist_data_DATA = data/adapter.jf
DISTCLEANFILES = data/adapter.jf
SUFFIXES = .jf .fa

.fa.jf:
	$(AM_V_GEN)mkdir -p $(@D); $(JELLYFISH) count -m 24 -s 5k -C -t 10 -o $@ $<

perl_scripts = quorum
libs_scripts = $(patsubst %,.libs/%,$(perl_scripts))

bin_SCRIPTS = $(libs_scripts)
noinst_SCRIPTS = $(perl_scripts)
CLEANFILES += $(libs_scripts) $(perl_scripts
EXTRA_DIST += $(patsubst %,src/%.in,$(perl_scripts))

if RELATIVE_PATHS
quorum_bindir=$$dir
else
quorum_bindir=$(bindir)
endif

if HAVE_RELATIVE_JF_PATH
quorum_jellyfish=$$dir/$(RELATIVE_JF_PATH)
local_quorum_jellyfish = $(shell realpath "$(bindir)/$(RELATIVE_JF_PATH)")
else
quorum_jellyfish=$(JELLYFISH)
local_quorum_jellyfish = $(quorum_jellyfish)
endif


inst_subst = sed -e 's|[@]JELLYFISH[@]|$(quorum_jellyfish)|g' \
                 -e 's|[@]bindir[@]|$(quorum_bindir)|g' \
	         -e 's|[@]JF_LIB_PATH[@]|$(JF_LIB_PATH)|g' \
	         -e 's|[@]PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	         -e 's|[@]JELLYFISH_VERSION[@]|$(JELLYFISH_VERSION)|g'

local_subst = sed -e 's|[@]JELLYFISH[@]|$(quorum_jellyfish)|g' \
	          -e 's|[@]bindir[@]|$(abs_builddir)|g' \
	          -e 's|[@]JF_LIB_PATH[@]|$(JF_LIB_PATH)|g' \
	          -e 's|[@]PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	          -e 's|[@]JELLYFISH_VERSION[@]|$(JELLYFISH_VERSION)|g'

.libs/%: src/%.in Makefile
	@mkdir -p $(dir $@) && $(inst_subst) < $< > $@
%: src/%.in Makefile
	$(AM_V_GEN)$(local_subst) < $< > $@ && chmod a+x $@

##############
# Unit tests #
##############
TESTS = all_tests
check_PROGRAMS = all_tests query_mer_database histo_mer_database

all_tests_SOURCES = unit_tests/test_mer_database.cc
all_tests_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/unit_tests/gtest/include -I$(srcdir)/unit_tests
all_tests_LDADD = libgtest_main.la $(LDADD)
noinst_HEADERS += unit_tests/test_misc.hpp

query_mer_database_SOURCES = src/query_mer_database.cc
histo_mer_database_SOURCES = src/histo_mer_database.cc

include gtest.mk
-include $(srcdir)/local.mk
