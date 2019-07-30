ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)
AM_CXXFLAGS = $(ALL_CXXFLAGS) -g -O3
AM_LDFLAGS =

bin_PROGRAMS = ufasta
noinst_HEADERS = $(YAGGO_SRCS) common.hpp
BUILT_SOURCES = $(YAGGO_SRCS)
CLEANFILES =
DISTCLEANFILES = $(BUILT_SOURCES)
EXTRA_DIST =

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO   " $@;
.yaggo.hpp:
#	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
	$(V_YAGGO)$(YAGGO) -o $@ $<

YAGGO_SRCS = # Append all file to be built by yaggo

##################
# Build binaries #
##################
SUB_CMDS = one.cc sizes.cc head.cc tail.cc rc.cc n50.cc extract.cc	\
           format.cc sort.cc dsort.cc split.cc

if HAVE_BOOST_REGEX
AM_LDFLAGS += -lboost_regex
SUB_CMDS += hgrep.cc dgrep.cc
endif

ufasta_SOURCES = ufasta.cc $(SUB_CMDS)
YAGGO_SRCS += $(SUB_CMDS:.cc=_cmdline.hpp)

#########
# Tests #
#########
tests/testsh: tests/testsh.in
	@mkdir -p $(dir $@)
	$(AM_V_GEN)sed -e 's,[@]SRCDIR[@],$(abs_srcdir),g' \
	               -e 's,[@]BUILDDIR[@],$(abs_builddir),g' $< > $@
	@chmod a+x $@

check_SCRIPTS = tests/testsh
CLEANFILES += $(check_SCRIPTS)
EXTRA_DIST += tests/testsh.in

TESTS = tests/test_head.sh tests/test_tail.sh tests/test_rc.sh		\
        tests/test_one.sh tests/test_hgrep.sh tests/test_n50.sh		\
        tests/test_sizes.sh tests/test_extract.sh			\
        tests/test_format.sh tests/test_sort.sh tests/test_dsort.sh	\
        tests/test_dgrep.sh tests/test_split.sh
EXTRA_DIST += $(TESTS) tests/test1.fasta
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = ./tests/testsh
