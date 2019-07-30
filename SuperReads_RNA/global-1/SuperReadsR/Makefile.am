ACLOCAL_AMFLAGS = -I m4
CXXFLAGS_ALL = -I$(srcdir)/include -I$(srcdir) -Wall -std=c++0x $(VALGRIND_CFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_ALL) $(JELLYFISH2_0_CFLAGS)
LDADD_ALL = -lrt -lpthread libutils.a
LDADD = $(LDADD_ALL) $(JELLYFISH2_0_LIBS)
AM_LDFLAGS = $(STATIC_FLAGS)
AM_DEV_CXXFLAGS = $(DEV_CXXFLAGS)

noinst_HEADERS = $(YAGGO_SOURCES)
bin_PROGRAMS =
dist_bin_SCRIPTS =
bin_SCRIPTS =
data_DATA =
BUILT_SOURCES = $(YAGGO_SOURCES)
CLEANFILES =
DISTCLEANFILES = $(BUILT_SOURCES)

EXTRA_DIST =

##############################
# jflib
##############################
JFLIB_I = include/jflib
noinst_HEADERS += $(JFLIB_I)/multiplexed_io.hpp $(JFLIB_I)/pool.hpp	\
                  $(JFLIB_I)/multiplexed_parser.hpp

##############################
# utils lib
##############################
UTILS_I = include
noinst_HEADERS += $(UTILS_I)/charb.hpp $(UTILS_I)/charbuf.hpp		\
                  $(UTILS_I)/exp_buffer.hpp				\
                  $(UTILS_I)/gzip_stream.hpp $(UTILS_I)/misc.hpp	\
                  $(UTILS_I)/heap.hpp $(UTILS_I)/reallocators.hpp	\
                  $(UTILS_I)/skip_list_set.hpp				\
                  $(UTILS_I)/skip_list_map.hpp				\
                  $(UTILS_I)/gcc_builtins.hpp $(UTILS_I)/rb_tree.hpp	\
                  $(UTILS_I)/rb_node.hpp				\
                  $(UTILS_I)/binary_search_tree_iterator.hpp		\
                  $(UTILS_I)/skip_list_common.hpp			\
                  $(UTILS_I)/multi_thread_skip_list_set.hpp		\
                  $(UTILS_I)/multi_thread_skip_list_map.hpp		\
                  $(UTILS_I)/tmpstream.hpp $(UTILS_I)/exp_vector.hpp	\
                  $(UTILS_I)/thread_pool.hpp $(UTILS_I)/fork_exec.hpp $(UTILS_I)/timer.h

noinst_LIBRARIES = libutils.a
libutils_a_SOURCES = lib/fork_exec.cpp
# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
YAGGO_SOURCES = # Append all file to be built by yaggo

##############################
# src pipeline
##############################
bin_PROGRAMS += reverse_complement assign_reads reduce_sr	\
                sorted_merge create_k_unitigs_large_k 			\
                create_k_unitigs_large_k2 homo_trim bloom_query		\
                rename_filter_fastq semaphore

dist_bin_SCRIPTS += src/getNumBasesPerReadInFastaFile.perl		\
                    src/extractreads.pl src/extract_unjoined_pairs.pl	\
                    src/extractreads_not.pl  				\
                    src/expand_fastq src/createSuperReads_RNA

noinst_HEADERS += src/aligned_simple_array.hpp				\
                  src/diskBasedUnitigger.h src/fibonacci_coding.hpp	\
                  src/sr_names.hpp src/bloom_filter.hpp			\
                  src/MurmurHash3.h src/bloom_hash.hpp			\
                  src/bloom_counter2.hpp src/read_parser.hpp		\
                  src/create_k_unitigs_common.hpp src/superread_index.hpp

assign_reads_SOURCES = src/assign_reads.cpp src/superread_index.cpp
homo_trim_SOURCES = src/homo_trim.cc
reverse_complement_SOURCES = src/reverse_complement.cc
reduce_sr_SOURCES = src/reduce_sr.cc
YAGGO_SOURCES    += src/reduce_sr_cmdline.hpp
sorted_merge_SOURCES = src/sorted_merge.cc
YAGGO_SOURCES       += src/sorted_merge_cmdline.hpp
create_k_unitigs_large_k_SOURCES = src/create_k_unitigs_large_k.cc
YAGGO_SOURCES += src/create_k_unitigs_large_k_cmdline.hpp
create_k_unitigs_large_k2_SOURCES = src/create_k_unitigs_large_k2.cc
YAGGO_SOURCES += src/create_k_unitigs_large_k2_cmdline.hpp
bloom_query_SOURCES = src/bloom_query.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src/bloom_query_cmdline.hpp
rename_filter_fastq_SOURCES = src/rename_filter_fastq.cc
YAGGO_SOURCES += src/rename_filter_fastq_cmdline.hpp

########################
# masurca script files #
########################
bin_SCRIPTS += src/MasurcaPaths.pm
CLEANFILES += src/MasurcaPaths.pm
EXTRA_DIST += src/MasurcaPaths.pm.in
dist_bin_SCRIPTS += src/MasurcaConf.pm		\
                    src/MasurcaCommon.pm 	\
                    src/MasurcaSuperReads.pm

src/MasurcaPaths.pm: src/MasurcaPaths.pm.in Makefile
	mkdir -p $(dir $@)
	sed -e 's,[@]bindir[@],$(bindir),g' \
	    -e 's,[@]libdir[@],$(libdir),g' \
	    -e 's,[@]PERL_EXT_DIR[@],$(PERL_EXT_LIB),g' < $< > $@
	chmod +x $@


##############################
# src2 pipeline
##############################
bin_PROGRAMS += createFastaSuperReadSequences				\
                findMatchesBetweenKUnitigsAndReads			\
                createKUnitigMaxOverlaps				\
                eliminateBadSuperReadsUsingList				\
                getSuperReadInsertCountsFromReadPlacementFile		\
                getSuperReadInsertCountsFromReadPlacementFileTwoPasses	\
                joinKUnitigs_v3 outputRecordsNotOnList \
		create_sr_frg splitFileAtNs

dist_bin_SCRIPTS += src2/getLengthStatisticsForKUnitigsFile.perl	\
                    src2/createSuperReadsForDirectory.perl src2/add_missing_mates.pl \
                    src2/outputAlekseysJellyfishReductionFile.perl \
                    src2/translateReduceFile.perl src/mergeSuperReadsUniquely.pl 

createFastaSuperReadSequences_SOURCES = src2/createFastaSuperReadSequences.cc

findMatchesBetweenKUnitigsAndReads_SOURCES = src/read_parser.cc				\
                                             src2/findMatchesBetweenKUnitigsAndReads.cc	\
                                             src/gcc_builtins.cc
create_sr_frg_SOURCES = src2/create_sr_frg.cc
YAGGO_SOURCES += src2/findMatchesBetweenKUnitigsAndReads_cmdline.hpp

createKUnitigMaxOverlaps_SOURCES = src2/createKUnitigMaxOverlaps.cc
eliminateBadSuperReadsUsingList_SOURCES = src2/eliminateBadSuperReadsUsingList.cc
getSuperReadInsertCountsFromReadPlacementFile_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFile.cc	\
src/fibonacci_coding.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFile_cmdline.hpp
getSuperReadInsertCountsFromReadPlacementFileTwoPasses_SOURCES =	\
src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses.cc		\
src/fibonacci_coding.cc src/MurmurHash3.cpp
YAGGO_SOURCES += src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses_cmdline.hpp
joinKUnitigs_v3_SOURCES = src2/joinKUnitigs_v3.cc
YAGGO_SOURCES          += src2/joinKUnitigs_v3_cmdline.hpp
outputRecordsNotOnList_SOURCES = src2/outputRecordsNotOnList.cc
splitFileAtNs_SOURCES = src2/splitFileAtNs.cc

semaphore_SOURCES = src/semaphore.cc
YAGGO_SOURCES += src/semaphore_cmdline.hpp

##############################
# unittests
##############################
TESTS = all_tests
EXTRA_DIST += $(TESTS)
check_PROGRAMS = all_tests

all_tests_SOURCES = unittests/test_main.cc unittests/test_charb.cc	\
                    unittests/test_exp_buffer.cc			\
                    unittests/test_charbuf.cc unittests/test_heap.cc	\
                    unittests/test_sr_names.cc				\
                    unittests/test_bloom_filter.cc			\
                    unittests/test_bloom_counter2.cc			\
                    unittests/test_skip_list_set.cc			\
                    unittests/test_skip_list_map.cc			\
                    unittests/test_rbtree.cc				\
                    unittests/test_multi_thread_skip_list_set.cc	\
                    unittests/test_multi_thread_skip_list_map.cc	\
                    unittests/test_read_parser.cc			\
                    unittests/test_hash.cc				\
                    unittests/test_exp_vector.cc			\
                    unittests/test_thread_pool.cc			\
                    unittests/test_fork_exec.cc
all_tests_SOURCES += src/fibonacci_coding.cc src/MurmurHash3.cpp	\
                     src/gcc_builtins.cc src/read_parser.cc		\
                     unittests/misc.cc
all_tests_CXXFLAGS = $(AM_CXXFLAGS) $(JELLYFISH2_0_CFLAGS)	\
-I$(srcdir)/unittests/gtest/include				\
-I$(srcdir)unittests/gtest -I$(srcdir)/include		\
-I$(srcdir)
all_tests_LDFLAGS = $(AM_LDFLAGS) $(JELLYFISH2_0_LIBS)
all_tests_LDADD = libutils.a libgtest.a
YAGGO_SOURCES += unittests/test_main_cmdline.hpp
noinst_HEADERS += unittests/test_main_cmdline.hpp unittests/misc.hpp

include gtest.mk

##############################
# Doxygen
##############################
#include doxygen.am
# MOSTLYCLEANFILE = doxygen-doc
# EXTRA_DIST += doxygen.cfg
# .PHONY: doxygen-doc
# doxygen-doc: doxygen.cfg
# 	rm -rf $@
# 	doxygen $<

##############################
# Extra configuration for development tree
##############################
-include $(srcdir)/development.mk
-include $(srcdir)/personal.mk

##############################
# Info
##############################
print-%:
	@echo $($*)

##############################
# data generation
##############################
#dist_data_DATA = data/adapter_0
#CLEANFILES += data/adapter_0
#SUFFIXES = _0 .fa
#.fa_0:
#	$(AM_V_GEN)mkdir -p $(@D); jellyfish count -m 24 -s 5k -C -t 10 -r -o $$(dirname $@)/$$(basename $@ _0) $<

