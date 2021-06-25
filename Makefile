#-- for now these MUST point to the included "samtools-0.x.x" and "gclib" sub-directories
HTSLIB  := ./htslib
#-- 
LIBDEFLATE := ${HTSLIB}/xlibs/lib/libdeflate.a
LIBBZ2 := ${HTSLIB}/xlibs/lib/libbz2.a
LIBLZMA := ${HTSLIB}/xlibs/lib/liblzma.a

GDIR := ./gclib
#--

INCDIRS := -I. -I${GDIR} -I${HTSLIB}

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++11 -fno-strict-aliasing -fno-exceptions -fno-rtti
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

# LDFLAGS += -L${BAM}

LIBS    := ${HTSLIB}/libhts.a ${LIBBZ2} ${LIBLZMA} ${LIBDEFLATE} -lz -lm

ifneq (,$(filter %nothreads %prof %profile, $(MAKECMDGOALS)))
 NOTHREADS=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
#ifdef WINDOWS ##<-- use MSYS
# RM = del /Q
#else
RM = rm -f
#endif

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS := -pthread ${LIBS}
   BASEFLAGS += -pthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

# Compiling for Windows with MinGW?
#ifneq ($(findstring -mingw,$(shell $(CC) -dumpmachine 2>/dev/null)),)
#LIBS += -lregex -lws2_32
#endif
# File endings
ifdef WINDOWS
 EXE = .exe
 LIBS += -lregex -lws2_32
else
 EXE =
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static %static-cpp, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     ifneq (,$(filter %prof %profile, $(MAKECMDGOALS)))
     ## profiling build
       CXXFLAGS := -DNDEBUG $(BASEFLAGS) -g -pg
       LDFLAGS += -g -pg
     else
        #just plain debug build
        DEBUG_BUILD=1
        CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
        ifneq (, $(findstring darwin, $(DMACH)))
           CXXFLAGS += -gdwarf-3
        endif
        CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
     endif
  endif
endif

ifdef RELEASE_BUILD
 ifneq ($(findstring static,$(MAKECMDGOALS)),) 
  # static or static-cpp found
  ifneq ($(findstring static-cpp,$(MAKECMDGOALS)),) 
     #not a full static build, only c/c++ libs
     LDFLAGS := -static-libgcc -static-libstdc++ ${LDFLAGS}
  else
     #full static build
     LDFLAGS := -static -static-libgcc -static-libstdc++ ${LDFLAGS}
  endif
 endif
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of stringtie.. )
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GSam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

OBJS += rlink.o tablemaker.o tmerge.o

all release static static-cpp debug: stringtie${EXE}
memcheck memdebug tsan tcheck thrcheck: stringtie${EXE}
memuse memusage memtrace: stringtie${EXE}
prof profile: stringtie${EXE}
nothreads: stringtie${EXE}

stringtie.o : $(GDIR)/GBitVec.h $(GDIR)/GHashMap.hh $(GDIR)/GSam.h
rlink.o : rlink.h tablemaker.h $(GDIR)/GSam.h $(GDIR)/GBitVec.h
tmerge.o : rlink.h tmerge.h
tablemaker.o : tablemaker.h rlink.h

##${BAM}/libbam.a: 
##	cd ${BAM} && make lib

${HTSLIB}/libhts.a:
	cd ${HTSLIB} && ./build_lib.sh

stringtie${EXE}: ${HTSLIB}/libhts.a $(OBJS) stringtie.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
test demo tests: stringtie${EXE}
	@./run_tests.sh
.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

#	echo $(PATH)
clean:
	${RM} stringtie${EXE} stringtie.o*  $(OBJS)
	${RM} core.*
##allclean cleanAll cleanall:
##	cd ${BAM} && make clean
##	${RM} stringtie${EXE} stringtie.o* $(OBJS)
##	${RM} core.*
