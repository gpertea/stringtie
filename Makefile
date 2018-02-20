BAM := ./samtools-0.1.18
#path to the directory where the samtools package was built (in place)
#so libbam.a and *.h files MUST be in here

GDIR      :=./gclib

INCDIRS   := -I. -I${GDIR} -I${BAM}

CC        := g++

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-strict-aliasing -fno-exceptions -fno-rtti

LINKER    := g++

LDFLAGS   := -g -L${BAM} $(LDFLAGS)

LIBS      := -lbam -lz

ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
 NOTHREADS=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CC} -dumpmachine)))
 WINDOWS=1
endif

# MinGW32 GCC 4.5 link problem fix
#ifdef WINDOWS
ifneq (,$(findstring 4.5.,$(shell g++ -dumpversion)))
 STATIC_CLIB=1
endif
#endif

# Misc. system commands
ifdef WINDOWS
 RM = del /Q
else
 RM = rm -f
endif

# File endings
ifdef WINDOWS
 EXE = .exe
else
 EXE =
endif

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS += -lpthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

#ifneq (,$(findstring release,$(MAKECMDGOALS)))

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CFLAGS := -DNDEBUG -g $(BASEFLAGS) $(CFLAGS) -O3
else
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     MEMCHECK_BUILD=1
     GCCVER49 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CFLAGS := -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
     GCCVER5 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector $(CFLAGS)
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
    ifeq (,$(findstring clean,$(MAKECMDGOALS)))
     #just plain debug build
     DEBUG_BUILD=1
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
    endif
  endif
endif

ifdef RELEASE_BUILD
 ifneq (,$(findstring static, $(MAKECMDGOALS)))
    STATIC_CLIB=1
 endif
endif

ifdef STATIC_CLIB
 LDFLAGS += -static-libstdc++ -static-libgcc
endif

ifdef DEBUG_BUILD
  $(warning Building DEBUG version [much slower], use 'make release' for a fast, optimized program.)
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make release" for a fast, optimized program.'
endif


OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 


ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif


%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

OBJS += rlink.o tablemaker.o tmerge.o

all release static debug: stringtie${EXE}
memcheck memdebug: stringtie${EXE}
memuse memusage memtrace: stringtie${EXE}
nothreads: stringtie${EXE}

${GDIR}/GBam.o : $(GDIR)/GBam.h
stringtie.o : $(GDIR)/GBitVec.h $(GDIR)/GHash.hh $(GDIR)/GBam.h
rlink.o : rlink.h tablemaker.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
tmerge.o : rlink.h tmerge.h
tablemaker.o : tablemaker.h rlink.h
${BAM}/libbam.a: 
	cd ${BAM} && make lib
stringtie: ${BAM}/libbam.a $(OBJS) stringtie.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}

.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

clean:
	${RM} stringtie${EXE} stringtie.o*  $(OBJS)
	${RM} core.*
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} stringtie${EXE} stringtie.o* $(OBJS)
	${RM} core.*
