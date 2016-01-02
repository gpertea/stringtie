BAM := ./samtools-0.1.18
#path to the directory where the samtools package was built (in place)
#so libbam.a and *.h files MUST be in here

GDIR :=./gclib

SEARCHDIRS := -I. -I${GDIR} -I${BAM}

#CC := clang++
CC      := g++


ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
 NOTHREADS=1
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CC} -dumpmachine)))
 WINDOWS=1
endif

LFLAGS = 
# MinGW32 GCC 4.5 link problem fix
#ifdef WINDOWS
ifneq (,$(findstring 4.5.,$(shell g++ -dumpversion)))
 LFLAGS += -static-libstdc++ -static-libgcc
endif
#endif

# Misc. system commands
#ifdef WINDOWS
#RM = del /Q
#else
RM = rm -f
#endif

# File endings
ifdef WINDOWS
EXE = .exe
else
EXE =
endif


BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-strict-aliasing -fno-exceptions -fno-rtti -D_DARWIN_C_SOURCE

# C/C++ linker

#LINKER := clang++
LINKER  := g++

LIBS := -lbam -lz

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS += -lpthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

###----- generic build rule

ifneq (,$(findstring release,$(MAKECMDGOALS)))
  CFLAGS := -O3 -DNDEBUG -g $(BASEFLAGS)
  LDFLAGS := -g -L${BAM} ${LFLAGS}
else
  #make memcheck : use the statically linked address sanitizer in gcc 4.9.x
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     CFLAGS := -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector $(CFLAGS)
     LDFLAGS := -g -L${BAM}
     #LIBS := -Wl,-Bstatic -lasan -lubsan -Wl,-Bdynamic -ldl $(LIBS)
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
   ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
       BASEFLAGS += -DGMEMTRACE
       GMEMTRACE=1
   endif
   #just plain debug build
    CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
    LDFLAGS := -g -L${BAM}
  endif
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifdef GMEMTRACE
 OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif



OBJS += rlink.o tablemaker.o tmerge.o
 
.PHONY : all debug clean cleanall cleanAll allclean release nothreads
all:     stringtie
release: stringtie
debug:   stringtie
memcheck: stringtie
memdebug: stringtie
memusage:  stringtie
memuse:    stringtie
memtrace:  stringtie
nothreads: stringtie

${GDIR}/GBam.o : $(GDIR)/GBam.h
stringtie.o : $(GDIR)/GBitVec.h $(GDIR)/GHash.hh $(GDIR)/GBam.h
rlink.o : rlink.h tablemaker.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
tmerge.o : rlink.h tmerge.h
tablemaker.o : tablemaker.h rlink.h
${BAM}/libbam.a: 
	cd ${BAM} && make lib
stringtie: ${BAM}/libbam.a $(OBJS) stringtie.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

clean:
	@${RM} stringtie stringtie.o* stringtie.exe $(OBJS)
	@${RM} core.*
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	@${RM} stringtie stringtie.o* stringtie.exe $(OBJS)
	@${RM} core.*


