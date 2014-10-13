BAM := ./samtools-0.1.18
#path to the directory where the samtools package was built (in place)
#so libbam.a and *.h files MUST be in here

GDIR := ./gclib

SEARCHDIRS := -I. -I${GDIR} -I${BAM}

#SYSTYPE :=     $(shell uname)

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


BASEFLAGS  := -Wall -Wextra ${SEARCHDIRS} $(MARCH) -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-strict-aliasing -fno-exceptions -fno-rtti

# C/C++ linker

LINKER  := g++

LIBS := -lbam -lz

# Non-windows systems need pthread
ifndef WINDOWS
 #ifndef NOTHREADS
 # unfortunately the new samtools requires pthread anyway
   LIBS += -lpthread
 #endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

###----- generic build rule

ifneq (,$(findstring release,$(MAKECMDGOALS)))
  CFLAGS := -O2 -DNDEBUG -g $(BASEFLAGS)
  LDFLAGS := -g -L${BAM} ${LFLAGS}
else
  CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  LDFLAGS := -g -L${BAM}
  GDEBUG = 1
endif

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

ifdef GDEBUG
 OBJS += ${GDIR}/proc_mem.o
endif

OBJS += rlink.o
 
.PHONY : all debug clean release nothreads
all:     stringtie
release: stringtie
debug:   stringtie
nothreads: stringtie

${GDIR}/GBam.o : $(GDIR)/GBam.h
stringtie.o : $(GDIR)/GBitVec.h $(GDIR)/GHash.hh $(GDIR)/GBam.h
rlink.o : rlink.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h
${BAM}/libbam.a: 
	cd ${BAM} && make lib
stringtie: ${BAM}/libbam.a $(OBJS) stringtie.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

# target for removing all object files

clean:
	cd ${BAM} && make clean
	@${RM} stringtie stringtie.o* stringtie.exe $(OBJS)
	@${RM} core.*


