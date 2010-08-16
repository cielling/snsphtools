# Try "make CC=icc ARCH=[your usual ARCH value]"

CPPFLAGS = -I$(TREEHOME)/include -wd1338,810
CFLAGS = -g -O0 -Wall
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

#PROGS = zonetopartid
#PROGS= readsection
#PROGS= addabundance
#PROGS= mergeSDFs
PROGS= SDFtoASCII-batch
#PROGS= SDFtoASCII
#PROGS= maketraj

.PHONY: all clean

all: $(PROGS)

#zonetopartid: zonetopartid.o
#readsection: readsection.o
#addabundance: addabundance.o
#mergeSDFs: mergeSDFs.o
SDFtoASCII-batch: SDFtoASCII-batch.o
#SDFtoASCII: SDFtoASCII.o
#maketraj: maketraj.o

clean:
	-$(RM) $(PROGS) *.o *~
