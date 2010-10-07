# Try "make CC=icc ARCH=[your usual ARCH value]"

CPPFLAGS = -I$(TREEHOME)/include -wd1338,810
CFLAGS = -g -O0 -Wall -D_FILE_OFFSET_BITS=64
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

#PROGS = zonetopartid
#PROGS= readsection
#PROGS= addabundance
#PROGS= mergeSDFs
#PROGS= SDFtoASCII-batch
#PROGS= SDFtoASCII
#PROGS= maketraj
PROGS= readanSDF

.PHONY: all clean

all: $(PROGS)

#zonetopartid: zonetopartid.o
#readsection: readsection.o
#addabundance: addabundance.o
#mergeSDFs: mergeSDFs.o
#SDFtoASCII-batch: SDFtoASCII-batch.o
#SDFtoASCII: SDFtoASCII.o
#maketraj: maketraj.o
readanSDF: readanSDF.o

clean:
	-$(RM) $(PROGS) *.o *~
