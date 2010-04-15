# Try "make CC=icc ARCH=[your usual ARCH value]"

CPPFLAGS = -I$(TREEHOME)/include -wd1338,810
CFLAGS = -g -O0 -Wall
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

#PROGS = todoublepos2
#PROGS= readsection
#PROGS= addabundance
#PROGS= mergeSDFs
PROGS= SDFtoASCII-batch
#PROGS= maketraj

.PHONY: all clean

all: $(PROGS)

#todoublepos2: todoublepos2.o
#readsection: readsection.o
#addabundance: addabundance.o
#mergeSDFs: mergeSDFs.o
SDFtoASCII-batch: SDFtoASCII-batch.o
#maketraj: maketraj.o

clean:
	-$(RM) $(PROGS) *.o *~
