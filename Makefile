# Try "make CC=icc ARCH=[your usual ARCH value]"

CPPFLAGS = -I$(TREEHOME)/include -wd1338,810
CFLAGS = -g -O0 -Wall -D_FILE_OFFSET_BITS=64
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

#PROGS = zonetopartid
#PROGS= readsection3
#PROGS= addabundance
#PROGS= mergeSDFs2
#PROGS= SDFtoASCII-batch
#PROGS= SDFtoASCII
#PROGS= maketraj
#PROGS= readanSDF
#PROGS= convertff
#PROGS= get_bndry
PROGS= getyield

.PHONY: all clean

all: $(PROGS)

#zonetopartid: zonetopartid.o
#readsection3: readsection3.o
#addabundance: addabundance.o
#mergeSDFs2: mergeSDFs2.o
#SDFtoASCII-batch: SDFtoASCII-batch.o
#SDFtoASCII: SDFtoASCII.o
#maketraj: maketraj.o
#readanSDF: readanSDF.o
#convertff: convertff.o
#get_bndry: get_bndry.o
getyield: getyield.o

clean:
	-$(RM) $(PROGS) *.o *~
