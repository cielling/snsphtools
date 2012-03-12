# Try "make CC=icc ARCH=[your usual ARCH value] PROGS=[program name (w/o '.c)]"
# ARCH:
# on Lonestar: tacc-gcc, no CC flag
# on Saguaro: c2icc
# on Mapache/Conejo: cjicc or mp-icc

CPPFLAGS = -I$(TREEHOME)/include -wd1338,810
CFLAGS = -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw -lm

.PHONY: all clean

all: $(PROGS)

$(PROGS): $(PROGS).o nrutil.o

nrutil.o: nrutil.c 

clean:
	-$(RM) $(PROGS) *.o *~
