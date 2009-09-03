# Try "make ARCH=[your usual ARCH value]"

CPPFLAGS = -I$(TREEHOME)/include
CFLAGS = -g -O2 -Wall
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

PROGS = todoublepos

.PHONY: all clean

all: $(PROGS)

todoublepos: todoublepos.o

clean:
	-$(RM) $(PROGS) *.o *~
