ARCH = g5

CC = gcc
CPPFLAGS = -I$(TREEHOME)/include
CFLAGS = -g -Wall
LDFLAGS = -L$(TREEHOME)/Objfiles/$(ARCH)
LDLIBS = -lsw

PROGS = todoublepos

.PHONY: all clean

all: $(PROGS)

clean:
	-$(RM) $(PROGS) *.o *~
