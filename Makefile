UNAME_S := $(shell uname -s)
CURRDIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
include $(CURRDIR)/User.make

MATLABINCLUDEDIR = $(MATLABROOT)/extern/include/


CFLAGS = -Wall -fpic
CFLAGS += -DGMP -DTIMES -DSIGNALS -DB64 -DLRS_QUIET -DNOINFO -DDIRECT


# Set DEBUG macro for functions to print GMPmat and mpx_t type data in debugger, 
# set NOINFO macro to avoid printing out number of rows/vertices/rays found during computation.
# CFLAGS += -DDEBUG
# DFLAG = -g


# Linker flags, do not modify!
LFLAGS = -shared -lmx -lmex -lmat -lgmp -L$(MATLABLIB)

ifeq ($(UNAME_S),Darwin)
	LFLAGS += -Wl,-no_pie
	EXTENTION = dylib
	MATLABLIB = $(MATLABROOT)/bin/maci64/
endif
ifeq ($(UNAME_S),Linux)	
	LFLAGS += -Wl,-rpath,$(MATLABLIB)
	EXTENTION = so
	MATLABLIB = $(MATLABROOT)/bin/glnxa64/
endif

OBJECTS = mainFunctions.o translation_functions.o lrslib.o lrsgmp.o
CC = gcc

all: libgeocalc.$(EXTENTION)
	mv libgeocalc.$(EXTENTION) $(INSTALLDIR)/libgeocalc.$(EXTENTION)
	cp *.m $(INSTALLDIR)
	rm *.o

libgeocalc.$(EXTENTION): $(OBJECTS)
	$(CC) $(LFLAGS) $(DFLAG) $^ -o libgeocalc.$(EXTENTION) 

.c.o:
	$(CC) $(DFLAG) $(CFLAGS) -I$(MATLABINCLUDEDIR) $< -o $@ -c

clean:
	rm -f *.o *.$(EXTENTION) testcase