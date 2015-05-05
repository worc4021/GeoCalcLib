CFLAGS = -Wall -fpic
# DFLAG = -g
CFLAGS += -DTIMES -DSIGNALS -DB64 -DLRS_QUIET -DNOINFO

# Set DEBUG macro for functions to print GMPmat and mpx_t type data in debugger, 
# set NOINFO macro to avoid printing out number of rows/vertices/rays found during computation.

# CFLAGS += -DDEBUG

# Path of GMP shared library, -lgmp does not work! It has to be built with an 64 bit ABI, 
# i.e. build with ./configure ABI=64 !
GMPdir = /Users/Manuel/Documents/Development/GMP
GMPlib = $(GMPdir)/lib/libgmp.10.dylib
GMPinc = $(GMPdir)/include/


# Linker flags, do not modify!
LFLAGS = -shared -Wl,-no_pie $(GMPlib) -lmx -lmex -lmat

MATLABROOT = /Applications/MATLAB_R2015a.app
MATLABINCLUDEDIR = $(MATLABROOT)/extern/include/
MATLABLIB = $(MATLABROOT)/bin/maci64/

# Path to which everything should be installed, has to be on Matlab path!
INSTALLDIR = /Users/Manuel/Documents/MATLAB/Funktionen/

OBJECTS = mainFunctions.o translation_functions.o lrslib.o lrsgmp.o

EXTENTION = dylib


all: libgeocalc.$(EXTENTION)
	mv libgeocalc.$(EXTENTION) $(INSTALLDIR)libgeocalc.$(EXTENTION)
	cp *.m $(INSTALLDIR)
	rm *.o


libgeocalc.$(EXTENTION): $(OBJECTS)
	$(CC) $(LFLAGS) -L$(MATLABLIB) $^ -o libgeocalc.$(EXTENTION) 

.c.o:
	$(CC) $(DFLAG) $(CFLAGS) -I$(MATLABINCLUDEDIR) -I$(GMPinc) $< -o $@ -c

clean:
	rm -f *.o *.$(EXTENTION)