CFLAGS = -Wall -fpic
DFLAG = -g
CFLAGS += -DTIMES -DSIGNALS -DB64 -DLRS_QUIET
# Set DEBUG macro for functions to print GMPmat and mpx_t type data in debugger, 
# set NOINFO macro to avoid printing out number of rows/vertices/rays found during computation.
CFLAGS += -DDEBUG -DNOINFO
LFLAGS = -shared -Wl,-no_pie $(GMP) -lmx -lmex -lmat

# Path of GMP shared library, -lgmp does not work! It has to be built with an 64 bit ABI, 
# i.e. build with ./configure ABI=64 !
GMP = /Users/Manuel/Documents/Development/GMPFiles/MexExp/lib/libgmp.10.dylib

MATLABINCLUDEDIR = /Applications/MATLAB_R2015a.app/extern/include/
MATLABLIB = /Applications/MATLAB_R2015a.app/bin/maci64/

# Path to which everything should be installed, has to be on Matlab path!
INSTALLDIR = /Users/Manuel/Documents/MATLAB/Funktionen/

OBJECTS = mainFunctions.o translation_functions.o lrslib.o lrsgmp.o

all: libgeocalc.dylib
	mv libgeocalc.dylib $(INSTALLDIR)libgeocalc.dylib
	cp *.m $(INSTALLDIR)
	rm *.o


libgeocalc.dylib: $(OBJECTS)
	$(CC) $(LFLAGS) -L$(MATLABLIB) $^ -o libgeocalc.dylib 

.c.o:
	$(CC) $(DFLAG) $(CFLAGS) -I$(MATLABINCLUDEDIR) $< -o $@ -c

clean:
	rm -f *.o *.dylib