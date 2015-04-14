CFLAGS = -Wall -fpic
DFLAG = -g
CFLAGS += -DTIMES -DSIGNALS -DB64 -DLRS_QUIET -DDEBUG
LFLAGS = -shared -Wl,-no_pie $(GMP) -lmx -lmex -lmat

# Path of GMP shared library, -lgmp does not work!
GMP = /Users/Manuel/Documents/Development/GMPFiles/MexExp/lib/libgmp.10.dylib
# GMP = /opt/local/lib/engines/libgmp.dylib

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