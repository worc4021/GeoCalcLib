UNAME_S := $(shell uname -s)

MATLABROOT = /Applications/MATLAB_R2016a.app
MATLABINCLUDEDIR = $(MATLABROOT)/extern/include/

# Modify for your distribution:
MATLABLIB = $(MATLABROOT)/bin/maci64/


CFLAGS = -Wall -fpic
CFLAGS += -DGMP -DTIMES -DSIGNALS -DB64 -DLRS_QUIET -DNOINFO


# Set DEBUG macro for functions to print GMPmat and mpx_t type data in debugger, 
# set NOINFO macro to avoid printing out number of rows/vertices/rays found during computation.
# CFLAGS += -DDEBUG
DFLAG = -g


# Linker flags, do not modify!
LFLAGS = -shared -lmx -lmex -lmat -lgmp -L$(MATLABLIB)

ifeq ($(UNAME_S),Darwin)
	LFLAGS += -Wl,-no_pie
	EXTENTION = dylib
endif
ifeq ($(UNAME_S),Linux)
	LFLAGS += -Wl,-rpath,$(MATLABLIB)
	EXTENTION = so
endif

# Path to which everything should be installed, has to be on Matlab path!
INSTALLDIR = /Users/Manuel/Documents/MATLAB/Funktionen/

OBJECTS = mainFunctions.o translation_functions.o lrslib.o lrsgmp.o
TESTOBJECTS = testcase.o translation_functions.o lrslib.o lrsgmp.o


all: libgeocalc.$(EXTENTION)
	mv libgeocalc.$(EXTENTION) $(INSTALLDIR)libgeocalc.$(EXTENTION)
	cp *.m $(INSTALLDIR)
	rm *.o


libgeocalc.$(EXTENTION): $(OBJECTS)
	$(CC) $(LFLAGS) $(DFLAG) $^ -o libgeocalc.$(EXTENTION) 

.c.o:
	$(CC) $(DFLAG) $(CFLAGS) -I$(MATLABINCLUDEDIR) $< -o $@ -c

clean:
	rm -f *.o *.$(EXTENTION) testcase

testcase: CFLAGS += -DNOMATLAB -DDEBUG

testcase: $(TESTOBJECTS)
	$(CC) -lgmp -lm $(DFLAG) $^ -o testcase

