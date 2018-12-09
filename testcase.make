UNAME_S := $(shell uname -s)
# CURRDIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
# include $(CURRDIR)/User.make

# MATLABINCLUDEDIR = -I$(MATLABROOT)/extern/include/


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

TESTOBJECTS = testcase.o translation_functions.o lrslib.o lrsgmp.o

.c.o:
	$(CC) $(DFLAG) $(CFLAGS) $(MATLABINCLUDEDIR) $< -o $@ -c

clean:
	rm -f *.o *.$(EXTENTION) testcase

testcase: CFLAGS += -DNOMATLAB -DDEBUG -DALTERNATIVEREDUCTION

testcase: $(TESTOBJECTS)
	$(CC) -lgmp -lm $(DFLAG) $^ -o testcase

