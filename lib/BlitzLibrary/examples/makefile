
# Path where Blitz is installed
# IMPORTANT: you may have to change this
BZDIR = /site/Blitz++-0.9_64

CXX = g++

# Flags for optimizing executables
# CXXFLAGS = -02 -I$(BZDIR) -ftemplate-depth-30

# Flags for debugging
CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG -I$(BZDIR)/include

LDFLAGS = 

LIBS = -L$(BZDIR)/lib -lblitz -lm

TARGETS = blitz_test

.SUFFIXES: .o.cpp

.cpp.o:
		$(CXX) $(CXXFLAGS) -c $*.cpp

$(TARGETS):     
		$(CXX)  $(LDFLAGS)  $@.o -o $@ $(LIBS)


all:        
		$(TARGETS)

blitz_test:		blitz_test.o


clean: 
		rm -f *.o 
