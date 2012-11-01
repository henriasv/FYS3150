CC 			= g++
MYPROGRAM 	= project3 
LDFLAGS 	= -larmadillo -lblas -llapack -openmp
SOURCES 	= main.cpp Helium_Solver.cpp
OBJECTS		= main.o Helium_Solver.o

#$(MYPROGRAM): $(CLASSES).o

$(MYPROGRAM): $(OBJECTS)

	$(CC) -o $(MYPROGRAM) $(OBJECTS) $(LDFLAGS)

Helium_Solver.o: Helium_Solver.h
main.o: Helium_Solver.h

$(OBJECTS): $(SOURCES)
	$(CC) -c $(SOURCES)

clean:

	rm *.o $(MYPROGRAM)
	
