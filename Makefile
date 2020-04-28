CC = g++

CXXFLAGS=-std=c++11 -O3 -Wall 

pcgSolver: main.o mult.o
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS) -lm

.c.o:
	$(CC) $(CXXFLAGS) $(INCLUDE) -c $< -lm

.PHONY: clean
clean:
	rm -f *.o *~ pcgSolver 
