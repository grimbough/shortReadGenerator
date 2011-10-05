CPP=g++
LIBRARIES = -lboost_program_options
CPPFLAGS = -g -O2 -Wall 
OBJ=shortReadGenerator.o genome.o reads.o

all: $(OBJ)
	$(CPP) $(CPPFLAGS) $(LIBRARIES) -o shortReadGenerator $(OBJ)

clean:
	rm -rf $(OBJ) readLengths