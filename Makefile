CPP=g++
FLAGS=-g -lgsl -lm

OUT_EXEC = main
objects= main.o

## Headers, and their dependencies:
##
matrices_hpp = matrices.hpp
## Exectuable:
##
$(OUT_EXEC) : $(objects)
	$(CPP) $^ -o $@ $(FLAGS)

## Object Files
##
main.o: main.cpp $(matrices_hpp)
	$(CPP) -c $< $(FLAGS)

## Cleanup
##
clean :
	rm -f $(objects) $(OUT_EXEC)
neat :
	rm -f $(objects)
