CPP=g++
FLAGS=-g -lgsl -lm

OUT_EXEC = main2
objects= functions.o main2.o

## Headers, and their dependencies:
##
matrices_hpp = matrices.hpp

## Exectuable:
##
$(OUT_EXEC) : $(objects)
	$(CPP) $^ -o $@ $(FLAGS)

## Object Files
##
functions.o: functions.cpp $(matrices_hpp)
	$(CPP) -c $< $(FLAGS)

main2.o: main2.cpp $(matrices_hpp)
	$(CPP) -c $< $(FLAGS)

## Cleanup
##
clean :
	rm -f $(objects) $(OUT_EXEC)
neat :
	rm -f $(objects)
