CPP=g++
FLAGS=-g -lgsl -lm

OUT_EXEC = main
objects= functions.o main2.o

## Headers, and their dependencies:
##
matrices_hpp = matrices.hpp
functions_hpp = functions.hpp $(matrices_hpp)
global_hpp = global.hpp
## Exectuable:
##
$(OUT_EXEC) : $(objects)
	$(CPP) $^ -o $@ $(FLAGS)

## Object Files
##
functions.o: functions.cpp $(functions_hpp) $(global_hpp)
	$(CPP) -c $< $(FLAGS)

main2.o: main2.cpp $(functions_hpp) $(global_hpp)
	$(CPP) -c $< $(FLAGS)

## Cleanup
##
clean :
	rm -f $(objects) $(OUT_EXEC)
neat :
	rm -f $(objects)
