CPP=g++
FLAGS=

OUT_EXEC = main
objects = main.o 

## Headers, and their dependencies:
##

## Exectuable:
##
$(OUT_EXEC) : $(objects)
	$(CPP) $^ -o $@ $(FLAGS)

## Object Files
##
main.o: main.cpp
	$(CPP) -c $< $(FLAGS)

## Cleanup
##
clean :
	rm -f $(objects) $(OUT_EXEC)
neat :
	rm -f $(objects)
