objects = main.o kbd.o command.o display.o \
	h_sb.o h_s.o h_b.o bath_terms.o \
	initial_state.o 

objects = main.o utils.o h_sb.o h_res.o

options = -ggdb -fopenmp

spin_boson : $(objects)
	g++ $(options) -o spin_boson $(objects)

main.o : main.cpp utils.hpp
	g++ $(options) -c -o main.o main.cpp

h_res.o : utils.hpp h_funcs.hpp h_res.cpp
	g++ $(options) -c -o h_res.o h_res.cpp
h_sb.o : utils.hpp h_funcs.hpp h_sb.cpp	
	g++ $(options) -c -o h_sb.o h_sb.cpp

utils.o : utils.cpp utils.hpp 
	g++ $(options) -c -o utils.o utils.cpp

.PHONY : clean

clean :
	-rm spin_boson $(objects)
