objects = main.o kbd.o command.o display.o \
	h_sb.o h_s.o h_b.o bath_terms.o \
	initial_state.o 

objects = main.o params.o

spin_boson : $(objects)
	g++ -o spin_boson $(objects)


.PHONY : clean

clean :
	-rm spin_boson $(objects)
