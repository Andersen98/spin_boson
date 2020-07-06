objects = main.o kbd.o command.o display.o \
	h_sb.o h_s.o h_b.o bath_terms.o \
	initial_state.o 

spin_boson : main.o
	g++ -o spin_boson main.o


.PHONY : clean

clean :
	-rm spin_boson
