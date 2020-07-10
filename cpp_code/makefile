objects = main.o kbd.o command.o display.o \
	h_sb.o h_s.o h_b.o bath_terms.o \
	initial_state.o 

objects = main.o params.o psi_naught.o

spin_boson : $(objects)
	g++ -o spin_boson $(objects)

main.o : psi_naught.hpp params.hpp

psi_naught.o : psi_naught.hpp params.hpp

params.o : params.hpp

.PHONY : clean

clean :
	-rm spin_boson $(objects)
