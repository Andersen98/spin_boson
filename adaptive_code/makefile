
BOOST_INCLUDE=/home/ethan/Documents/boost_1_74_0

BOOST_LIBS=/home/ethan/Documents/boost_1_74_0/stage/lib

objects=main.o

spin_boson : $(objects)
	g++ -o spin_boson $(objects) 

main.o: main.cpp
	g++ -c main.cpp -L $(BOOST_LIBS) -I $(BOOST_INCLUDE)


.PHONY : clean
clean :
	rm spin_boson $(objects)
