
objects = main.o func.o

main : $(objects)
	g++ -o main $(objects) -lgfortran

main.o : main.cpp
	g++ -c -o main.o main.cpp

func.o : func.for
	gfortran -c -o func.o func.for

.PHONY : clean

clean :
	-rm main $(objects)
