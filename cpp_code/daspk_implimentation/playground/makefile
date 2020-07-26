
objects = testCmain.o testF.o

testC : $(objects)
	g++ -o testC $(objects) -lgfortran 

testCmain.o : testCmain.cpp
	g++ -c -o testCmain.o testCmain.cpp

testF.o : testF.for
	gfortran -c -o testF.o testF.for


