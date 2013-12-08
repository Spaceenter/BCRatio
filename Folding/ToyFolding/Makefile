CXX = g++
INCLUDE = `root-config --cflags` 
LIBRARY = `root-config --glibs` -lRooFit -lRooFitCore -lFoam -lMinuit 

default: toy

toy: toy.o
	$(CXX) toy.o $(LIBRARY) -o toy.exe

toy.o: toy.cxx
	$(CXX) $(INCLUDE) -c toy.cxx

clean:
	rm -rf *.o
	rm -rf *.exe
