CXX = g++
INCLUDE = `root-config --cflags` 
LIBRARY = `root-config --glibs` -lRooFit -lRooFitCore -lFoam -lMinuit 

default: fold

fold: fold.o
	$(CXX) fold.o $(LIBRARY) -o fold.exe

fold.o: fold.cxx
	$(CXX) $(INCLUDE) -c fold.cxx

clean:
	rm -rf *.o
	rm -rf *.exe
