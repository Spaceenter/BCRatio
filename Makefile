CXX = g++
INCLUDE = `root-config --cflags` -I/afs/cern.ch/user/w/weisun/Software/cfitsio/include/
LIBRARY = `root-config --glibs` -lRooFit -lRooFitCore -lFoam -lMinuit -L/afs/cern.ch/user/w/weisun/Software/cfitsio/lib/ -lcfitsio

default: BCRatio

BCRatio: BCRatio.o AnalysisM.o
	$(CXX) BCRatio.o AnalysisM.o $(LIBRARY) -o BCRatio.exe

BCRatio.o: BCRatio.cxx
	$(CXX) $(INCLUDE) -c BCRatio.cxx

AnalysisM.o: AnalysisM.cxx
	$(CXX) $(INCLUDE) -c AnalysisM.cxx

ResMatrix: ResMatrix.o
	$(CXX) ResMatrix.o $(LIBRARY) -o ResMatrix.exe

ResMatrix.o: ResMatrix.cxx
	$(CXX) $(INCLUDE) -c ResMatrix.cxx

SurProb: SurProb.o
	$(CXX) SurProb.o $(LIBRARY) -o SurProb.exe

SurProb.o: SurProb.cxx
	$(CXX) $(INCLUDE) -c SurProb.cxx

IntMC: IntMC.o
	$(CXX) IntMC.o $(LIBRARY) -o IntMC.exe

IntMC.o: IntMC.cxx
	$(CXX) $(INCLUDE) -c IntMC.cxx

clean:
	rm -rf *.o
	rm -rf *.exe

allclean:
	rm -rf *.o
	rm -rf *.exe
	rm -rf *.root
	rm -rf *.txt
	rm -rf *.job
	rm -rf core.*
