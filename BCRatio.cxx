#include "AnalysisM.h"

using namespace std;

int main(int argc, char *argv[])
{
  AnalysisM *pBC = new AnalysisM();

  pBC->ReadCuts(argv[1]);
  pBC->GenerateEffCorr();
  pBC->CalcBCRatio();
  pBC->SaveResults(argv[2]);
  //pBC->StudyCuts();
  //pBC->PlotChargeID();
  //pBC->CalcTrigEffR();

  delete pBC; 

  return 0;
}
