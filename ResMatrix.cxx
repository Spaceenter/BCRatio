#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

//**********************************
//Read data
//**********************************

//TTree variables
Int_t RunN;
Int_t EventN;
Float_t Chi2T;
Float_t Chi2C;
Float_t TOFCl;
Float_t Beta;
Int_t ZTOF;
Float_t ZProb;
Float_t Rig[2];
Float_t Chisq[2];
Float_t QTr1;
Float_t QTr9;
Float_t QTrInner;
Int_t TrYPat;
Float_t TrCl;
Float_t TrInAsyE;
Int_t NACC;
Float_t MomMC;
Int_t TrGeoPat;

//Read data
TChain *chain;
void ReadData()
{
  chain = new TChain("QTree");
  chain->Add("/afs/cern.ch/work/w/weisun/public/data/mc/unmche.root");

  chain->SetBranchAddress("RunN",&RunN);
  chain->SetBranchAddress("EventN",&EventN);
  chain->SetBranchAddress("Chi2T",&Chi2T);
  chain->SetBranchAddress("Chi2C",&Chi2C);
  chain->SetBranchAddress("TOFCl",&TOFCl);
  chain->SetBranchAddress("Beta",&Beta);
  chain->SetBranchAddress("ZTOF",&ZTOF);
  chain->SetBranchAddress("ZProb",&ZProb);
  chain->SetBranchAddress("Rig",&Rig);
  chain->SetBranchAddress("Chisq",&Chisq);
  chain->SetBranchAddress("QTr1",&QTr1);
  chain->SetBranchAddress("QTr9",&QTr9);
  chain->SetBranchAddress("QTrInner",&QTrInner);
  chain->SetBranchAddress("TrYPat",&TrYPat);
  chain->SetBranchAddress("TrCl",&TrCl);
  chain->SetBranchAddress("TrInAsyE",&TrInAsyE);
  chain->SetBranchAddress("NACC",&NACC);
  chain->SetBranchAddress("MomMC",&MomMC);
  chain->SetBranchAddress("TrGeoPat",&TrGeoPat);
}

//**********************************
//Define variables 
//**********************************

int YL[9], InNHitsY;
bool YGood, TOFCut, TrCut, SelCut;
double Rg, Rm, Rm9;

void DefineVariables()
{
  //TrTrack Y pattern
  InNHitsY=0;
  for(int i=0;i<9;i++)
  {
    if((TrYPat&(1<<i))==(1<<i))
    {
      YL[i]=1;
      if(i!=0 && i!=8) InNHitsY++;
    }
    else YL[i]=0;
  }
  if(YL[0] && YL[1] && (YL[2]||YL[3]) && (YL[4]||YL[5]) && (YL[6]||YL[7])) YGood = true; 
  else YGood = false;

  //Selection cuts
  TOFCut = Chi2T<30 && Chi2C<30 && ZProb>0.7;
  TrCut = Chisq[1]<20 && YGood && InNHitsY>=5 && TrCl>0.4 && TrInAsyE<0.8;
  SelCut = TOFCut && TrCut;

  //Rigidity
  Rg = MomMC/2.0;
  Rm = Rig[1];
  Rm9 = Rig[0];
}

//**********************************
//Define histograms and save results
//**********************************

//Binning of Rigidity - for unfolding
const int NBinRigUn = 43;
double XBinRigUn[NBinRigUn+1] = {
  1.091, 1.259, 1.466, 1.718, 2.013, 2.339, 2.688, 3.098, 3.561, 4.066, 4.657, 5.326, 6.048, 6.885, 7.847, 8.915, 10.123, 11.535, 13.148, 14.955, 17.105, 19.614, 22.448, 25.872, 29.947, 34.641, 40.417, 47.399, 55.747, 66.231, 79.050, 95.110, 114.847, 139.226, 170.497, 210.249, 263.172, 333.828, 431.453, 577.082, 819.239, 1341.318, 3271.122, 10000
};
double XBinRigUnInv[2*NBinRigUn+1];

TH1D *hres[10], *hres9[10];
TH2D *hmatrix, *hmatrix9;
TGraph *gr, *gr9;

//Define histograms
void DefHist()
{
  char name[100];
  for(int i=0; i<10; i++)
  {
    sprintf(name,"hres_%d",(int)pow(2.0,i+1));
    hres[i] = new TH1D(name,name,200,-5,5);
    sprintf(name,"hres9_%d",(int)pow(2.0,i+1));
    hres9[i] = new TH1D(name,name,200,-5,5);
  }

  for(int i=0; i<NBinRigUn; i++) XBinRigUnInv[i] = -1.0/XBinRigUn[i];
  XBinRigUnInv[NBinRigUn] = 0;
  for(int i=0; i<NBinRigUn; i++) XBinRigUnInv[NBinRigUn+1+i] = -XBinRigUnInv[NBinRigUn-1-i];
  hmatrix = new TH2D("hmatrix","hmatrix",2*NBinRigUn,XBinRigUnInv,2*NBinRigUn,XBinRigUnInv);
  hmatrix9 = new TH2D("hmatrix9","hmatrix9",2*NBinRigUn,XBinRigUnInv,2*NBinRigUn,XBinRigUnInv);
}

//Save results
void SaveResults()
{
  //Open file
  TFile *ff=new TFile("matrix.root","RECREATE");
  ff->cd();

  //Write hisograms
  for(int i=0; i<10; i++) 
  {
    hres[i]->Write();
    hres9[i]->Write();
  }
  hmatrix->Write();
  hmatrix9->Write();
  gr->Write();
  gr9->Write();

  //Close file
  ff->Write();
  ff->Close();
}

//**********************************
//Event loop
//**********************************

void LoopResHist()
{
  for(int index=0;index<chain->GetEntries();index++)
  {
    chain->GetEntry(index);

    //Define variables
    DefineVariables();

    //Selection cuts
    if(SelCut==false) continue;

    for(int i=0; i<10; i++)
    {
      double RPoint = pow(2.0,i+1);
      if(Rg>RPoint*0.8 && Rg<RPoint*1.2)
      {
	hres[i]->Fill((1/Rm-1/Rg)*Rg);
	break;
      }
    }

    if(YL[8]==0) continue;

    for(int i=0; i<10; i++)
    {
      double RPoint = pow(2.0,i+1);
      if(Rg>RPoint*0.8 && Rg<RPoint*1.2)
      {
	hres9[i]->Fill((1/Rm9-1/Rg)*Rg);
	break;
      }
    }
  }

  double scale;
  for(int i=0; i<10; i++)
  {
    scale = hres[i]->Integral();
    hres[i]->Scale(1.0/scale);
    scale = hres9[i]->Integral();
    hres9[i]->Scale(1.0/scale);
  }
}

//**********************************
//Generate resolution matrix
//**********************************

void GenerateResMatrix(double Expansion)
{
  double rigP[10] = {2,4,8,16,32,64,128,256,512,1024};
  double sigma1P[10] = {0.1099,0.0937,0.0897,0.0923,0.1104,0.1561,0.2118,0.2965,0.4478,0.7631};
  double sigma1_256 = 0.2965; 
  double sigma2_256 = 0.7242;
  double weight21_256 = 0.308;

  gr = new TGraph(10,rigP,sigma1P);
  gr->SetName("gr");

  TF1 *ff = new TF1("ff","(TMath::Exp(-x*x/2/[0]/[0])/[0]+[2]/[1]*TMath::Exp(-x*x/2/[1]/[1]))",-10,10);

  for(int i=0; i<NBinRigUn; i++) //Positive side of the matrix
  {
    int iBin = NBinRigUn+i+1;
    double AveRigInv = 0.5*(XBinRigUnInv[iBin-1]+XBinRigUnInv[iBin]);

    double sigma1 = gr->Eval(1.0/AveRigInv,0,"S");
    double sigma2 = sigma1/sigma1_256*sigma2_256;
    double weight21 = weight21_256;
    ff->SetParameter(0,sigma1*Expansion);
    ff->SetParameter(1,sigma2*Expansion);
    ff->SetParameter(2,weight21);

    double LeftInt, RightInt;
    for(int j=0; j<2*NBinRigUn; j++)
    {
      LeftInt = XBinRigUnInv[j]/AveRigInv-1;
      RightInt = XBinRigUnInv[j+1]/AveRigInv-1;
      hmatrix->SetBinContent(iBin,j+1,ff->Integral(LeftInt,RightInt));
    }
  }

  for(int i=0; i<NBinRigUn; i++) //Negative side of the matrix
  {
    for(int j=0; j<2*NBinRigUn; j++) 
    {
      hmatrix->SetBinContent(i+1,j+1,hmatrix->GetBinContent(2*NBinRigUn-i,2*NBinRigUn-j));
    }
  }

  //Normalization
  for(int i=0; i<2*NBinRigUn; i++) 
  {
    double scale = hmatrix->Integral(i+1,i+1,1,2*NBinRigUn);
    for(int j=0; j<2*NBinRigUn; j++) 
    {
      hmatrix->SetBinContent(i+1,j+1,hmatrix->GetBinContent(i+1,j+1)/scale);
    }
  }
}

void GenerateResMatrix9(double Expansion)
{
  double rigP[10] = {2,4,8,16,32,64,128,256,512,1024};
  double sigma1P[10] = {0.1286,0.1009,0.0924,0.0939,0.1029,0.1204,0.1401,0.1667,0.2057,0.3029};
  double sigma1_512 = 0.1862; 
  double sigma2_512 = 0.5913;
  double weight21_512 = 0.1395;

  gr9 = new TGraph(10,rigP,sigma1P);
  gr9->SetName("gr9");

  TF1 *ff = new TF1("ff","(TMath::Exp(-x*x/2/[0]/[0])/[0]+[2]/[1]*TMath::Exp(-x*x/2/[1]/[1]))",-10,10);

  for(int i=0; i<NBinRigUn; i++) //Positive side of the matrix
  {
    int iBin = NBinRigUn+i+1;
    double AveRigInv = 0.5*(XBinRigUnInv[iBin-1]+XBinRigUnInv[iBin]);

    double sigma1 = gr9->Eval(1.0/AveRigInv,0,"S")/0.2057*0.1862;
    double sigma2 = sigma1/sigma1_512*sigma2_512;
    double weight21 = weight21_512;
    ff->SetParameter(0,sigma1*Expansion);
    ff->SetParameter(1,sigma2*Expansion);
    ff->SetParameter(2,weight21);

    double LeftInt, RightInt;
    for(int j=0; j<2*NBinRigUn; j++)
    {
      LeftInt = XBinRigUnInv[j]/AveRigInv-1;
      RightInt = XBinRigUnInv[j+1]/AveRigInv-1;
      hmatrix9->SetBinContent(iBin,j+1,ff->Integral(LeftInt,RightInt));
    }
  }

  for(int i=0; i<NBinRigUn; i++) //Negative side of the matrix
  {
    for(int j=0; j<2*NBinRigUn; j++) 
    {
      hmatrix9->SetBinContent(i+1,j+1,hmatrix9->GetBinContent(2*NBinRigUn-i,2*NBinRigUn-j));
    }
  }

  //Normalization
  for(int i=0; i<2*NBinRigUn; i++) 
  {
    double scale = hmatrix9->Integral(i+1,i+1,1,2*NBinRigUn);
    for(int j=0; j<2*NBinRigUn; j++) 
    {
      hmatrix9->SetBinContent(i+1,j+1,hmatrix9->GetBinContent(i+1,j+1)/scale);
    }
  }
}

//**********************************
//Main program
//**********************************

int main(int argc,char **argv)
{
  ReadData();
  DefHist();
  LoopResHist();
  GenerateResMatrix(1.05);
  GenerateResMatrix9(1.05);
  SaveResults();

  return 0;
}
