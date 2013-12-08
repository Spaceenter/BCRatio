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
int RunN;
int EventN;
int NTOFHit;
double Chi2T;
double Chi2C;
double TOFCl;
double Beta;
int NCounter[4];
double TOFEdepL[4];
int ZTOF;
double ZProb;
double TOFAsyUL;
double TOFAsyMM;
int IsMatch;
double Rig[3];
double Chisq[3];
double QTr1;
double QTr9;
double QTrInner;
int TrYPat;
double TrCl[2];
double TrInAsyE;
int NACC;
int NMCEv;
double MomMC;
double QMC;
double MassMC;
double DirMC[3];
int SelInfo;
int PhyT;
int TrGeoPat;
double IntCoo[3];
double IntDir[3];
int IntPID;
double IntMom;

//Read data
TChain *chain;
void ReadData()
{
  chain = new TChain("QTree");
  chain->Add("/afs/cern.ch/work/w/weisun/public/data/mc/dpmjet.root");

  chain->SetBranchAddress("RunN",&RunN);
  chain->SetBranchAddress("EventN",&EventN);
  chain->SetBranchAddress("NTOFHit",&NTOFHit);
  chain->SetBranchAddress("Chi2T",&Chi2T);
  chain->SetBranchAddress("Chi2C",&Chi2C);
  chain->SetBranchAddress("TOFCl",&TOFCl);
  chain->SetBranchAddress("Beta",&Beta);
  chain->SetBranchAddress("NCounter",&NCounter);
  chain->SetBranchAddress("TOFEdepL",&TOFEdepL);
  chain->SetBranchAddress("ZTOF",&ZTOF);
  chain->SetBranchAddress("ZProb",&ZProb);
  chain->SetBranchAddress("TOFAsyUL",&TOFAsyUL);
  chain->SetBranchAddress("TOFAsyMM",&TOFAsyMM);
  chain->SetBranchAddress("IsMatch",&IsMatch);
  chain->SetBranchAddress("Rig",&Rig);
  chain->SetBranchAddress("Chisq",&Chisq);
  chain->SetBranchAddress("QTr1",&QTr1);
  chain->SetBranchAddress("QTr9",&QTr9);
  chain->SetBranchAddress("QTrInner",&QTrInner);
  chain->SetBranchAddress("TrYPat",&TrYPat);
  chain->SetBranchAddress("TrCl",&TrCl);
  chain->SetBranchAddress("TrInAsyE",&TrInAsyE);
  chain->SetBranchAddress("NACC",&NACC);
  chain->SetBranchAddress("NMCEv",&NMCEv);
  chain->SetBranchAddress("MomMC",&MomMC);
  chain->SetBranchAddress("QMC",&QMC);
  chain->SetBranchAddress("MassMC",&MassMC);
  chain->SetBranchAddress("DirMC",&DirMC);
  chain->SetBranchAddress("SelInfo",&SelInfo);
  chain->SetBranchAddress("PhyT",&PhyT);
  chain->SetBranchAddress("TrGeoPat",&TrGeoPat);
  chain->SetBranchAddress("IntCoo",&IntCoo);
  chain->SetBranchAddress("IntDir",&IntDir);
  chain->SetBranchAddress("IntPID",&IntPID);
  chain->SetBranchAddress("IntMom",&IntMom);
}

//**********************************
//Define variables 
//**********************************

//Binning
double qbin[3] = {200,0,10};
double xbin[3] = {100,0,10};
double rbin[3] = {18,0.2,2};
char name[100];

//Defined variables
double LogRig;
int iBinRig;
bool pb10, pb11, pc12, po16;
bool TrigCut, RecCut, PreCut, TempCut, TOFCut;

//Define variables
void DefineVariables()
{
  //Binning
  LogRig = TMath::Log10(MomMC/QMC);
  iBinRig = (int)(10*LogRig)-2;

  //Particle
  pb10 = QMC==5 && MassMC<10;
  pb11 = QMC==5 && MassMC>10;
  pc12 = QMC==6;
  po16 = QMC==8;

  //Cut
  //TrigCut = (PhyT&0x2) || (PhyT&0x4) || (PhyT&0x10);
  TrigCut = (PhyT&0x4) || (PhyT&0x8);
  RecCut = SelInfo==0;
  PreCut = NTOFHit==4 && IsMatch;
  TempCut = IntCoo[2]<-50 && QTrInner>QMC*0.75;
  TOFCut = Chi2T<15 && Chi2C<15 && Beta>0.4 && Beta<1.2 && ZProb>0.9;
  TOFCut &= TOFCl>0.7 && TMath::Abs(TOFAsyUL)<0.3 && TOFAsyMM<0.4;
}


//**********************************
//Define histograms and save results
//**********************************

//Survival probability
TH1D *htot[4], *hrec[4][18], *htmp[4][18];

//BetaH reconstruction efficiency
TH1D *hetof[4];

//Define histograms
void DefHist()
{
  for(int i=0;i<4;i++)
  {
    sprintf(name,"htot_%d",i);
    htot[i] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
    for(int j=0;j<18;j++)
    {
      sprintf(name,"hrec_%d_%d",i,j+2);
      hrec[i][j] = new TH1D(name,name,qbin[0],qbin[1],qbin[2]);
      sprintf(name,"htmp_%d_%d",i,j+2);
      htmp[i][j] = new TH1D(name,name,qbin[0],qbin[1],qbin[2]);
    }
    sprintf(name,"hetof_%d",i);
    hetof[i] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
  }
}

//Save results
void SaveResults()
{
  //Open file
  TFile *ff=new TFile("int_mc_dpmjet.root","RECREATE");
  ff->cd();

  //Write Hisogram
  for(int i=0;i<4;i++)
  {
    htot[i]->Write();
    for(int j=0;j<18;j++)
    {
      hrec[i][j]->Write();
      htmp[i][j]->Write();
    }
    hetof[i]->Write();
  }

  //Close file
  ff->Write();
  ff->Close();
}

//**********************************
//Event loop
//**********************************

void DoEventLoop()
{
  for(int index=0;index<chain->GetEntries();index++)
  {
    chain->GetEntry(index);

    //Define variables
    DefineVariables();

    //Prelimilary cuts - within acceptance of L1 to L8
    if(LogRig<0.2 || LogRig>=2) continue; //Rigidity range
    if(TrigCut==false) continue;

    //Total
    if(pb10) htot[0]->Fill(LogRig);
    if(pb11) htot[1]->Fill(LogRig);
    if(pc12) htot[2]->Fill(LogRig);
    if(po16) htot[3]->Fill(LogRig);

    //QTrInner spectrum and template
    if(RecCut && PreCut)
    {
      //With TrTrack reconstructed
      if(pb10) hrec[0][iBinRig]->Fill(QTrInner);
      if(pb11) hrec[1][iBinRig]->Fill(QTrInner);
      if(pc12) hrec[2][iBinRig]->Fill(QTrInner);
      if(po16) hrec[3][iBinRig]->Fill(QTrInner);

      //Templates
      if(TempCut)
      {
	if(pb10) htmp[0][iBinRig]->Fill(QTrInner);
	if(pb11) htmp[1][iBinRig]->Fill(QTrInner);
	if(pc12) htmp[2][iBinRig]->Fill(QTrInner);
	if(po16) htmp[3][iBinRig]->Fill(QTrInner);
      }
    }

    //BetaH reconstruction efficiency
    if(SelInfo!=2)
    {
      if(pb10) hetof[0]->Fill(LogRig);
      if(pb11) hetof[1]->Fill(LogRig);
      if(pc12) hetof[2]->Fill(LogRig);
      if(po16) hetof[3]->Fill(LogRig);
    }
  }
}


//**********************************
//Retouch
//**********************************

void Retouch()
{
  float a, b, ea, eb;
  for(int i=0;i<4;i++)
  {
    for(int r=1;r<=rbin[0];r++)
    {
      a=hetof[i]->GetBinContent(r);
      b=htot[i]->GetBinContent(r);
      ea=TMath::Sqrt(a);
      eb=TMath::Sqrt(b);
      hetof[i]->SetBinContent(r,a/b);
      hetof[i]->SetBinError(r,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));
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
  DoEventLoop();
  Retouch();
  SaveResults();

  return 0;
}
