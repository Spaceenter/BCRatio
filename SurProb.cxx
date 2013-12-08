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
float Cutoff;
int TOFNHit;
float Chi2T;
float Chi2C;
float TOFCl;
float Beta;
int ZTOF;
float ZProb;
float TOFAsyUL;
float TOFAsyMM;
float qtofl1;
float qtofl2;
float qtofl3;
float qtofl4;
float QTrInner;
float QTr1[6];
float dsl1[6];
float qtrd;
float qtrdu;
float qtrdl;
float AsyE;

//Read data
TChain *chain;
void ReadData()
{
  chain = new TChain("sptree");
  chain->Add("/afs/cern.ch/work/w/weisun/public/data/surp/surp.root");

  chain->SetBranchAddress("RunN",&RunN);
  chain->SetBranchAddress("EventN",&EventN);
  chain->SetBranchAddress("Cutoff",&Cutoff);
  chain->SetBranchAddress("TOFNHit",&TOFNHit);
  chain->SetBranchAddress("Chi2T",&Chi2T);
  chain->SetBranchAddress("Chi2C",&Chi2C);
  chain->SetBranchAddress("TOFCl",&TOFCl);
  chain->SetBranchAddress("Beta",&Beta);
  chain->SetBranchAddress("ZTOF",&ZTOF);
  chain->SetBranchAddress("ZProb",&ZProb);
  chain->SetBranchAddress("TOFAsyUL",&TOFAsyUL);
  chain->SetBranchAddress("TOFAsyMM",&TOFAsyMM);
  chain->SetBranchAddress("qtofl1",&qtofl1);
  chain->SetBranchAddress("qtofl2",&qtofl2);
  chain->SetBranchAddress("qtofl3",&qtofl3);
  chain->SetBranchAddress("qtofl4",&qtofl4);
  chain->SetBranchAddress("QTrInner",&QTrInner);
  chain->SetBranchAddress("QTr1",&QTr1);
  chain->SetBranchAddress("dsl1",&dsl1);
  chain->SetBranchAddress("qtrd",&qtrd);
  chain->SetBranchAddress("qtrdu",&qtrdu);
  chain->SetBranchAddress("qtrdl",&qtrdl);
  chain->SetBranchAddress("AsyE",&AsyE);
}

//**********************************
//Define variables 
//**********************************

//Binning
float rbin[3] = {10,0.4,1.4};
float ubin[3] = {90,0,9};
float lbin[3] = {110,-2,9};
char name[100];

//Defined variables
float LogCutoff;
float value;
float weight;

//Cuts
bool cut12, cut12_h, cut12_hh, tofcut;
bool sel5[6], sel6[6];

//Define variables
void DefineVariables()
{
  //Defined variables
  LogCutoff = TMath::Log10(Cutoff);

  //Cuts
  cut12 = QTrInner<2.6;
  cut12 &= qtofl1+qtofl2<7;
  cut12 &= qtrd<3.5;
  cut12 &= qtrdu+qtrdl<7;
  if(QTrInner<0) cut12 &= qtofl3+qtofl4<7;

  cut12_h = QTrInner<2.6;
  cut12_h &= qtofl1+qtofl2<7;
  cut12_h &= qtrd<4;
  cut12_h &= qtrdu+qtrdl<8;
  if(QTrInner<0) cut12_h &= qtofl3+qtofl4<7;

  cut12_hh = QTrInner<2.6;
  cut12_hh &= qtofl1+qtofl2<8;
  cut12_hh &= qtrd<5;
  cut12_hh &= qtrdu+qtrdl<10;
  if(QTrInner<0) cut12_hh &= qtofl3+qtofl4<8;

  for(int i=0;i<6;i++)
  {
    sel5[i] = QTr1[i]>4.9 && QTr1[i]<5.1;
    sel6[i] = QTr1[i]>5.9 && QTr1[i]<6.1;
  }

  tofcut = TOFNHit==4 && ZProb>0.95 && Chi2C<10 && Chi2T<10 && TOFCl>0.7;
}

//**********************************
//Define histograms and save results
//**********************************

//L1 templates (Z=1,2) and spectrum
TH1D *ht[6][10], *hs0[6][10], *hmt[6][10], *hms[6][10];

//Detector - TrInner
TH1D *hqin5[6][10], *hqin6[6][10];

//Survival probability and ratio
TH1D *hmsp5[6], *hmsp6[6], *hmsp8[6], *hmspr56[6];

//Survival probability and ratio for different window size at LogCutoff>0.4
TH1D *hwspr56, *hmwspr56;

//Charge distributions for cuts study
TH1D *hctrin, *hctrd;
TH2D *hctofu, *hctofl, *hctofl2, *hctrdh;

//Track reconstruction efficiency
TH1D *hetr[4], *hetrr;

//Define histograms
void DefHist()
{
  for(int ws=0;ws<6;ws++)
  {
    for(int i=0;i<rbin[0];i++)
    {
      sprintf(name,"ht_%d_%d",ws,i+4);
      ht[ws][i] = new TH1D(name,name,ubin[0],ubin[1],ubin[2]);
      sprintf(name,"hs0_%d_%d",ws,i+4);
      hs0[ws][i] = new TH1D(name,name,ubin[0],ubin[1],ubin[2]);
      sprintf(name,"hmt_%d_%d",ws,i+4);
      hmt[ws][i] = new TH1D(name,name,ubin[0],ubin[1],ubin[2]);
      sprintf(name,"hms_%d_%d",ws,i+4);
      hms[ws][i] = new TH1D(name,name,ubin[0],ubin[1],ubin[2]);
      sprintf(name,"hqin5_%d_%d",ws,i+4);
      hqin5[ws][i] = new TH1D(name,name,lbin[0],lbin[1],lbin[2]);
      sprintf(name,"hqin6_%d_%d",ws,i+4);
      hqin6[ws][i] = new TH1D(name,name,lbin[0],lbin[1],lbin[2]);
    }
    sprintf(name,"hmsp5_%d",ws);
    hmsp5[ws] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
    sprintf(name,"hmsp6_%d",ws);
    hmsp6[ws] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
    sprintf(name,"hmspr56_%d",ws);
    hmspr56[ws] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
  }
  sprintf(name,"hctrin");
  hctrin = new TH1D(name,name,lbin[0],lbin[1],lbin[2]);
  sprintf(name,"hctrd");
  hctrd = new TH1D(name,name,ubin[0],ubin[1],ubin[2]);
  sprintf(name,"hctofu");
  hctofu = new TH2D(name,name,ubin[0],ubin[1],ubin[2],ubin[0],ubin[1],ubin[2]);
  sprintf(name,"hctofl");
  hctofl = new TH2D(name,name,ubin[0],ubin[1],ubin[2],ubin[0],ubin[1],ubin[2]);
  sprintf(name,"hctofl2");
  hctofl2 = new TH2D(name,name,ubin[0],ubin[1],ubin[2],ubin[0],ubin[1],ubin[2]);
  sprintf(name,"hctrdh");
  hctrdh = new TH2D(name,name,ubin[0],ubin[1],ubin[2],ubin[0],ubin[1],ubin[2]);
  sprintf(name,"hwspr56");
  hwspr56 = new TH1D(name,name,6,9.5,15.5);
  sprintf(name,"hmwspr56");
  hmwspr56 = new TH1D(name,name,6,9.5,15.5);
  for(int i=0; i<4; i++)
  {
    sprintf(name,"hetr_%d",i);
    hetr[i] = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
  }
  sprintf(name,"hetrr");
  hetrr = new TH1D(name,name,rbin[0],rbin[1],rbin[2]);
}

//Save results
void SaveResults()
{
  //Open file
  TFile *ff=new TFile("sur_prob.root","RECREATE");
  ff->cd();

  //Write hisograms
  for(int ws=0;ws<6;ws++)
  {
    for(int i=0;i<rbin[0];i++)
    {
      ht[ws][i]->Write();
      hs0[ws][i]->Write();
      hmt[ws][i]->Write();
      hms[ws][i]->Write();
      hqin5[ws][i]->Write();
      hqin6[ws][i]->Write();
    }
    hmsp5[ws]->Write();
    hmsp6[ws]->Write();
    hmspr56[ws]->Write();
  }
  hctrin->Write();
  hctrd->Write();
  hctofu->Write();
  hctofl->Write();
  hctofl2->Write();
  hctrdh->Write();
  hmwspr56->Write();
  for(int i=0; i<4; i++) hetr[i]->Write();
  hetrr->Write();

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

    //Remove very low energy stuff
    if(LogCutoff<0.4) continue;

    //Determine weight
    if(QTr1[5]<1.6) weight = 500;
    else if(QTr1[5]<3.5) weight = 50;
    else weight = 1;

    //For defining cuts
    hctrin->Fill(QTrInner,weight);
    hctofu->Fill(qtofl1,qtofl2,weight);
    hctofl->Fill(qtofl3,qtofl4,weight);
    if(QTrInner<0) hctofl2->Fill(qtofl3,qtofl4,weight);
    hctrd->Fill(qtrd,weight);
    hctrdh->Fill(qtrdu,qtrdl,weight);

    //For estimating survival probability
    for(int ws=0;ws<6;ws++)
    {
      for(int i=0;i<rbin[0];i++)
      {
	double LogCutoffCut = 0.4+i*0.1;
	if(LogCutoff<LogCutoffCut) continue;

	//Charge 1 2 template on L1
	if(cut12 && LogCutoffCut<1) ht[ws][i]->Fill(QTr1[ws],weight);
	if(cut12_h && LogCutoffCut>=1 && LogCutoffCut<1.2) ht[ws][i]->Fill(QTr1[ws],weight);
	if(cut12_hh && LogCutoffCut>=1.2) ht[ws][i]->Fill(QTr1[ws],weight);

	//Spectrum on L1
	hs0[ws][i]->Fill(QTr1[ws],weight);

	//QTrInner
	if(sel5[ws]) hqin5[ws][i]->Fill(QTrInner);
	if(sel6[ws]) hqin6[ws][i]->Fill(QTrInner);
      }
    }

    //Track reconstruction efficiency
    if(tofcut)
    {
      if(ZTOF==5 && QTr1[5]>4.5 && QTr1[5]<5.2)
      {
	hetr[0]->Fill(LogCutoff);
	if(QTrInner>4 && QTrInner<6) hetr[1]->Fill(LogCutoff);
      }
      if(ZTOF==6 && QTr1[5]>5.5 && QTr1[5]<6.5)
      {
	hetr[2]->Fill(LogCutoff);
	if(QTrInner>5 && QTrInner<7) hetr[3]->Fill(LogCutoff);
      }
    }
  }
}

//**********************************
//Retouch
//**********************************

void Retouch()
{
  //Spectrum on L1 subtracting charge 1 2
  for(int ws=0;ws<6;ws++)
  {
    for(int i=0;i<rbin[0];i++)
    {
      double scale = hs0[ws][i]->Integral(9,12)/ht[ws][i]->Integral(9,12);
      ht[ws][i]->Scale(scale);

      //Modified version with cubic spline interpolation
      float xx[6] = {5.1, 5.3, 5.5, 6.5, 6.7, 6.9};
      float yy[6];
      yy[0] = 0.5*ht[ws][i]->Integral(51,52);
      yy[1] = 0.5*ht[ws][i]->Integral(53,54);
      yy[2] = 0.5*ht[ws][i]->Integral(55,56);
      yy[3] = 0.5*ht[ws][i]->Integral(65,66);
      yy[4] = 0.5*ht[ws][i]->Integral(66,67);
      yy[5] = 0.5*ht[ws][i]->Integral(67,68);
      TGraph *grc = new TGraph(6,xx,yy);
      for(int j=1;j<=ubin[0];j++)
      {
	if(j<=56 || j>=65) hmt[ws][i]->SetBinContent(j,ht[ws][i]->GetBinContent(j));
	else hmt[ws][i]->SetBinContent(j,grc->Eval(j/10.0-0.05,0,"S"));
	value = hs0[ws][i]->GetBinContent(j)-hmt[ws][i]->GetBinContent(j);
	hms[ws][i]->SetBinContent(j,value);
      }
    }
  }

  //Survival probability without considering contamination of Z>=3
  float a, b, ea, eb;
  for(int ws=0;ws<6;ws++)
  {
    for(int i=0;i<rbin[0];i++)
    {
      //B - Modified
      a=hqin5[ws][i]->Integral(66,75);
      b=hms[ws][i]->Integral(50,51)-hqin5[ws][i]->Integral(76,85);
      ea=TMath::Sqrt(a);
      eb=TMath::Sqrt(b);
      hmsp5[ws]->SetBinContent(i+1,a/b);
      hmsp5[ws]->SetBinError(i+1,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));

      //C - Modified
      a=hqin6[ws][i]->Integral(76,85);
      b=hms[ws][i]->Integral(60,61);
      ea=TMath::Sqrt(a);
      eb=TMath::Sqrt(b);
      hmsp6[ws]->SetBinContent(i+1,a/b);
      hmsp6[ws]->SetBinError(i+1,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));

      //Ratio of B and C - Modified
      a=hmsp5[ws]->GetBinContent(i+1);
      b=hmsp6[ws]->GetBinContent(i+1);
      ea=hmsp5[ws]->GetBinError(i+1);
      eb=hmsp6[ws]->GetBinError(i+1);
      hmspr56[ws]->SetBinContent(i+1,a/b);
      hmspr56[ws]->SetBinError(i+1,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));
    }

    //Different window size at LogCutoff>0.4
    hmwspr56->SetBinContent(ws+1,hmspr56[ws]->GetBinContent(1));
    hmwspr56->SetBinError(ws+1,hmspr56[ws]->GetBinError(1));
  }

  //Track reconstruction efficiency 
  for(int i=1; i<=rbin[0]; i++)
  {
    //B
    a=hetr[1]->GetBinContent(i);
    b=hetr[0]->GetBinContent(i);
    ea=TMath::Sqrt(a);
    eb=TMath::Sqrt(b);
    hetr[1]->SetBinContent(i,a/b);
    hetr[1]->SetBinError(i,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));

    //C
    a=hetr[3]->GetBinContent(i);
    b=hetr[2]->GetBinContent(i);
    ea=TMath::Sqrt(a);
    eb=TMath::Sqrt(b);
    hetr[3]->SetBinContent(i,a/b);
    hetr[3]->SetBinError(i,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));

    //Ratio of B and C
    a=hetr[1]->GetBinContent(i);
    b=hetr[3]->GetBinContent(i);
    ea=hetr[1]->GetBinError(i);
    eb=hetr[3]->GetBinError(i);
    hetrr->SetBinContent(i,a/b);
    hetrr->SetBinError(i,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));
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
