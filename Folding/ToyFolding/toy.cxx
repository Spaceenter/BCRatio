#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"

#define NSTAT 9999999999

using namespace std;

/////////////////////////////////////
//calculate histogram ratio 
/////////////////////////////////////

void CalcHistRatio(TH1D *ha, TH1D *hb, TH1D *hc, bool HasError)
{
  double a, b, ea, eb;
  for(int i=1;i<=ha->GetNbinsX();i++)
  {
    a = ha->GetBinContent(i);
    b = hb->GetBinContent(i);
    if(HasError)
    {    
      ea = ha->GetBinError(i);
      eb = hb->GetBinError(i);
    }    
    else 
    {    
      ea = sqrt(a);
      eb = sqrt(b);
    }    
    if(b!=0) 
    {    
      hc->SetBinContent(i,a/b);
      hc->SetBinError(i,a/b*sqrt(ea*ea/a/a+eb*eb/b/b));
    }
    else
    {
      hc->SetBinContent(i,1);
      hc->SetBinError(i,1);
    }
  }
}

/////////////////////////////////////
//definition of histograms
/////////////////////////////////////

//binning
int NBin = 30;
double XBin[31];

//histograms
TH1D *hbt, *hct, *hbct; //true spectrum
TH1D *hbf, *hcf, *hbcf; //folded spectrum
TH1D *hbcft, *hbcftm; //between hbcf and hbct
TH2D *hmat; //resolution matrix
TGraph *gr; //graph for sigma1

void DefHist()
{
  //binning
  for(int i=0; i<=NBin; i++) XBin[i] = pow(10,i/10.0);

  //definitions
  hbt = new TH1D("hbt","hbt",NBin,XBin);
  hct = new TH1D("hct","hct",NBin,XBin);
  hbct = new TH1D("hbct","hbct",NBin,XBin);
  hbf = new TH1D("hbf","hbf",NBin,XBin);
  hcf = new TH1D("hcf","hcf",NBin,XBin);
  hbcf = new TH1D("hbcf","hbcf",NBin,XBin);
  hbcft = new TH1D("hbcft","hbcft",NBin,XBin);
  hbcftm = new TH1D("hbcftm","hbcftm",NBin,XBin);
  hmat = new TH2D("hmat","hmat",NBin,XBin,NBin,XBin);

  //spectrum
  TF1 *fb = new TF1("fb","0.4*TMath::Power(x,-3)",1,1000);
  TF1 *fc = new TF1("fc","TMath::Power(x,-2.7)",1,1000);
  for(int i=1; i<=NBin; i++)
  {
    hbt->SetBinContent(i,NSTAT*fb->Integral(XBin[i-1],XBin[i]));
    hct->SetBinContent(i,NSTAT*fc->Integral(XBin[i-1],XBin[i]));
  }
  CalcHistRatio(hbt,hct,hbct,false);

  //resolution matrix
  double rigP[10] = {2,4,8,16,32,64,128,256,512,1024};
  double sigma1P[10] = {0.1099,0.0937,0.0897,0.0923,0.1104,0.1561,0.2118,0.2965,0.4478,0.7631};
  double sigma1_256 = 0.2965; 
  double sigma2_256 = 0.7242;
  double weight21_256 = 0.308;
  gr = new TGraph(10,rigP,sigma1P);
  TF1 *ff = new TF1("ff","(TMath::Exp(-x*x/2/[0]/[0])/[0]+[2]/[1]*TMath::Exp(-x*x/2/[1]/[1]))",-10,10);
  for(int i=1; i<=NBin; i++) //Positive side of the matrix
  {
    double XAve = 0.5*(XBin[i-1]+XBin[i]);

    double sigma1 = gr->Eval(XAve,0,"S");
    double sigma2 = gr->Eval(XAve,0,"S")/sigma1_256*sigma2_256;
    double weight21 = weight21_256;
    ff->SetParameter(0,sigma1);
    ff->SetParameter(1,sigma2);
    ff->SetParameter(2,weight21);
    double norm_ff = ff->Integral(-10,10); //matrix normalization 
    for(int j=1; j<=NBin; j++)
    {   
      double LeftInt = XAve/XBin[j]-1;
      double RightInt = XAve/XBin[j-1]-1;
      hmat->SetBinContent(i,j,ff->Integral(LeftInt,RightInt)/norm_ff);
    }   
  }
}

/////////////////////////////////////
//folding
/////////////////////////////////////

void DoFolding()
{
  for(int i=1; i<=NBin; i++)
  {
    for(int j=1; j<=NBin; j++)
    {
      hbf->AddBinContent(j,hbt->GetBinContent(i)*hmat->GetBinContent(i,j));
      hcf->AddBinContent(j,hct->GetBinContent(i)*hmat->GetBinContent(i,j));
    }
  }
  CalcHistRatio(hbf,hcf,hbcf,false);
  CalcHistRatio(hbcf,hbct,hbcft,true);
  for(int i=1; i<=NBin; i++)
  {
    hbcftm->SetBinContent(i,hbcf->GetBinContent(i)-hbct->GetBinContent(i));
  }
}

/////////////////////////////////////
//save results
/////////////////////////////////////

void SaveResults()
{
  TFile *f = new TFile("result.root","RECREATE");
  f->cd();

  hbt->Write();
  hct->Write();
  hbct->Write();
  hbf->Write();
  hcf->Write();
  hbcf->Write();
  hbcft->Write();
  hbcftm->Write();
  hmat->Write();

  f->Write();
  f->Close();
}

/////////////////////////////////////
//main function
/////////////////////////////////////

int main(int argc, char **argv)
{
  DefHist();
  DoFolding();
  SaveResults();
  return 0;
}
