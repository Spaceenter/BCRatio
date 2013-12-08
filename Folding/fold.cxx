#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TMinuit.h"

#define NSTAT 9999999999
#define MB10 9.32699
#define MB11 10.2551 
#define MC12 11.1779
#define FITA 12.5 
#define FITB 1000.0
#define FIXDOF 2 
#define BREAK 258.004

using namespace std;

//calculate histogram ratio 
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

//binning
int NBinRig = 21;
double XBinRig[22] = {
    1.259, 1.718, 2.339, 3.098, 4.066, 5.326, 6.885, 8.915, 11.535, 14.955, 19.614, 25.872, 34.641, 47.399, 66.231, 95.110, 139.226, 210.249, 333.828, 577.082, 1341.318, 3271.122
};
int NBin = 42;
const double XBin[43] = {
  -0.794281, -0.582072, -0.427533, -0.322789, -0.245942, -0.187758, -0.145243, -0.11217, -0.0866927, -0.0668673, -0.050984, -0.0386518, -0.0288675, -0.0210975, -0.0150987, -0.0105141, -0.00718257, -0.00475627, -0.00299555, -0.00173286, -0.000745535, 0, 0.000745535, 0.00173286, 0.00299555, 0.00475627, 0.00718257, 0.0105141, 0.0150987, 0.0210975, 0.0288675, 0.0386518, 0.050984, 0.0668673, 0.0866927, 0.11217, 0.145243, 0.187758, 0.245942, 0.322789, 0.427533, 0.582072, 0.794281
};

//matrix parameters
double rigP[10] = {2,4,8,16,32,64,128,256,512,1024};
double sigma1P[10] = {0.1099,0.0937,0.0897,0.0923,0.1104,0.1561,0.2118,0.2965,0.4478,0.7631};
double sigma1_256 = 0.2965; 
double sigma2_256 = 0.7242;
double weight21_256 = 0.308;

//histograms
TH1D *hbo, *hco, *hbco; //original spectrum
TH1D *hbt, *hct, *hbct; //true spectrum
TH1D *hbf, *hcf, *hbcf; //folded spectrum
TH2D *hmat; //resolution matrix
TGraph *gr; //graph for sigma1
TFile fr("../Draw/result.root");
TFile fr1("../Draw/ExpT.root");
TF1 *fb, *fc; //power law function
TH1D *heffbc; //selection efficiency
TH1D *hrt; //exposure time
TH1D *hsigma; //deviation

//define histograms
void DefHist()
{
  //original spectrum
  hbo = (TH1D*)fr.Get("hunb");
  hco = (TH1D*)fr.Get("hunc");
  hbco = (TH1D*)hbo->Clone("hbco");
  CalcHistRatio(hbo,hco,hbco,false);

  //selection efficiency
  heffbc = (TH1D*)fr.Get("her1");

  //exposure time
  hrt = (TH1D*)fr1.Get("hrt");
  double scale = hrt->GetBinContent(NBinRig);
  for(int i=1; i<=NBinRig; i++) hrt->SetBinContent(i,hrt->GetBinContent(i)/scale);

  //other definitions
  hbt = new TH1D("hbt","hbt",NBin,XBin);
  hct = new TH1D("hct","hct",NBin,XBin);
  hbct = new TH1D("hbct","hbct",NBin,XBin);
  hbf = new TH1D("hbf","hbf",NBin,XBin);
  hcf = new TH1D("hcf","hcf",NBin,XBin);
  hbcf = new TH1D("hbcf","hbcf",NBin,XBin);
  hmat = new TH2D("hmat","hmat",NBin,XBin,NBin,XBin);
  hsigma = new TH1D("hsigma","hsigma",NBin,XBin);
}

//solar modulation
double RigSolarM(double Rig, double Potential, int type)
{
  double Mass, Z;
  if(type==0) 
  {
    Mass = 0.7*MB11 + 0.3*MB10;
    Z = 5;
  }
  else if(type==1) 
  {
    Mass = MC12;
    Z = 6;
  }
  else return 0;
  double E = sqrt(Rig*Rig*Z*Z + Mass*Mass);
  double EM = E + Z*Potential;
  double RigM = sqrt(EM*EM - Mass*Mass)/Z;
  return RigM;
}

//model
void SetTrueSpectra(double *par)
{
  /*
  //B function - extra source
  fb = new TF1("fb","[0]*(TMath::Power(x,[1])+[2]*TMath::Power(x,[3]))",1,5000);
  fb->SetParameter(0,par[0]);
  fb->SetParameter(1,par[1]);
  fb->SetParameter(2,par[4]);
  fb->SetParameter(3,par[2]);
  */

  //B function - break
  fb = new TF1("fb","[0]*(TMath::Power(x,[1])*(1+TMath::Floor(([2]-x)/100000000.0)) + TMath::Power([2],-[3])*TMath::Power(x,[1]+[3])*(1+TMath::Floor((x-[2])/100000000.0)))",1,5000);
  fb->SetParameter(0,par[0]);
  fb->SetParameter(1,par[1]);
  fb->SetParameter(2,BREAK);
  fb->SetParameter(3,par[7]);

  //C function
  fc = new TF1("fc","TMath::Power(x,[0])",1,5000);
  fc->SetParameter(0,par[2]);

  for(int i=NBinRig+1; i<=2*NBinRig; i++)
  {
    //live time correction
    double BinLowEdge = hbt->GetBinLowEdge(i);
    double BinWidth = hbt->GetBinWidth(i);
    double LeftInt = 1.0/(BinLowEdge + BinWidth);
    if(BinLowEdge == 0) BinLowEdge = 1.0/XBinRig[NBinRig];
    double RightInt = 1.0/BinLowEdge;
    double BinCenterRig = (LeftInt + RightInt)/2.0;
    int iBinRig = hrt->FindBin(BinCenterRig);
    double CorrB = hrt->GetBinContent(iBinRig);
    double CorrC = hrt->GetBinContent(iBinRig);

    //solar modulation
    double Potential_0 = par[5]; //B
    double Potential_1 = par[6]; //C
    double LeftIntMB = RigSolarM(LeftInt,Potential_0,0); //B
    double RightIntMB = RigSolarM(RightInt,Potential_0,0); //B
    double BinCenterRigMB = RigSolarM(BinCenterRig,Potential_0,0); //B
    double LeftIntMC = RigSolarM(LeftInt,Potential_1,1); //C
    double RightIntMC = RigSolarM(RightInt,Potential_1,1); //C
    double BinCenterRigMC = RigSolarM(BinCenterRig,Potential_1,1); //C
    CorrB = CorrB*pow(BinCenterRig/BinCenterRigMB,2);
    CorrC = CorrC*pow(BinCenterRig/BinCenterRigMC,2);

    hbt->SetBinContent(i,NSTAT*fb->Integral(LeftIntMB,RightIntMB)*CorrB);
    hct->SetBinContent(i,NSTAT*fc->Integral(LeftIntMC,RightIntMC)*CorrC);
  }
  CalcHistRatio(hbt,hct,hbct,false);
}

//resolution matrix
void SetMatrix(double *par)
{
  gr = new TGraph(10,rigP,sigma1P);
  TF1 *ff = new TF1("ff","(TMath::Exp(-x*x/2/[0]/[0])/[0]+[2]/[1]*TMath::Exp(-x*x/2/[1]/[1]))",-10,10);
  for(int i=0; i<NBinRig; i++) //Positive side of the matrix
  {
    int iBin = NBinRig+i+1;
    double AveRigInv = 0.5*(XBin[iBin-1]+XBin[iBin]);

    double sigma1 = gr->Eval(1.0/AveRigInv,0,"S");
    double sigma2 = gr->Eval(1.0/AveRigInv,0,"S")/sigma1_256*sigma2_256;
    double weight21 = weight21_256; //fixed
    ff->SetParameter(0,sigma1*par[3]);
    ff->SetParameter(1,sigma2*par[3]);
    ff->SetParameter(2,weight21);

    for(int j=0; j<2*NBinRig; j++)
    {
      double LeftInt = XBin[j]/AveRigInv-1;
      double RightInt = XBin[j+1]/AveRigInv-1;
      hmat->SetBinContent(iBin,j+1,ff->Integral(LeftInt,RightInt));
    }
  }
  for(int i=0; i<NBinRig; i++) //Negative side of the matrix
  {
    for(int j=0; j<2*NBinRig; j++)
    {
      hmat->SetBinContent(i+1,j+1,hmat->GetBinContent(2*NBinRig-i,2*NBinRig-j));
    }
  }
  for(int i=0; i<2*NBinRig; i++) //Normalization
  {
    double scale = hmat->Integral(i+1,i+1,1,2*NBinRig);
    for(int j=0; j<2*NBinRig; j++)
    {
      hmat->SetBinContent(i+1,j+1,hmat->GetBinContent(i+1,j+1)/scale);
    }
  }
}

//folding
void DoFolding()
{
  for(int i=1; i<=NBin; i++)
  {
    hbf->SetBinContent(i,0);
    hcf->SetBinContent(i,0);
  }
  for(int i=1; i<=NBin; i++)
  {
    for(int j=1; j<=NBin; j++)
    {
      hbf->AddBinContent(j,hbt->GetBinContent(i)*hmat->GetBinContent(i,j));
      hcf->AddBinContent(j,hct->GetBinContent(i)*hmat->GetBinContent(i,j));
    }
  }
  CalcHistRatio(hbf,hcf,hbcf,false);

  //corrections
  for(int i=NBinRig+1; i<=2*NBinRig; i++)
  {
    double value = hbcf->GetBinContent(i);

    //survival probability and top of instrument corrections
    value = value*1.04+0.005;
  
    //efficiency correction
    double BinLowEdge = hbt->GetBinLowEdge(i);
    double BinWidth = hbt->GetBinWidth(i);
    double LeftInt = 1.0/(BinLowEdge + BinWidth);
    if(BinLowEdge == 0) BinLowEdge = 1.0/XBinRig[NBinRig];
    double RightInt = 1.0/BinLowEdge;
    double BinCenterRig = (LeftInt + RightInt)/2.0;
    int iBinRig = heffbc->FindBin(BinCenterRig);
    value = value*heffbc->GetBinContent(iBinRig);

    hbcf->SetBinContent(i,value);
  }
}

//fcn 
void FCN(Int_t &npar,Double_t *gin,Double_t &f,Double_t *par,Int_t iflag)
{
  SetTrueSpectra(par);
  SetMatrix(par);
  DoFolding();

  double Chi2, Delta, Sigma;
  int StartBin = hbco->FindBin(1/FITB);
  int EndBin = hbco->FindBin(1/FITA);
  Chi2 = 0;
  for(int i=StartBin; i<=EndBin; i++)
  {
    Delta = hbcf->GetBinContent(i) - hbco->GetBinContent(i);
    Sigma = hbco->GetBinError(i);
    Chi2 += Delta*Delta/Sigma/Sigma;
  }
  f = Chi2;

  cout<<par[7]<<": "<<Chi2<<endl;
}

//minimization
void DoMinimization()
{
  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->SetFCN(FCN);
  gMinuit->SetErrorDef(1); //Chi2 

  //variables
  double arglist[10];
  int ierr = 0;
  double start[8] = {0.7375, -3.072, -2.71, 1.07, 0, 1.1, 0.85, 0};
  double step[8] = {0.01, 0.01, 0.01, 0.01, 0.001, 0.01, 0.01, 0.001};
  double low[8] = {0.4, -3.2, -2.8, 0.9, 0, 0.4, 0.4, -2};
  double up[8] = {1, -2.8, -2.6, 1.3, 0.2, 1.4, 1.4, 2};

  //set parameters
  gMinuit->mnparm(0,"ratio",start[0],step[0],low[0],up[0],ierr);
  gMinuit->mnparm(1,"b_index",start[1],step[1],low[1],up[1],ierr);
  gMinuit->mnparm(2,"c_index",start[2],step[2],low[2],up[2],ierr);
  gMinuit->mnparm(3,"mat_expansion",start[3],step[3],low[3],up[3],ierr);
  gMinuit->mnparm(4,"b_source",start[4],step[4],low[4],up[4],ierr);
  gMinuit->mnparm(5,"potential_0",start[5],step[5],low[5],up[5],ierr);
  gMinuit->mnparm(6,"potential_1",start[6],step[6],low[6],up[6],ierr);
  gMinuit->mnparm(7,"b_break_index",start[7],step[7],low[7],up[7],ierr);
  gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);
  gMinuit->FixParameter(4);
  gMinuit->FixParameter(5);
  gMinuit->FixParameter(6);
//  gMinuit->FixParameter(7);

  //execution
  arglist[0] = 100000;
  gMinuit->mnexcm("MIGRAD",arglist,1,ierr);

  //get results
  double fitp[8],fitperr[8];
  cout<<"Fitted results from TMinuit:"<<endl;
  for(int i=0; i<8; i++)
  {
    gMinuit->GetParameter(i,fitp[i],fitperr[i]);
    cout<<"x"<<i<<": "<<fitp[i]<<" +- "<<fitperr[i]<<endl;
  }

  //save results into histograms
  SetTrueSpectra(fitp);
  SetMatrix(fitp);
  DoFolding();
}

//get Chi2
void CalcChi2()
{
  double Delta, Sigma, dChi2;
  int StartBin = hbco->FindBin(1/FITB);
  int EndBin = hbco->FindBin(1/FITA);
  double Chi2 = 0;
  double DOF = 0;
  for(int i=StartBin; i<=EndBin; i++)
  {
    Delta = hbcf->GetBinContent(i) - hbco->GetBinContent(i);
    Sigma = hbco->GetBinError(i);
    dChi2 = Delta*Delta/Sigma/Sigma;
    hsigma->SetBinContent(i,sqrt(dChi2));
    Chi2 += dChi2;
    DOF += 1;
  }
  DOF = DOF-FIXDOF;
  cout<<"DOF = "<<DOF<<endl;
  cout<<"Chi2 = "<<Chi2<<endl;
  cout<<"Chi2/DOF = "<<Chi2/DOF<<endl;
}

//save results
void SaveResults()
{
  TFile *fw = new TFile("result.root","RECREATE");
  fw->cd();

  hbo->Write();
  hco->Write();
  hbco->Write();
  hbt->Write();
  hct->Write();
  hbct->Write();
  hbf->Write();
  hcf->Write();
  hbcf->Write();
  hmat->Write();
  gr->Write();
  fb->Write();
  fc->Write();
  heffbc->Write();
  hsigma->Write();

  fw->Write();
  fw->Close();
}

/////////////////////////////////////
//main function
/////////////////////////////////////

int main(int argc, char **argv)
{
  DefHist();
  DoMinimization();
  CalcChi2();
  //Test();
  SaveResults();
  return 0;
}
