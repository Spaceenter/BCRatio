#include "AnalysisM.h"

using namespace std;
using namespace RooFit;

/////////////////////////////////////////////////////////////////////

//Binning of Rigidity - for efficiency correction
double AnalysisM::XBinRig[NBinRig+1] = {
  1.259, 1.718, 2.339, 3.098, 4.066, 5.326, 6.885, 8.915, 
  11.535, 14.955, 19.614, 25.872, 34.641, 47.399, 66.231, 95.110, 
  139.226, 210.249, 333.828, 577.082, 1341.318, 3271.122
};

//Binning of Rigidity - for unfolding 
double AnalysisM::XBinRigUn[NBinRigUn+1] = {
  1.259, 1.718, 2.339, 3.098, 4.066, 5.326, 6.885, 8.915, 11.535, 14.955, 19.614, 25.872, 34.641, 47.399, 66.231, 95.110, 139.226, 210.249, 333.828, 577.082, 1341.318, 3271.122
};

//Binning of Ek/A Based on Rigidity - in ICRC publication, even
double AnalysisM::XBinE[NBinE+1] = {
  0.192765, 0.335611, 0.563700, 0.875894, 
  1.304780, 1.889560, 2.634640, 3.622370, 4.910960, 6.603990, 8.919660, 
  12.038200, 16.413800, 22.786200, 32.196900, 46.632800, 68.687600, 
  104.197000, 165.985000, 287.611000, 669.728000, 1634.630000
}; 

/*
//Binning of Ek/A Based on Rigidity - my original one, odd, and works better
double AnalysisM::XBinE[NBinE+1] = {
  0.14791, 0.253927, 0.439944, 0.703865, 
  1.07782, 1.57625, 2.23265, 3.101, 4.21492, 5.70801, 7.6715, 
  10.331, 14.0708, 19.2986, 26.9574, 38.6047, 56.4995, 84.322, 
  130.658, 214.797, 408.689, 1634.63
};
*/

//Binning of Ek/A Based on RICH Beta
double AnalysisM::XBinERICH[NBinERICH+1] = {
  2.634640, 3.622370, 4.910960
};

//Exposure time
double AnalysisM::ExpT_B10[18] = {
  0.0750171, 0.0925195, 0.112976, 0.137171, 0.166206, 0.201828, 0.248216, 0.31489, 0.382437, 0.446007, 0.510877, 0.581642, 0.665007, 0.784206, 0.906214, 0.977915, 0.999779, 1
};
double AnalysisM::ExpT_B11[18] = {
  0.0867169, 0.105735, 0.128031, 0.154411, 0.186378, 0.22651, 0.284269, 0.352255, 0.415031, 0.477734, 0.543948, 0.618576, 0.712156, 0.8368, 0.949825, 0.996949, 0.999997, 1
};
double AnalysisM::ExpT_C12[18] = {
  0.0749522, 0.0924521, 0.112902, 0.137095, 0.166126, 0.201742, 0.248108, 0.31477, 0.382346, 0.44593, 0.510812, 0.581579, 0.664947, 0.784141, 0.90616, 0.9779, 0.999778, 1
};

double AnalysisM::UnfoldingCorr[NBinRig] = { // using 1.05 as expansion coefficient
  0.933972, 1.00823, 1.00308, 0.999412, 0.999272, 1.00137, 1.00218, 0.99823, 0.9952, 0.992551, 0.994422, 0.99004, 0.982354, 0.991582, 0.972845, 0.954405, 0.86847, 0.851341, 0.746588, 0.747967, 0.747967 
}; // last one is added manually

/////////////////////////////////////////////////////////////////////

AnalysisM::AnalysisM()
{
  chain = new TChain("qtree");
  ReadData();
  DefineHistograms();
}

/////////////////////////////////////////////////////////////////////

AnalysisM::~AnalysisM()
{
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::ReadCuts(char *InputFile)
{
  fstream File;
  string InputCutValue;
  File.open(InputFile,fstream::in);
  for(int i=0;i<9;i++)
  {
    getline(File,InputCutValue);
    CutValue[i] = atof(InputCutValue.c_str());
  }
  File.close();
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::GenerateEffCorr()
{
  LoopEffCorr();
  DoEffCorr();
//  DoPurityTest();
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::CalcBCRatio()
{
  LoopBCRatio();
  DoBCRatio();
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::SaveResults(char *OutputFile)
{
  //Open file
  TFile *f = new TFile(OutputFile,"RECREATE");
  f->cd();

  //Rigidity plots for binning
  hhrerr->Write();
  frerr->Write();

  //Selection efficiency corrections
  for(int i=0;i<4;i++)
  {
    he[i]->Write();
    hetr[i]->Write();
    hetof[i]->Write();
    hep1[i]->Write();
    hep9[i]->Write();
    hepr[i]->Write();
  }
  her1->Write();
  her2->Write();
  hertr->Write();
  hertof->Write();
  herp1->Write();
  herp9->Write();
  herpr->Write();
  her1s->Write();

  //Spectra on L1 and templates on L2
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<NBinRig;j++)
    {
      hspec1[i][j]->Write();
    }
  }
  for(int i=0;i<4;i++)
  {
    htmp[i]->Write();
  }

  /*
  //Purity for B and C
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<NBinRig;j++)
    {
      qframe[i][j]->Write();
    }
  }
  hpb->Write();
  hpc->Write();
  */

  //Ek/A spectra of B, C, B/C
  hebrig->Write();
  hecrig->Write();
  hebcrig->Write();
  hebrig9->Write();
  hecrig9->Write();
  hebcrig9->Write();
  hebbetar->Write();
  hecbetar->Write();
  hebcbetar->Write();

  //Unfolding spectrum
  hunb->Write();
  hunc->Write();
  hunb9->Write();
  hunc9->Write();

  //Close file
  f->Write();
  f->Close();
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::ReadData()
{
  chain->Add("/afs/cern.ch/work/w/weisun/public/data/th.zg2.l9.more/*.root");

  chain->SetBranchAddress("RunN",&RunN);
  chain->SetBranchAddress("EventN",&EventN);
  chain->SetBranchAddress("Cutoff",&Cutoff);
  chain->SetBranchAddress("Chi2T",&Chi2T);
  chain->SetBranchAddress("Chi2C",&Chi2C);
  chain->SetBranchAddress("TOFCl",&TOFCl);
  chain->SetBranchAddress("Beta",&Beta);
  chain->SetBranchAddress("ZTOF",&ZTOF);
  chain->SetBranchAddress("ZProb",&ZProb);
  chain->SetBranchAddress("qtof",&qtof);
  chain->SetBranchAddress("Rig",&Rig);
  chain->SetBranchAddress("Chisq",&Chisq);
  chain->SetBranchAddress("RigErr",&RigErr);
  chain->SetBranchAddress("QTr1",&QTr1);
  chain->SetBranchAddress("QTr2",&QTr2);
  chain->SetBranchAddress("QTr9",&QTr9);
  chain->SetBranchAddress("QTrInner",&QTrInner);
  chain->SetBranchAddress("TrYPat",&TrYPat);
  chain->SetBranchAddress("TrCl",&TrCl);
  chain->SetBranchAddress("TrCl9",&TrCl9);
  chain->SetBranchAddress("NHit9",&NHit9);
  chain->SetBranchAddress("YRes9",&YRes9);
  chain->SetBranchAddress("TrInAsyE",&TrInAsyE);
  chain->SetBranchAddress("NSecTr",&NSecTr);
  chain->SetBranchAddress("TrGeoPat",&TrGeoPat);
  chain->SetBranchAddress("QRICH",&QRICH);
  chain->SetBranchAddress("BetaRICH",&BetaRICH);
  chain->SetBranchAddress("QTRD",&QTRD);
  //chain->SetBranchAddress("PhyT",&PhyT);
  chain->SetBranchAddress("NACC",&NACC);
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DefineHistograms()
{
  cout<<"AnalysisM::DefineHistograms() -- Started!"<<endl;

  //Rigidity plots for binning
  hhrerr = new TH2D("hhrerr","hhrerr",32,0,3.2,120,0,1.2);
  frerr = new TF1("frerr","pol6",0,3.2);

  //Selection efficiency corrections
  for(int i=0;i<4;i++)
  {
    sprintf(name,"he_%d",i);
    he[i] = new TH1D(name,name,NBinRig,XBinRig);
    sprintf(name,"hetr_%d",i);
    hetr[i] = new TH1D(name,name,NBinRig,XBinRig);
    sprintf(name,"hetof_%d",i);
    hetof[i] = new TH1D(name,name,NBinRig,XBinRig);
    sprintf(name,"hep1_%d",i);
    hep1[i] = new TH1D(name,name,NBinRig,XBinRig);
    sprintf(name,"hep9_%d",i);
    hep9[i] = new TH1D(name,name,NBinRig,XBinRig);
    sprintf(name,"hepr_%d",i);
    hepr[i] = new TH1D(name,name,NBinRig,XBinRig);
  }
  her1 = new TH1D("her1","her1",NBinRig,XBinRig);
  her2 = new TH1D("her2","her2",NBinRig,XBinRig);
  hertr = new TH1D("hertr","hertr",NBinRig,XBinRig);
  hertof = new TH1D("hertof","hertof",NBinRig,XBinRig);
  herp1 = new TH1D("herp1","herp1",NBinRig,XBinRig);
  herp9 = new TH1D("herp9","herp9",NBinRig,XBinRig);
  herpr = new TH1D("herpr","herpr",NBinRig,XBinRig);

  //Spectra on L1 and templates on L2
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<NBinRig;j++)
    {
      sprintf(name,"hspec1_%d_%d",i+5,j);
      hspec1[i][j] = new TH1D(name,name,NBinQFit,XBinQFitLow,XBinQFitUp);
    }
  }
  for(int i=0;i<4;i++)
  {
    sprintf(name,"htmp_%d",i+5);
    htmp[i] = new TH1D(name,name,NBinQFit,XBinQFitLow,XBinQFitUp);
  }

  /*
  //Purity for B and C
  hpb = new TH1D("hpb","hpb",NBinRig,XBinRig);
  hpc = new TH1D("hpc","hpc",NBinRig,XBinRig);
  */

  //Ek/A spectra of B, C, B/C
  hebrig = new TH1D("hebrig","hebrig",NBinE,XBinE);
  hecrig = new TH1D("hecrig","hecrig",NBinE,XBinE);
  hebcrig = new TH1D("hebcrig","hebcrig",NBinE,XBinE);
  hebrig9 = new TH1D("hebrig9","hebrig9",NBinE,XBinE);
  hecrig9 = new TH1D("hecrig9","hecrig9",NBinE,XBinE);
  hebcrig9 = new TH1D("hebcrig9","hebcrig9",NBinE,XBinE);
  hebbetar = new TH1D("hebbetar","hebbetar",NBinERICH,XBinERICH);
  hecbetar = new TH1D("hecbetar","hecbetar",NBinERICH,XBinERICH);
  hebcbetar = new TH1D("hebcbetar","hebcbetar",NBinERICH,XBinERICH);

  //Unfolding spectrum
  double XBinRigUnInv[2*NBinRigUn+1];
  for(int i=0; i<NBinRigUn; i++) XBinRigUnInv[i] = -1.0/XBinRigUn[i];
  XBinRigUnInv[NBinRigUn] = 0;
  for(int i=0; i<NBinRigUn; i++) XBinRigUnInv[NBinRigUn+1+i] = -XBinRigUnInv[NBinRigUn-1-i];
  hunb = new TH1D("hunb","hunb",2*NBinRigUn,XBinRigUnInv);
  hunc = new TH1D("hunc","hunc",2*NBinRigUn,XBinRigUnInv);
  hunb9 = new TH1D("hunb9","hunb9",2*NBinRigUn,XBinRigUnInv);
  hunc9 = new TH1D("hunc9","hunc9",2*NBinRigUn,XBinRigUnInv);

  cout<<"AnalysisM::DefineHistograms() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

bool AnalysisM::IsInRange(double X, double XLow, double XUp)
{
  if(X>=XLow && X<XUp) return true;
  else return false;
}

/////////////////////////////////////////////////////////////////////

int AnalysisM::GetIBinRig(double Rigidity)
{
  if(Rigidity<0) Rigidity = 0-Rigidity;
  if(Rigidity<XBinRig[1]) return 1;
  if(Rigidity>XBinRig[NBinRig-1]) return NBinRig;
  for(int i=1; i<NBinRig-1;i++)
  {
    if(Rigidity>XBinRig[i] && Rigidity<XBinRig[i+1]) return i+1;
  }
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DefineVariables()
{
  //Bin of rigidity
  IBinRig = GetIBinRig(Rig[3]);

  //Ek/A
  EFromRigB10 = GetEkOverAFromRigidity(Rig[3],1);
  EFromRigB11 = GetEkOverAFromRigidity(Rig[3],2);
  EFromRigC12 = GetEkOverAFromRigidity(Rig[3],3);
  EFromRig9B10 = GetEkOverAFromRigidity(Rig[0],1);
  EFromRig9B11 = GetEkOverAFromRigidity(Rig[0],2);
  EFromRig9C12 = GetEkOverAFromRigidity(Rig[0],3);
  EFromBetaR = GetEkOverAFromBeta(BetaRICH);

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
  InYGood=0;
  if(YL[1] && (YL[2]||YL[3]) && (YL[4]||YL[5]) && (YL[6]||YL[7])) InYGood=1; //Inner

  //TrTrack geometry pattern
  for(int i=0;i<9;i++)
  {
    if((TrGeoPat&(1<<i))==(1<<i))
    {
      GeoL[i]=1;
    }
    else GeoL[i]=0;
  }

  //Variables used for cuts
  BetaR = 1.0/sqrt(1+4*0.939*0.939/Rig[3]/Rig[3]);
  BRAgr = TMath::Abs(Beta-BetaR)/(Beta+BetaR);
  TOFAsyUL = TMath::Abs(qtof[0]+qtof[1]-qtof[2]-qtof[3])/(qtof[0]+qtof[1]+qtof[2]+qtof[3]);

  //Absolute value of L1+Inner rigidity
  AbsRigL1 = fabs(Rig[3]);
  AbsRigIn = fabs(Rig[2]);
  AbsRig9 = fabs(Rig[0]);
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DefineCuts()
{
  //TOF
  //TOFCut = true;
  TOFCut = Chi2T<CutValue[0] && Chi2C<CutValue[1] && Beta>0.4 && Beta<1.2;
  TOFCut &= ZProb>CutValue[2] && TOFCl>CutValue[3];

  //Tr
  //TrCut = true;
  TrCut = Chisq[3]<CutValue[4] && BRAgr<CutValue[5];
  TrCut &= TrCl>CutValue[6] && TrInAsyE<CutValue[7] && NSecTr==0;
  TrCut &= InYGood && InNHitsY>=5;
  TrCut9 = Chisq[0]<CutValue[4] && BRAgr<CutValue[5];
  TrCut9 &= TrCl>CutValue[6] && TrInAsyE<CutValue[7] && NSecTr==0;
  TrCut9 &= InYGood && InNHitsY>=5;

  //SelCut and TempCut
  SelCut = TOFCut && TrCut;
  SelCut9 = TOFCut && TrCut9;
  TempCut = TOFCut && TrCut;
  TempCut &= NACC<=2 && ZProb>0.95 && Chi2T<8 && Chi2C<8;
  TempCut &= TOFAsyUL<0.1;

  //Charge
  ztof3 = ZTOF==3;
  ztof4 = ZTOF==4;
  ztof5 = ZTOF==5;
  ztof6 = ZTOF==6;
  ztof7 = ZTOF==7;
  ztof8 = ZTOF==8;
  ztrin3 = TMath::Abs(QTrInner-3)<0.5;
  ztrin4 = TMath::Abs(QTrInner-4)<0.5;
  ztrin5 = TMath::Abs(QTrInner-5)<0.5;
  ztrin6 = TMath::Abs(QTrInner-6)<0.5;
  ztrin7 = TMath::Abs(QTrInner-7)<0.5;
  ztrin8 = TMath::Abs(QTrInner-8)<0.5;
  ztrtop5 = QTr1>4 && QTr1<5.5;
  ztrtop6 = QTr1>5 && QTr1<6.8;
  ztrlow5 = QTr9>4.3 && QTr9<5.7;
  ztrlow6 = QTr9>5.3 && QTr1<6.7;
  qtrtop5 = TMath::Abs(QTr1-5.029)<2*0.239;
  qtrtop6 = TMath::Abs(QTr1-6.039)<2*0.250;
  qtrtop7 = TMath::Abs(QTr1-7.025)<2*0.267;
  qtrtop8 = TMath::Abs(QTr1-8.036)<2*0.291;
  qtrd5 = TMath::Abs(QTRD-5.02)<1.5*0.300;
  qtrd6 = TMath::Abs(QTRD-5.96)<1.5*0.379;
  qtrd7 = TMath::Abs(QTRD-7.05)<1.5*0.582;
  qtrd8 = TMath::Abs(QTRD-8.16)<1.5*0.600;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::CalcHistRatio(TH1D *ha, TH1D *hb, TH1D *hc, bool HasError)
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

/////////////////////////////////////////////////////////////////////

void AnalysisM::CalcHistProduct(TH1D *ha, TH1D *hb, TH1D *hc, bool HasError)
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
      hc->SetBinContent(i,a*b);
      hc->SetBinError(i,a*b*sqrt(ea*ea/a/a+eb*eb/b/b));
    }
    else 
    {
      hc->SetBinContent(i,1);
      hc->SetBinError(i,1);
    }
  }
}

/////////////////////////////////////////////////////////////////////

double AnalysisM::GetEkOverAFromRigidity(double Rigidity, int Particle)
{
  if(Particle==1) return (sqrt(MB10*MB10+Rigidity*Rigidity*25)-MB10)/10.0; //B10
  if(Particle==2) return (sqrt(MB11*MB11+Rigidity*Rigidity*25)-MB11)/11.0; //B11
  if(Particle==3) return (sqrt(MC12*MC12+Rigidity*Rigidity*36)-MC12)/12.0; //C12
  return 0;
}

/////////////////////////////////////////////////////////////////////

double AnalysisM::GetEkOverAFromBeta(double Beta)
{
  if(Beta>=1) return 100000;
  return (1/sqrt(1-Beta*Beta)-1)*MC12/12.0; //Same for B10, B11, C12
}

/////////////////////////////////////////////////////////////////////

double AnalysisM::GetExpT(double Energy, int Particle)
{
  double ExpT;
  int iBinE = (int)(10*log10(Energy))+4;
  if(iBinE>=0 && iBinE<18) 
  {
    if(Particle==0) return ExpT_B10[iBinE];
    if(Particle==1) return ExpT_B11[iBinE];
    if(Particle==2) return ExpT_C12[iBinE];
  }
  if(iBinE>=18) return 1;
  if(iBinE<0) 
  {
    if(Particle==0) return ExpT_B10[0];
    if(Particle==1) return ExpT_B11[0];
    if(Particle==2) return ExpT_C12[0];
  }
}

/////////////////////////////////////////////////////////////////////

bool AnalysisM::IsBetaAboveCutoff(double Beta)
{
  double RigC = 1.3*Cutoff;
  double BetaC1 = 1/sqrt(1+MB10*MB10/25.0/RigC/RigC);
  double BetaC2 = 1/sqrt(1+MB11*MB11/25.0/RigC/RigC);

  if(Beta>BetaC1 && Beta>BetaC2) return true;
  else return false;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DefineCorrections()
{
  double XRig, high, low;
  int XBin;

  //Selection efficiency
  XRig = AbsRigL1;
  high = 2000;
  low = XBinRig[1];
  if(XRig>high) XRig = high;
  if(XRig<low) XRig = low;
  XBin = her1s->FindBin(XRig); //use smoothed her1
  CorrSelEff = her1s->GetBinContent(XBin);

  //L9 pickup efficiency
  XRig = AbsRigL1;
  high = 2000;
  low = XBinRig[1];
  if(XRig>high) XRig = high;
  if(XRig<low) XRig = low;
  XBin = herp9->FindBin(XRig); 
  CorrL9Pickup = herp9->GetBinContent(XBin);

  //Unfolding correction factor
  XRig = AbsRigL1;
  high = XBinRig[NBinRig];
  low = (XBinRig[1]+XBinRig[2])/2.0;
  if(XRig>high) XRig = high;
  if(XRig<low) XRig = low;
  XBin = her1s->FindBin(XRig);
  CorrUnfolding = UnfoldingCorr[XBin-1];
 
  //CorrQTrL1
  double QTrL1EffB = htmp[0]->Integral(41,55)/htmp[0]->Integral();
  double QTrL1EffC = htmp[1]->Integral(51,68)/htmp[1]->Integral();
  CorrQTrL1 = QTrL1EffB/QTrL1EffC;

  //CorrRICH
  CorrRICH = herpr->GetBinContent(11);
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::LoopEffCorr()
{
  cout<<"AnalysisM::LoopEffCorr() -- Started!"<<endl;

  clock_t start_time = clock();

  for(int index=0;index<chain->GetEntries()/PORTION;index++)
  {
    if(index%1000000==0) cout<<"Processed: "<<index<<"/"<<chain->GetEntries()<<endl;

    chain->GetEntry(index);

    //Define variables and cuts 
    DefineVariables();
    DefineCuts();

    //Rigidity plots for binning
    if(SelCut && YL[0] && YL[8])
    {
      if(ztof6 && ztrin6 && ztrtop6) //Carbon
      {
	hhrerr->Fill(log10(fabs(Rig[0])),RigErr/Rig[0]);
      }
    }

    //Rigidity range for efficiency corrections
    if(!IsInRange(AbsRigL1,XBinRig[0],XBinRig[NBinRig])) continue;

    //Efficiency - Total at once
    if(ztof5 && ztrin5 && ztrtop5)
    {
      he[0]->Fill(AbsRigL1);
      if(SelCut) he[1]->Fill(AbsRigL1);
    }
    if(ztof6 && ztrin6 && ztrtop6)
    {
      he[2]->Fill(AbsRigL1);
      if(SelCut) he[3]->Fill(AbsRigL1);
    }

    //Efficiency - Tracker
    if(TOFCut)
    {
      if(ztof5 && ztrtop5)
      {
	hetr[0]->Fill(AbsRigL1);
	if(TrCut) hetr[1]->Fill(AbsRigL1);
      }
      if(ztof6 && ztrtop6)
      {
	hetr[2]->Fill(AbsRigL1);
	if(TrCut) hetr[3]->Fill(AbsRigL1);
      }
    }

    //Efficiency - TOF
    if(TrCut)
    {
      if(ztrin5 && ztrtop5)
      {
	hetof[0]->Fill(AbsRigL1);
	if(TOFCut) hetof[1]->Fill(AbsRigL1);
      }
      if(ztrin6 && ztrtop6)
      {
	hetof[2]->Fill(AbsRigL1);
	if(TOFCut) hetof[3]->Fill(AbsRigL1);
      }
    }

    //Efficiency - L1 XY hit pick up - on top of Inner + L1 geometry
    if(SelCut && GeoL[0])
    {
      if(ztrin5 && ztof5 && qtrd5)
      {
	hep1[0]->Fill(AbsRigIn);
	if(YL[0] && QTr1>0) hep1[1]->Fill(AbsRigIn);
      }
      if(ztrin6 && ztof6 && qtrd6)
      {
	hep1[2]->Fill(AbsRigIn);
	if(YL[0] && QTr1>0) hep1[3]->Fill(AbsRigIn);
      }
    }

    //Efficiency - L9 XY hit pick up - on top of Inner + L1 XY
    if(SelCut && YL[0] && QTr1>0)
    {
      if(ztrin5 && ztrtop5)
      {
	hep9[0]->Fill(AbsRigL1);
	if(YL[8] && QTr9>0) hep9[1]->Fill(AbsRigL1);
      }
      if(ztrin6 && ztrtop6)
      {
	hep9[2]->Fill(AbsRigL1);
	if(YL[8] && QTr9>0) hep9[3]->Fill(AbsRigL1);
      }
    }

    //Efficiency - RICH 
    if(SelCut && YL[0] && QTr1>0)
    {
      if(ztrin5 && ztrtop5)
      {
	hepr[0]->Fill(AbsRigL1);
	if(QRICH>0) hepr[1]->Fill(AbsRigL1);
      }
      if(ztrin6 && ztrtop6)
      {
	hepr[2]->Fill(AbsRigL1);
	if(QRICH>0) hepr[3]->Fill(AbsRigL1);
      }
    }

    //For fitting purity on L1
    if(SelCut && GeoL[0])
    {
      if(ztrin5) hspec1[0][IBinRig-1]->Fill(QTr1);
      if(ztrin6) hspec1[1][IBinRig-1]->Fill(QTr1);
    }
    if(TempCut)
    {
      if(ztof5 && ztrin5) htmp[0]->Fill(QTr2);
      if(ztof6 && ztrin6) htmp[1]->Fill(QTr2);
      if(ztof7 && ztrin7) htmp[2]->Fill(QTr2);
      if(ztof8 && ztrin8) htmp[3]->Fill(QTr2);
    }
  }

  double duration_time = (clock()-start_time)/(double)CLOCKS_PER_SEC;
  cout<<duration_time<<" seconds are used for this loop."<<endl;

  cout<<"AnalysisM::LoopEffCorr() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DoEffCorr()
{
  cout<<"AnalysisM::DoEffCorr() -- Started!"<<endl;

  //Rigidity plots for binning
  hhrerr->ProfileX()->Fit("frerr");

  //Efficiency - Total at once
  CalcHistRatio(he[1],he[0],he[1],false); //B
  CalcHistRatio(he[3],he[2],he[3],false); //C
  CalcHistRatio(he[1],he[3],her1,true); //B/C
  her1s = (TH1D*)her1->Clone("her1s");
  her1s->Smooth(1);

  //Efficiency - Tracker + TOF
  CalcHistRatio(hetr[1],hetr[0],hetr[1],false); //B Tr
  CalcHistRatio(hetr[3],hetr[2],hetr[3],false); //C Tr
  CalcHistRatio(hetr[1],hetr[3],hertr,true); //B/C Tr
  CalcHistRatio(hetof[1],hetof[0],hetof[1],false); //B TOF
  CalcHistRatio(hetof[3],hetof[2],hetof[3],false); //C TOF
  CalcHistRatio(hetof[1],hetof[3],hertof,true); //B/C TOF
  CalcHistProduct(hertr,hertof,her2,true); //B/C Tr + TOF

  //L1 pickup efficiency
  CalcHistRatio(hep1[1],hep1[0],hep1[1],false); //B
  CalcHistRatio(hep1[3],hep1[2],hep1[3],false); //C
  CalcHistRatio(hep1[1],hep1[3],herp1,true); //B/C

  //L9 pickup efficiency
  CalcHistRatio(hep9[1],hep9[0],hep9[1],false); //B
  CalcHistRatio(hep9[3],hep9[2],hep9[3],false); //C
  CalcHistRatio(hep9[1],hep9[3],herp9,true); //B/C

  //RICH pickup efficiency
  CalcHistRatio(hepr[1],hepr[0],hepr[1],false); //B
  CalcHistRatio(hepr[3],hepr[2],hepr[3],false); //C
  CalcHistRatio(hepr[1],hepr[3],herpr,true); //B/C

  cout<<"AnalysisM::DoEffCorr() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DoPurityTest()
{
  cout<<"AnalysisM::DoPurityTest() -- Started!"<<endl;

  char name[100];
  RooRealVar *Charge = new RooRealVar("Charge","Charge",0,10);

  //Templates
  RooDataHist *Hist[4];
  RooHistPdf *Pdf[4];
  for(int i=0;i<4;i++)
  {
    sprintf(name,"Hist_%d",i+5);
    Hist[i] = new RooDataHist(name,name,RooArgSet(*Charge),htmp[i]);
    sprintf(name,"Pdf_%d",i+5);
    Pdf[i] = new RooHistPdf(name,name,RooArgSet(*Charge),*Hist[i],0);
  }

  //Fit
  RooDataHist *HistB[NBinRig], *HistC[NBinRig];
  RooRealVar *FracB1[NBinRig], *FracB2[NBinRig], *FracB3[NBinRig]; 
  RooRealVar *FracC1[NBinRig], *FracC2[NBinRig];
  RooAddPdf *PdfB[NBinRig], *PdfC[NBinRig];
  for(int i=0; i<NBinRig; i++)
  {
    //B fit
    sprintf(name,"HistB_%d",i);
    HistB[i] = new RooDataHist(name,name,RooArgSet(*Charge),hspec1[0][i]);
    sprintf(name,"FracB1_%d",i);
    FracB1[i] = new RooRealVar(name,name,0.93,0,1);
    sprintf(name,"FracB2_%d",i);
    FracB2[i] = new RooRealVar(name,name,0.06,0,1);
    sprintf(name,"FracB3_%d",i);
    FracB3[i] = new RooRealVar(name,name,0.005,0,1);
    sprintf(name,"PdfB_%d",i);
    PdfB[i] = new RooAddPdf(name,name,RooArgList(*Pdf[0],*Pdf[1],*Pdf[2],*Pdf[3]),RooArgList(*FracB1[i],*FracB2[i],*FracB3[i]));
    PdfB[i]->fitTo(*HistB[i],Range(4.5,8.5));

    //B plot
    qframe[0][i] = Charge->frame();
    sprintf(name,"qframeB_%d",i);
    qframe[0][i]->SetNameTitle(name,name);
    HistB[i]->plotOn(qframe[0][i]);
    PdfB[i]->plotOn(qframe[0][i]);
    PdfB[i]->plotOn(qframe[0][i],Components(RooArgSet(*Pdf[0])),LineColor(kBlack),LineStyle(2));
    PdfB[i]->plotOn(qframe[0][i],Components(RooArgSet(*Pdf[1])),LineColor(kRed),LineStyle(2));
    PdfB[i]->plotOn(qframe[0][i],Components(RooArgSet(*Pdf[2])),LineColor(kOrange),LineStyle(2));
    PdfB[i]->plotOn(qframe[0][i],Components(RooArgSet(*Pdf[3])),LineColor(kGreen),LineStyle(2));
    delete HistB[i];
    delete FracB2[i];
    delete FracB3[i];
    delete PdfB[i];

    //C fit
    sprintf(name,"HistC_%d",i);
    HistC[i] = new RooDataHist(name,name,RooArgSet(*Charge),hspec1[1][i]);
    sprintf(name,"FracC1_%d",i);
    FracC1[i] = new RooRealVar(name,name,0.98,0,1);
    sprintf(name,"FracC2_%d",i);
    FracC2[i] = new RooRealVar(name,name,0.01,0,1);
    sprintf(name,"PdfC_%d",i);
    PdfC[i] = new RooAddPdf(name,name,RooArgList(*Pdf[1],*Pdf[2],*Pdf[3]),RooArgList(*FracC1[i],*FracC2[i]));
    PdfC[i]->fitTo(*HistC[i],Range(5,9));

    //C plot
    qframe[1][i] = Charge->frame();
    sprintf(name,"qframeC_%d",i);
    qframe[1][i]->SetNameTitle(name,name);
    HistC[i]->plotOn(qframe[1][i]);
    PdfC[i]->plotOn(qframe[1][i]);
    PdfC[i]->plotOn(qframe[1][i],Components(RooArgSet(*Pdf[1])),LineColor(kBlack),LineStyle(2));
    PdfC[i]->plotOn(qframe[1][i],Components(RooArgSet(*Pdf[2])),LineColor(kRed),LineStyle(2));
    PdfC[i]->plotOn(qframe[1][i],Components(RooArgSet(*Pdf[3])),LineColor(kOrange),LineStyle(2));
    delete HistC[i];
    delete FracC2[i];
    delete PdfC[i];
  }

  //Write results into histograms
  for(int i=0; i<NBinRig; i++)
  {
    hpb->SetBinContent(i+1,FracB1[i]->getVal());
    hpb->SetBinError(i+1,FracB1[i]->getError());
    delete FracB1[i];
    hpc->SetBinContent(i+1,FracC1[i]->getVal());
    hpc->SetBinError(i+1,FracC1[i]->getError());
    delete FracC1[i];
  }

  //Clean the left things
  delete Charge;
  for(int i=0; i<4; i++) 
  {
    delete Hist[i];
    delete Pdf[i];
  }

  cout<<"AnalysisM::DoPurityTest() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::LoopBCRatio()
{
  cout<<"AnalysisM::LoopBCRatio() -- Started!"<<endl;

  clock_t start_time = clock();

  for(int index=0;index<chain->GetEntries()/PORTION;index++)
  {
    if(index%1000000==0) cout<<"Processed: "<<index<<"/"<<chain->GetEntries()<<endl;

    chain->GetEntry(index);

    //Define variables and cuts 
    DefineVariables();
    DefineCuts();

    //Selection cuts
    if((!ztrin5) && (!ztrin6)) continue;

    //Define all kinds of corrections
    DefineCorrections();

    //L1+Inner 
    if(SelCut && YL[0])
    {
      if(ztrin5 && ztrtop5) //B
      {
	if(Rig[3]>0) //Exclude negative rigidity since unfolding correction is applied
	{
	  if(IsInRange(EFromRigB10,XBinE[0],XBinE[NBinE]))
	  {
	    value=1.0-B11FRAC; 
	    value=value/CorrSelEff;
	    value=value/CorrQTrL1;
	    value=value/GetExpT(EFromRigB10,0);
	    value=value*CorrUnfolding;
	    hebrig->Fill(EFromRigB10,value);
	  }
	  if(IsInRange(EFromRigB11,XBinE[0],XBinE[NBinE]))
	  {
	    value=B11FRAC; 
	    value=value/CorrSelEff;
	    value=value/CorrQTrL1;
	    value=value/GetExpT(EFromRigB11,1);
	    value=value*CorrUnfolding;
	    hebrig->Fill(EFromRigB11,value);
	  }
	}
	if(AbsRigL1<XBinRigUn[NBinRigUn]) hunb->Fill(1.0/Rig[3]); //Unfolding
      }
      if(ztrin6 && ztrtop6) //C
      {
	if(Rig[3]>0) //Exclude negative rigidity since unfolding correction is applied
	{
	  if(IsInRange(EFromRigC12,XBinE[0],XBinE[NBinE]))
	  {
	    value=1.0;
	    value=value/GetExpT(EFromRigC12,2);
	    hecrig->Fill(EFromRigC12,value);
	  }
	}
	if(AbsRigL1<XBinRigUn[NBinRigUn]) hunc->Fill(1.0/Rig[3]); //Unfolding
      }
    }

    //L1+Inner+L9
//    if(SelCut9 && YL[0] && YL[8])
//    if(SelCut9 && YL[0] && YL[8] && (Rig[3]/Rig[0]<0 || Rig[3]/Rig[0]>0.8))
    if(SelCut9 && YL[0] && YL[8] && YRes9<0.0004 && NHit9<30 && TrCl9>0.3)
    {
      if(ztrin5 && ztrtop5 && ztrlow5) //B
      {
	if(Rig[0]>0) //Exclude negative rigidity since unfolding correction is applied
	{
	  if(IsInRange(EFromRig9B10,XBinE[0],XBinE[NBinE]))
	  {
	    value=1.0-B11FRAC; 
	    value=value/CorrSelEff;
	    value=value/CorrL9Pickup;
	    value=value/CorrQTrL1;
	    value=value/GetExpT(EFromRig9B10,0);
	    //	  value=value*CorrUnfolding; // to be done
	    hebrig9->Fill(EFromRig9B10,value);
	  }
	  if(IsInRange(EFromRig9B11,XBinE[0],XBinE[NBinE]))
	  {
	    value=B11FRAC; 
	    value=value/CorrSelEff;
	    value=value/CorrL9Pickup;
	    value=value/CorrQTrL1;
	    value=value/GetExpT(EFromRig9B11,1);
	    //	  value=value*CorrUnfolding;
	    hebrig9->Fill(EFromRig9B11,value);
	  }
	}
	if(AbsRig9<XBinRigUn[NBinRigUn]) hunb9->Fill(1.0/Rig[0]); //Unfolding
      }
      if(ztrin6 && ztrtop6 && ztrlow6) //C
      {
	if(Rig[0]>0) //Exclude negative rigidity since unfolding correction is applied
	{
	  if(IsInRange(EFromRig9C12,XBinE[0],XBinE[NBinE]))
	  {
	    value=1.0;
	    value=value/GetExpT(EFromRig9C12,2);
	    hecrig9->Fill(EFromRig9C12,value);
	  }
	}
	if(AbsRig9<XBinRigUn[NBinRigUn]) hunc9->Fill(1.0/Rig[0]); //Unfolding
      }
    }

    //RICH+L1+Inner
    if(SelCut && YL[0] && QRICH!=0)
    {
      if(IsBetaAboveCutoff(BetaRICH) && BetaRICH<0.9999)
      {
	if(ztrin5 && ztrtop5) //B
	{
	  if(IsInRange(EFromBetaR,XBinERICH[0],XBinERICH[NBinERICH]))
	  {
	    value=1.0;
	    value=value/CorrSelEff;
	    value=value/CorrQTrL1;
	    value=value/CorrRICH;
	    value=value*CorrUnfolding;
	    hebbetar->Fill(EFromBetaR,value);
	  }
	}
	if(ztrin6 && ztrtop6) //C
	{
	  if(IsInRange(EFromBetaR,XBinERICH[0],XBinERICH[NBinERICH]))
	  {
	    value=1.0; 
	    hecbetar->Fill(EFromBetaR,value);
	  }
	}
      }
    }
  }

  double duration_time = (clock()-start_time)/(double)CLOCKS_PER_SEC;
  cout<<duration_time<<" seconds are used for this loop."<<endl;

  cout<<"AnalysisM::LoopBCRatio() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::DoBCRatio()
{
  //L1+Inner
  CalcHistRatio(hebrig,hecrig,hebcrig,false);

  //L1+Inner+L9
  CalcHistRatio(hebrig9,hecrig9,hebcrig9,false);

  //RICH
  CalcHistRatio(hebbetar,hecbetar,hebcbetar,false);
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::StudyCuts()
{
  cout<<"AnalysisM::StudyCuts() -- Started!"<<endl;

  clock_t start_time = clock();

  //Define histograms
  hsg[0] = new TH1D("hsg_0","hsg_0",100,0,100);
  hsg[1] = new TH1D("hsg_1","hsg_1",100,0,100);
  hsg[2] = new TH1D("hsg_2","hsg_2",100,0,1);
  hsg[3] = new TH1D("hsg_3","hsg_3",100,0,1);
  hsg[4] = new TH1D("hsg_4","hsg_4",100,0,100);
  hsg[5] = new TH1D("hsg_5","hsg_5",100,0,1);
  hsg[6] = new TH1D("hsg_6","hsg_6",100,0,1);
  hsg[7] = new TH1D("hsg_7","hsg_7",100,0,1);
  hsg[8] = new TH1D("hsg_8","hsg_8",2,0,2);
  hsg[9] = new TH1D("hsg_9","hsg_9",4,0,4);
  hsb[0] = new TH1D("hsb_0","hsb_0",100,0,100);
  hsb[1] = new TH1D("hsb_1","hsb_1",100,0,100);
  hsb[2] = new TH1D("hsb_2","hsb_2",100,0,1);
  hsb[3] = new TH1D("hsb_3","hsb_3",100,0,1);
  hsb[4] = new TH1D("hsb_4","hsb_4",100,0,100);
  hsb[5] = new TH1D("hsb_5","hsb_5",100,0,1);
  hsb[6] = new TH1D("hsb_6","hsb_6",100,0,1);
  hsb[7] = new TH1D("hsb_7","hsb_7",100,0,1);
  hsb[8] = new TH1D("hsb_8","hsb_8",2,0,2);
  hsb[9] = new TH1D("hsb_9","hsb_9",4,0,4);
  hsbeta[0] = new TH1D("hsbeta_0","hsbeta_0",150,0,1.5);
  hsbeta[1] = new TH1D("hsbeta_1","hsbeta_1",150,0,1.5);
  hsrig[0] = new TH1D("hsrig_0","hsrig_0",200,-1,1);
  hsrig[1] = new TH1D("hsrig_1","hsrig_1",200,-1,1);

  //Loop
  for(int index=0;index<chain->GetEntries()/PORTION;index++)
  {
    if(index%1000000==0) cout<<"Processed: "<<index<<"/"<<chain->GetEntries()<<endl;

    chain->GetEntry(index);

    //Define variables and cuts 
    DefineVariables();
    DefineCuts();

    //TOF cuts 
    if(ztrin6)
    {
      if(Beta<1.2 && Beta>0.4)
      {
	hsg[0]->Fill(Chi2T);
	hsg[1]->Fill(Chi2C);
	hsg[2]->Fill(ZProb);
	hsg[3]->Fill(TOFCl);
      }
      if(Beta>1.2)
      {
	hsb[0]->Fill(Chi2T);
	hsb[1]->Fill(Chi2C);
	hsb[2]->Fill(ZProb);
	hsb[3]->Fill(TOFCl);
      }
      hsbeta[0]->Fill(Beta);
      if(TOFCut) hsbeta[1]->Fill(Beta);
    }

    //Tracker cuts
    if(ztof6)
    {
      if(Rig[0]>0)
      {
	hsg[4]->Fill(Chisq[0]);
	hsg[5]->Fill(BRAgr);
	hsg[6]->Fill(TrCl);
	hsg[7]->Fill(TrInAsyE);
	hsg[8]->Fill(InYGood && InNHitsY>=5);
      }
      if(Rig[0]<0)
      {
	hsb[4]->Fill(Chisq[0]);
	hsb[5]->Fill(BRAgr);
	hsb[6]->Fill(TrCl);
	hsb[7]->Fill(TrInAsyE);
	hsb[8]->Fill(InYGood && InNHitsY>=5);
      }
      hsrig[0]->Fill(1/Rig[0]);
      if(TrCut) hsrig[1]->Fill(1/Rig[0]);
    }

    //Fragmentation cuts
    if(ztrin6)
    {
      if(QTr1>5.5 && QTr1<6.5) hsg[9]->Fill(NSecTr);
      if(QTr1>7 && QTr1<9) hsb[9]->Fill(NSecTr);
    }
  }

  //Save results
  TFile *f = new TFile("study_cuts.root","RECREATE");
  f->cd();
  for(int i=0; i<2; i++) 
  {
    hsbeta[i]->Write();
    hsrig[i]->Write();
  }
  for(int i=0; i<10; i++) 
  {
    hsg[i]->Write();
    hsb[i]->Write();
  }
  f->Write();
  f->Close();

  double duration_time = (clock()-start_time)/(double)CLOCKS_PER_SEC;
  cout<<duration_time<<" seconds are used for this loop."<<endl;

  cout<<"AnalysisM::StudyCuts() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::PlotChargeID()
{
  cout<<"AnalysisM::PlotChargeID() -- Started!"<<endl;

  clock_t start_time = clock();

  //Define histograms
  hqtrinall = new TH1D("hqtrinall","hqtrinall",120,2.5,8.5);
  hqtrin5 = new TH1D("hqtrin5","hqtrin5",120,2.5,8.5);
  hqtrin6 = new TH1D("hqtrin6","hqtrin6",120,2.5,8.5);

  //Loop
  for(int index=0;index<chain->GetEntries()/PORTION;index++)
  {
    if(index%1000000==0) cout<<"Processed: "<<index<<"/"<<chain->GetEntries()<<endl;

    chain->GetEntry(index);

    //Define variables and cuts 
    DefineVariables();
    DefineCuts();

    if(!TempCut) continue;

    hqtrinall->Fill(QTrInner);
    if(ztof5 && qtrd5 && qtrtop5) hqtrin5->Fill(QTrInner);
    if(ztof6 && qtrd6 && qtrtop6) hqtrin6->Fill(QTrInner);
  }

  //Save results
  TFile *f = new TFile("charge_id.root","RECREATE");
  f->cd();
  hqtrinall->Write();
  hqtrin5->Write();
  hqtrin6->Write();
  f->Write();
  f->Close();

  double duration_time = (clock()-start_time)/(double)CLOCKS_PER_SEC;
  cout<<duration_time<<" seconds are used for this loop."<<endl;

  cout<<"AnalysisM::PlotChargeID() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////

void AnalysisM::CalcTrigEffR()
{
  cout<<"AnalysisM::CalcTrigEffR() -- Started!"<<endl;

  clock_t start_time = clock();

  //Define histograms
  for(int i=0; i<4; i++)
  {
    sprintf(name,"hetrig_%d",i);
    hetrig[i] = new TH1D(name,name,NBinRig,XBinRig);
  }
  hetrigr = new TH1D("hetrigr","hetrigr",NBinRig,XBinRig);

  //Loop
  for(int index=0;index<chain->GetEntries()/PORTION;index++)
  {
    if(index%1000000==0) cout<<"Processed: "<<index<<"/"<<chain->GetEntries()<<endl;

    chain->GetEntry(index);

    //Define variables and cuts 
    DefineVariables();
    DefineCuts();

    //Selection cuts
    if(!SelCut) continue;
    if(YL[0]==0) continue; //L1 is needed for everyone
    if((!ztrin5) && (!ztrin6)) continue;

    //Additional bad time - special trigger setting runs
    if(RunN==1321198167) continue;
    if(RunN>=1307125541 && RunN<=1307218054) continue;

    //Z=5
    if(ztrin5 && qtrtop5)
    {
      if((PhyT&0x3E)!=0)
      {
	hetrig[0]->Fill(AbsRigL1);
	hetrig[1]->Fill(AbsRigL1);
      }
      else hetrig[0]->Fill(AbsRigL1,100);
    }

    //Z=5
    if(ztrin6 && qtrtop6)
    {
      if((PhyT&0x3E)!=0)
      {
	hetrig[2]->Fill(AbsRigL1);
	hetrig[3]->Fill(AbsRigL1);
      }
      else hetrig[2]->Fill(AbsRigL1,100);
    }
  }

  //Calculate the ratio
  CalcHistRatio(hetrig[1],hetrig[0],hetrig[1],false); //B
  CalcHistRatio(hetrig[3],hetrig[2],hetrig[3],false); //C
  CalcHistRatio(hetrig[1],hetrig[3],hetrigr,true); //B/C

  //Save results
  TFile *f = new TFile("trig_eff.root","RECREATE");
  f->cd();
  for(int i=0; i<4; i++)
  {
    hetrig[i]->Write();
  }
  hetrigr->Write();
  f->Write();
  f->Close();

  double duration_time = (clock()-start_time)/(double)CLOCKS_PER_SEC;
  cout<<duration_time<<" seconds are used for this loop."<<endl;

  cout<<"AnalysisM::CalcTrigEffR() -- Finished!"<<endl;
}

/////////////////////////////////////////////////////////////////////
