#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

class AnalysisM
{
  public:

    //Default constructor and destructor
    AnalysisM();
    ~AnalysisM();

    //Read cut values from input document
    void ReadCuts(char *InputFile);

    //Generate efficiency and purity corrections 
    void GenerateEffCorr();

    //Calculate B to C ratio
    void CalcBCRatio();

    //Save results into root file
    void SaveResults(char *OutputFile);

    //Study selection cuts
    void StudyCuts();

    //Charge ID plots
    void PlotChargeID();

    //Calculation trigger efficiency ratio
    void CalcTrigEffR();

  private:

    ////////////////////////////////
    //Functions
    ////////////////////////////////

    //Read data from TTree
    void ReadData();

    //Define histograms
    void DefineHistograms();

    //Check whether in range
    bool IsInRange(double X, double XLow, double XUp);

    //Find which bin in rigidity
    int GetIBinRig(double Rigidity);

    //Define variables and cuts for each event
    void DefineVariables();
    void DefineCuts();

    //Calculate ratio and product values and errors of histograms
    void CalcHistRatio(TH1D *ha, TH1D *hb, TH1D *hc, bool HasError);
    void CalcHistProduct(TH1D *ha, TH1D *hb, TH1D *hc, bool HasError);

    //Get Ek/A 
    double GetEkOverAFromRigidity(double Rigidity, int Particle);
    double GetEkOverAFromBeta(double Beta);

    //Get exposure time
    double GetExpT(double Energy, int Particle);

    //Check whether above cutoff for Beta of all isotopes
    bool IsBetaAboveCutoff(double Beta);

    //Define all kinds of corrections
    void DefineCorrections();

    //For GenerateEffCorr
    void LoopEffCorr();
    void DoEffCorr();
    void DoPurityTest();

    //For CalcBCRatio
    void LoopBCRatio();
    void DoBCRatio();

    ////////////////////////////////
    //Variables
    ////////////////////////////////

    //Constants
    static const double MB10 = 9.32699;
    static const double MB11 = 10.2551;
    static const double MC12 = 11.1779;
    static const double PORTION = 1.0;
    static const double B11FRAC = 0.7;

    //TChain
    TChain *chain;

    //Variables from TTree
    Int_t RunN;
    Int_t EventN;
    Float_t Cutoff;
    Float_t Chi2T;
    Float_t Chi2C;
    Float_t TOFCl;
    Float_t Beta;
    Short_t ZTOF;
    Float_t ZProb;
    Float_t qtof[4];
    Float_t Rig[4];
    Float_t Chisq[4];
    Float_t RigErr;
    Float_t QTr1;
    Float_t QTr2;
    Float_t QTr9;
    Float_t QTrInner;
    Short_t TrYPat;
    Float_t TrCl;
    Float_t TrCl9;
    Short_t NHit9;
    Float_t YRes9;
    Float_t TrInAsyE;
    Short_t NSecTr;
    Short_t TrGeoPat;
    Float_t QRICH;
    Float_t BetaRICH;
    Float_t QTRD;
    Short_t PhyT;
    Short_t NACC;

    //Cut values to be read from input file
    double CutValue[8];

    //Variables to be defined
    int IBinRig;
    double EFromRigB10, EFromRigB11, EFromRigC12;
    double EFromRig9B10, EFromRig9B11, EFromRig9C12;
    double EFromBetaR;
    double BetaR, BRAgr;
    int YL[9], InYGood, InNHitsY, GeoL[9];
    double TOFAsyUL;
    double AbsRigL1, AbsRigIn, AbsRig9;
    double value;
    char name[100];

    //Cuts to be defined 
    bool TOFCut, TrCut, TrCut9, SelCut, SelCut9, TempCut;
    bool ztof3, ztof4, ztof5, ztof6, ztof7, ztof8;
    bool ztrin3, ztrin4, ztrin5, ztrin6, ztrin7, ztrin8;
    bool ztrtop5, ztrtop6;
    bool ztrlow5, ztrlow6;
    bool qtrtop5, qtrtop6, qtrtop7, qtrtop8;
    bool qtrd5, qtrd6, qtrd7, qtrd8;

    //Exposure time
    static double ExpT_B10[18];
    static double ExpT_B11[18];
    static double ExpT_C12[18];

    ////////////////////////////////
    //Histograms
    ////////////////////////////////

    //Binning of Rigidity - for effciency correction
    static const int NBinRig = 21;
    static double XBinRig[NBinRig+1];

    //Unfolding correction factor
    static double UnfoldingCorr[NBinRig];

    //Binning of Rigidity - for unfolding
    static const int NBinRigUn = 21;
    static double XBinRigUn[NBinRigUn+1];

    //Binning of Ek/A Based on Rigidity
    static const int NBinE = 21;
    static double XBinE[NBinE+1]; 

    //Binning of Ek/A Based on RICH Beta 
    static const int NBinERICH = 2;
    static double XBinERICH[NBinERICH+1]; 

    //Binning of charge for template fits
    static const int NBinQFit = 100;
    static const double XBinQFitLow = 0;
    static const double XBinQFitUp = 10;

    //Rigidity plots for binning
    TH2D *hhrerr; 
    TF1 *frerr;

    //Selection efficiency corrections
    TH1D *he[4], *hetr[4], *hetof[4], *hep1[4], *hep9[4], *hepr[4]; //B, C
    TH1D *her1, *her2, *hertr, *hertof, *herp1, *herp9, *herpr; //B/C
    TH1D *her1s; //Smoothed her1
    double CorrSelEff, CorrQTrL1, CorrL9Pickup, CorrRICH, CorrUnfolding;

    //Spectra on L1 and templates on L2
    TH1D *hspec1[2][NBinRig], *htmp[4];

    //Purity for B and C
    RooPlot *qframe[2][NBinRig];
    TH1D *hpb, *hpc;

    //Ek/A spectra of B, C, B/C 
    TH1D *hebrig, *hecrig, *hebcrig; //L1+Inner
    TH1D *hebrig9, *hecrig9, *hebcrig9; //L1+Inner+L9
    TH1D *hebbetar, *hecbetar, *hebcbetar; //RICH

    //Study selection cuts
    TH1D *hsg[10], *hsb[10]; //10 plots for good and bad sample
    TH1D *hsbeta[2], *hsrig[2]; //Beta and Rigidity - before and after cuts

    //Charge ID plots
    TH1D *hqtrinall, *hqtrin5, *hqtrin6;

    //Trigger efficiency
    TH1D *hetrig[4], *hetrigr;

    //Unfolding spectrum
    TH1D *hunb, *hunc; //L1+Inner
    TH1D *hunb9, *hunc9; //L1+Inner+L9
};
