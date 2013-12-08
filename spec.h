#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "root_RVSP.h"
#include "amsdbc.h"
#include "amschain.h"

#include "Tofrec02_ihep.h"

#include "TkDBc.h"
#include "TrCharge.h"
#include "TrCluster.h"
#include "TrTrack.h"
#include "TrExtAlignDB.h"

#include "TrdKCluster.h"

/////////////////////////////////////////////////////////////////////

class spec : public AMSEventR
{
  public:
    spec(){}
    ~spec(){}
    void UBegin();
    bool UProcessCut();
    void UProcessFill();
    void UTerminate();

    //Self functions
    bool RTICut();
    bool DoTOFLayer();
    void DoTrackerLayer();
    void DoTrackerSecondary();
    bool DoTrackerRigidity();
    bool DoTrackerAcceptance();
    void DoRICH();
    bool RICHCut(int iRichRing);
    float GetTotalCollectedPhotoElectrons();
    void SetTRDVariable(int value);

//**********************************
//TTree variables
//**********************************

    //General
    Int_t RunN;
    Int_t EventN;
    Float_t Cutoff;

    //TOF
    Float_t Chi2T;
    Float_t Chi2C;
    Float_t TOFCl;
    Float_t Beta;
    Short_t ZTOF;
    Float_t ZProb;
    Float_t qtof[4];

    //Tracker
    Float_t Rig[4]; 
    Float_t Chisq[4];
    Float_t RigErr;
    Float_t QTr1;
    Float_t QTr2;
    Float_t QTr9;
    Float_t QTrInner;
    Short_t TrYPat;
    Float_t TrCl;
    Float_t TrInAsyE;
    Short_t NSecTr;
    Short_t TrGeoPat;

    //RICH
    Float_t QRICH;
    Float_t BetaRICH;

    //TRD
    Float_t QTRD;

    //Trigger and ACC
    Short_t PhyT;
    Short_t NACC;

//**********************************
//Other global variables
//**********************************

    //TOF
    BetaHR *pBH;
    int NLayer;
    float RMS,Prob;
    float TOFEdepL[4];

    //TRD
    TrdKCluster *pTRDKCluster;

    //Tracker
    TrTrackR *pTr;
    int fitcode;
    int fitid[4];
    float QTrL[9];
    int YL[9];

    //TTree
    TTree *qtree;
};

/////////////////////////////////////////////////////////////////////

void spec::UBegin()
{
  //Initialization for TTree
  qtree=new TTree("qtree","qtree");
  qtree->SetAutoSave(100000);

  qtree->Branch("RunN",&RunN,"RunN/I");
  qtree->Branch("EventN",&EventN,"EventN/I");
  qtree->Branch("Cutoff",&Cutoff,"Cutoff/F");

  qtree->Branch("Chi2T",&Chi2T,"Chi2T/F");
  qtree->Branch("Chi2C",&Chi2C,"Chi2C/F");
  qtree->Branch("TOFCl",&TOFCl,"TOFCl/F");
  qtree->Branch("Beta",&Beta,"Beta/F");
  qtree->Branch("ZTOF",&ZTOF,"ZTOF/S");
  qtree->Branch("ZProb",&ZProb,"ZProb/F");
  qtree->Branch("qtof",&qtof,"qtof[4]/F");

  qtree->Branch("Rig",&Rig,"Rig[4]/F");
  qtree->Branch("Chisq",&Chisq,"Chisq[4]/F");
  qtree->Branch("RigErr",&RigErr,"RigErr/F");
  qtree->Branch("QTr1",&QTr1,"QTr1/F");
  qtree->Branch("QTr2",&QTr2,"QTr2/F");
  qtree->Branch("QTr9",&QTr9,"QTr9/F");
  qtree->Branch("QTrInner",&QTrInner,"QTrInner/F");
  qtree->Branch("TrYPat",&TrYPat,"TrYPat/S");
  qtree->Branch("TrCl",&TrCl,"TrCl/F");
  qtree->Branch("TrInAsyE",&TrInAsyE,"TrInAsyE/F");
  qtree->Branch("NSecTr",&NSecTr,"NSecTr/S");
  qtree->Branch("TrGeoPat",&TrGeoPat,"TrGeoPat/S");

  qtree->Branch("QRICH",&QRICH,"QRICH/F");
  qtree->Branch("BetaRICH",&BetaRICH,"BetaRICH/F");

  qtree->Branch("QTRD",&QTRD,"QTRD/F");

  qtree->Branch("PhyT",&PhyT,"PhyT/S");
  qtree->Branch("NACC",&NACC,"NACC/S");

  //TRDKCluster
  pTRDKCluster=0;

  //Initialization for using RICH calibration
  RichRingR::setBetaCorrection(RichRingR::fullUniformityCorrection);

  //TrCluster linearity correction
  TrClusterR::SetLinearityCorrection();

  //Latest alignment
  TkDBc::UseLatest();
}

/////////////////////////////////////////////////////////////////////

bool spec::UProcessCut()
{
  //No DAQ error
  if(nDaqEvent()==0) return false;
  for(int i=0; i<nDaqEvent(); i++)
  {
    if(pDaqEvent(i)->HasHWError()) return false;
  }
  
  //Particle selection based on BetaH
  TofRecH::ReBuild();
  if(NBetaH()==0) return false;
  int index;
  float QTrTemp;
  float QTrMax=0;
  for(int i=0;i<NBetaH();i++)
  {
    if(pBetaH(i)->iTrTrack()==-1) continue;
    QTrTemp=pBetaH(i)->pTrTrack()->GetQ();
    if(QTrTemp>QTrMax)
    {
      QTrMax=QTrTemp;
      index=i;
    }
  }
  if(QTrMax==0) return false;
  pBH=pBetaH(index);
  Beta=pBH->GetBeta();
  if(Beta<0) return false;
  pTr=pBH->pTrTrack();

  //Select Z range
  QTrInner = pTr->GetInnerQ(Beta);
  if(QTrInner<2.3 || QTrInner>9) return false;

  //RTI 
  if(RTICut()==false) return false; 

  //BetaH
  if(pBH->IsGoodBeta()==false) return false;
  ZTOF=pBH->GetZ(NLayer,ZProb);
  Chi2T=pBH->GetNormChi2T();
  Chi2C=pBH->GetNormChi2C();
  if(DoTOFLayer()==false) return false;

  //Tracker
  if(DoTrackerRigidity()==false) return false; //Including L1 geometry requirement
  DoTrackerSecondary();
  DoTrackerLayer();

  //RICH
  DoRICH();
  
  //TRD
  if(NTrdRawHit()==0 || fitcode<0) QTRD = -20;
  else
  {
    if(pTRDKCluster) delete pTRDKCluster;
    pTRDKCluster=new TrdKCluster(this,pTr,fitcode);
    int QTRDStatus=pTRDKCluster->CalculateTRDCharge(0,Beta); 
    if(QTRDStatus<0) QTRD = -30;
    else if(pTRDKCluster->GetQNHit()<15) QTRD = -40;
    else QTRD=pTRDKCluster->GetTRDCharge(); 
  }
    
  return true;
}

/////////////////////////////////////////////////////////////////////

void spec::UProcessFill()
{
  RunN=Run();
  EventN=Event();

  //Trigger and ACC
  PhyT=pLevel1(0)->PhysBPatt;
  NACC=nAntiCluster();

  qtree->Fill();
}

/////////////////////////////////////////////////////////////////////

void spec::UTerminate()
{
  delete pTRDKCluster;
}

/////////////////////////////////////////////////////////////////////

bool spec::RTICut()
{
  AMSSetupR::RTI a;
  if(GetRTI(a)!=0) return false; //Not exist in RTI table

  bool cut[6];
  cut[0]=(a.ntrig/a.nev>0.98);
  cut[1]=(a.npart/a.ntrig>0.07/1600*a.ntrig && a.npart/a.ntrig<0.25);
  cut[2]=(a.lf>0.5);
  cut[3]=(a.zenith<40);
  cut[4]=(a.nerr>=0 && a.nerr/a.nev<0.1);
  cut[5]=(a.npart>0 && a.nev<1800);
  bool tcut=(cut[0] && cut[1] && cut[2] && cut[3] && cut[4] && cut[5]);
  if(tcut==false) return false;

  if(fabs(pTr->GetRigidity())<1.2*a.cf[3][1]) return false; //1.2*Cutoff (Positive 40 degree)

  Cutoff = a.cf[3][1];

  return true;
}

/////////////////////////////////////////////////////////////////////

bool spec::DoTOFLayer()
{
  //At least 3 bad path length, efficiency ~ 96%
  int NGPL=0;
  for(int i=0;i<4;i++) if(pBH->IsGoodQPathL(i)) NGPL++;
  if(NGPL<3) return false;

  float TOFESum=0;
  float TOFEOnTrack=0;
  for(int i=0;i<nTofClusterH();i++) TOFESum += pTofClusterH(i)->GetEdep();
  for(int i=0;i<4;i++)
  {
    TOFEOnTrack += pBH->GetEdepL(i);
    qtof[i] = pBH->GetQL(i);
  }

  TOFCl=TOFEOnTrack/TOFESum;

  return true;
}

/////////////////////////////////////////////////////////////////////

void spec::DoTrackerLayer()
{
  //QTrL and TrYPat
  TrYPat=0;
  for(int i=0;i<9;i++)
  {
    QTrL[i]=pTr->GetLayerJQ(i+1,Beta,fitcode,11.1779/6.0);
    YL[i]=0;
    if(QTrL[i])
    {
      if(pTr->GetHitLJ(i+1)->GetYCluster()) YL[i]=1;
      if((pTr->GetHitLJ(i+1)->GetQStatus()&0x1001FD)!=0) QTrL[i]=-1;
      if(pTr->GetHitLJ(i+1)->GetXCluster()==0) QTrL[i]=-2;
    }
    if(YL[i]==1) TrYPat |= (1<<i);
  }
  QTr1=QTrL[0];
  QTr2=QTrL[1];
  QTr9=QTrL[8];
  QTrInner=pTr->GetInnerQ(Beta,fitcode,11.1779/6.0);

  //Variables
  TrRecHitR *pTrHit;
  int TrLayer;
  float dx,dy;
  AMSPoint TrCooTemp;
  float TrEdep;
  float TrEdepSum[9]; //Within a circle of 2cm radius
  float TrEdepOnTrack[9];
  AMSPoint TrCoo[9];
  float QTr[9];
  float max,min;

  for(int i=0;i<9;i++) 
  {
    QTr[i]=pTr->GetLayerJQ(i+1,Beta);
    TrEdepSum[i]=0;
    TrEdepOnTrack[i]=0;
  }

  //Coordinate and Edep on track
  for(int i=0;i<pTr->NTrRecHit();i++)
  {
    pTrHit=pTr->pTrRecHit(i);
    TrLayer=pTrHit->GetLayerJ();
    TrEdep=pTrHit->GetSignalCombination(2,TrClusterR::DefaultEdepCorrOpt);
    if(TrEdep==0) QTr[TrLayer-1]=-2;
    if(QTr[TrLayer-1]<=0) continue;
    TrEdepOnTrack[TrLayer-1]=TrEdep;
    TrCoo[TrLayer-1]=pTrHit->GetGlobalCoordinateA();
  }

  //Total Edep within 2cm radius circle 
  for(int i=0;i<NTrRecHit();i++)
  {
    pTrHit=pTrRecHit(i);
    TrLayer=pTrHit->GetLayerJ();
    if(QTr[TrLayer-1]<=0) continue;
    TrEdep=pTrHit->GetSignalCombination(2,TrClusterR::DefaultEdepCorrOpt);
    if(TrEdep==0) QTr[TrLayer-1]=-2;
    if(TrEdep==0) continue;
    TrCooTemp=pTrHit->GetGlobalCoordinateA();
    dx=TrCooTemp.x()-TrCoo[TrLayer-1].x();
    dy=TrCooTemp.y()-TrCoo[TrLayer-1].y();
    if(dx*dx+dy*dy>4) continue;
    TrEdepSum[TrLayer-1]=TrEdepSum[TrLayer-1]+TrEdep;
  }

  //Tracker cleanliness and layer Edep consistency
  float TotalInnerTrEdepOnTrack=0;
  float TotalInnerTrEdepSum=0;
  max=0;
  min=999999;
  for(int i=1;i<8;i++) //Inner
  {
    if(QTr[i]<=0) continue;
    TotalInnerTrEdepOnTrack=TotalInnerTrEdepOnTrack+TrEdepOnTrack[i];
    TotalInnerTrEdepSum=TotalInnerTrEdepSum+TrEdepSum[i];
    if(TrEdepOnTrack[i]>max) max=TrEdepOnTrack[i];
    if(TrEdepOnTrack[i]<min) min=TrEdepOnTrack[i];
  }
  TrCl = TotalInnerTrEdepOnTrack/TotalInnerTrEdepSum; //Inner
  TrInAsyE = (max-min)/(max+min);
}

/////////////////////////////////////////////////////////////////////

void spec::DoTrackerSecondary()
{
  NSecTr = 0;
  for(int i=0; i<nTrTrack(); i++)
  {
    if(pTr==pTrTrack(i)) continue; //Not the main track
    for(int j=0; j<nBetaH(); j++)
    {
      if(pBetaH(j)->iTrTrack()==i) //Associated with a BetaH
      {
	float RigSec = pTrTrack(i)->GetRigidity();
	int NHitsYSec = pTrTrack(i)->GetNhitsY();
	if(RigSec>0.5 && NHitsYSec>=4) //Good secondary track
	{
	  NSecTr++;
	  break;
	}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////

bool spec::DoTrackerRigidity()
{
  fitid[0]=pTr->iTrTrackPar(1,0,3); //Choutko, full span
  fitcode=fitid[0]; //Choutko has higher efficiency
  if(fitcode<0) return false;
  Chisq[0]=pTr->GetChisq(fitid[0]);
  Rig[0]=pTr->GetRigidity(fitid[0]);
  RigErr = pTr->GetErrRinv(fitid[0])*Rig[0]*Rig[0];

  //Require L1 geometry
  if(DoTrackerAcceptance()==false) return false;

  fitid[1]=pTr->iTrTrackPar(3,0,3); //Chikanian, full span
  Chisq[1]=pTr->GetChisq(fitid[1]);
  Rig[1]=pTr->GetRigidity(fitid[1]);

  fitid[2]=pTr->iTrTrackPar(1,3,3); //Choutko, inner
  Chisq[2]=pTr->GetChisq(fitid[2]);
  Rig[2]=pTr->GetRigidity(fitid[2]);

  fitid[3]=pTr->iTrTrackPar(1,5,3); //Choutko, inner + L1
  Chisq[3]=pTr->GetChisq(fitid[3]);
  Rig[3]=pTr->GetRigidity(fitid[3]);

  return true;
}

/////////////////////////////////////////////////////////////////////

bool spec::DoTrackerAcceptance()
{
  double Z[9] = {158.920,53.060,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882};
  float Edge[9][4] = {
    {-62.14,  -47.40,   62.14,   47.40},
    {-62.14,  -40.10,   62.14,   40.10},
    {-49.70,  -43.75,   49.70,   43.75},
    {-49.72,  -43.75,   49.72,   43.75},
    {-49.71,  -36.45,   49.70,   36.45},
    {-49.72,  -36.45,   49.72,   36.45},
    {-49.72,  -43.75,   49.71,   43.75},
    {-49.72,  -43.75,   49.71,   43.75},
    {-45.62,  -29.48,   45.55,   29.53}
  };

  AMSPoint Point;
  float X[9], Y[9];
  TrGeoPat = 0;
  for(int ii=0; ii<9; ii++)
  {
    Point = pTr->InterpolateLayerJ(ii+1, fitcode);
    X[ii]=Point.x();
    Y[ii]=Point.y();

    bool IsInLayer=false;
    if(X[ii]>Edge[ii][0] && X[ii]<Edge[ii][2] && Y[ii]>Edge[ii][1] && Y[ii]<Edge[ii][3])
    {
      if(ii==8) IsInLayer=true;
      else if((X[ii]*X[ii])+(Y[ii]*Y[ii])<Edge[ii][2]*Edge[ii][2]) IsInLayer=true;
    }

    //Require L1 geometry
    if(ii==0 && IsInLayer==false) return false;

    if(IsInLayer) TrGeoPat |= (1<<ii);
  }

  return true;
}

/////////////////////////////////////////////////////////////////////

void spec::DoRICH()
{
  QRICH=0;
  BetaRICH=0;
  if(nRichRing())
  {
    int iRichRing=-1;
    for(int i=0;i<nRichRing();i++)
    {
      if(pRichRing(i)->pTrTrack()==pTr) 
      {
	iRichRing=i;
	break;
      }
    }
    if(iRichRing!=-1)
    {
      if(RICHCut(iRichRing))
      {
	QRICH=sqrt(pRichRing(iRichRing)->getCharge2Estimate());
	BetaRICH=pRichRing(iRichRing)->getBeta();
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////

bool spec::RICHCut(int iRichRing)
{
  bool IsGood=pRichRing(iRichRing)->IsGood();
  int Used=pRichRing(iRichRing)->Used;
  float Prob=pRichRing(iRichRing)->Prob;
  float UDist=pRichRing(iRichRing)->UDist;
  float NPERing=pRichRing(iRichRing)->getPhotoElectrons();
  float NPETotal=GetTotalCollectedPhotoElectrons();
  int NpExp=pRichRing(iRichRing)->NpExp;

  if(!IsGood) return false;
  if(Used<3) return false;
  if(Prob<0.01) return false;
  if(UDist>20) return false;
  if(NPERing<0.5*NPETotal) return false;
  if(pRichRing(iRichRing)->PmtCorrectionsFailed()) return false;
  if(NpExp<3) return false;

  return true;
}

/////////////////////////////////////////////////////////////////////

float spec::GetTotalCollectedPhotoElectrons()
{
  float counter=0;
  for(int i=0;i<nRichHit();i++)
  {
    if(!pRichHit(i)) continue;
    if(pRichHit(i)->IsCrossed()) continue;
    counter=counter+pRichHit(i)->Npe;
  }
  return counter;
}

/////////////////////////////////////////////////////////////////////

