{

/********************************************************************/

// CREAM II: H. S. Ahn, ... arXiv:0808.1718v1 (2008)

/********************************************************************/

// Energy GeV/n
Double_t   CIIEnergyLowEdge[7] = {1.,4.,16.,63.,251.,1000.,4000.};
Double_t   CIIEnergy[6] = {1.4,5.7,23.,91.,363.,1450.};
Double_t  eCIIEnergy[6] = {0.2,0.9, 3.,14., 54., 217.};

// B/C
Double_t   CIIRatioBC[6] = {0.320,0.225,0.155,0.101,0.071,0.033};
Double_t  eCIIRatioBC[6] = {
  sqrt(pow(0.007,2.) + pow(0.030,2.)),
  sqrt(pow(0.004,2.) + pow(0.020,2.)),
  sqrt(pow(0.005,2.) + pow(0.014,2.)),
  sqrt(pow(0.011,2.) + pow(0.010,2.)),
  sqrt(pow(0.025,2.) + pow(0.010,2.)),
  sqrt(pow(0.025,2.) + pow(0.010,2.))
};
TGraphErrors* CREAMII_BC = new TGraphErrors(6,CIIEnergy,CIIRatioBC,eCIIEnergy,eCIIRatioBC);

// N/O
Double_t   CIIRatioNO[6] = {0.299,0.246,0.210,0.174,0.124,0.050};
Double_t  eCIIRatioNO[6] = {
  sqrt(pow(0.006,2.) + pow(0.030,2.)),
  sqrt(pow(0.004,2.) + pow(0.025,2.)),
  sqrt(pow(0.009,2.) + pow(0.020,2.)),
  sqrt(pow(0.026,2.) + pow(0.020,2.)),
  sqrt(pow(0.072,2.) + pow(0.010,2.)),
  sqrt(pow(0.034,2.) + pow(0.010,2.))
};
TGraphErrors* CREAMII_NO = new TGraphErrors(6,CIIEnergy,CIIRatioNO,eCIIEnergy,eCIIRatioNO);

// C/O
Double_t   CIIRatioCO[6] = {1.10,1.11,1.16,1.25,1.25,0.66};
Double_t  eCIIRatioCO[6] = {
  sqrt(pow(0.01,2.) + pow(0.1,2.)),
  sqrt(pow(0.01,2.) + pow(0.1,2.)),
  sqrt(pow(0.03,2.) + pow(0.1,2.)),
  sqrt(pow(0.10,2.) + pow(0.1,2.)),
  sqrt(pow(0.32,2.) + pow(0.1,2.)),
  sqrt(pow(0.41,2.) + pow(0.1,2.))
};
TGraphErrors* CREAMII_CO = new TGraphErrors(6,CIIEnergy,CIIRatioCO,eCIIEnergy,eCIIRatioCO);

cout << "Loading CREAMII Nuclei data ...\n";
cout << "TGraphErrors are named CREAMII_xx where xx can be:\n";
cout << "BC, NO, CO\n";

}
