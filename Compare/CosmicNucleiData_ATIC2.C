{

/********************************************************************/

// ATIC-2: A. D. Panov, ... ICRC07 (2007)

/********************************************************************/


// Energy GeV/n
Double_t  A2Energy[5] = {19.9,38.3,74.3,149.,307.};
Double_t eA2Energy[5] = {  0.,  0.,  0.,  0.,  0.};


// B/C
Double_t   A2RatioBC[5] = {0.180,0.169,0.119,0.156,0.064};
Double_t  eA2RatioBC[5] = {0.011,0.015,0.029,0.053,0.063};
TGraphErrors* ATIC2_BC = new TGraphErrors(4,A2Energy,A2RatioBC,eA2Energy,eA2RatioBC);


// N/O
Double_t   A2RatioNO[5] = {0.219,0.199,0.184,0.172,0.144};
Double_t  eA2RatioNO[5] = {0.010,0.013,0.024,0.039,0.068};
TGraphErrors* ATIC2_NO = new TGraphErrors(4,A2Energy,A2RatioNO,eA2Energy,eA2RatioNO);


// C/O
Double_t   A2RatioCO[5] = {1.020,1.087,0.933,0.934,1.022};
Double_t  eA2RatioCO[5] = {0.026,0.043,0.060,0.105,0.227};
TGraphErrors* ATIC2_CO = new TGraphErrors(4,A2Energy,A2RatioCO,eA2Energy,eA2RatioCO);


cout << "Loading ATIC-2 Nuclei data ...\n";
cout << "TGraphErrors are named ATIC2_xx where xx can be:\n";
cout << "BC, NO, CO\n";

}
