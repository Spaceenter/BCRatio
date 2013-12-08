{

/********************************************************************/

// JJ Engelmann et al. 1990 - A&A 
// Franch-Danish experiment on board of the NASA HEAO-3 satellite flown in 1979 
// 5 Cerenkov counters + hodoscope by 4 tube arrays 
// 0.6 GeV/n < E < 35 GeV/n, 4 <= Z <= 28 no isotopes

/********************************************************************/

// Energy Bins GeV/n
Double_t  EnergyBins[14] = {0.55,0.70,0.91,1.11,1.40,1.82,2.35,2.96,3.79,4.89,6.42,8.60,12.0,17.8};

// Energy GeV/n
Double_t  EEnergy[14]    = { 0.62, 0.80, 1.00, 1.25, 1.60, 2.10, 2.65, 3.35, 4.30, 5.60, 7.50,10.60,16.20,35.00}; 
Double_t eEEnergy[14]    = {   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0.,  0.,   0.,   0.,   0.,   0.};

// Oxygen absolute Flux (m^2 sr s GeV/n)^-1 
Double_t  EDiffFluxO[14] = {2.695,2.205,1.675,1.504,1.091,0.753,0.496,0.335,0.220,0.132,0.073,0.0329,0.0106,0.0013};
Double_t eEDiffFluxO[14] = {0.042,0.035,0.030,0.025,0.018,0.013,0.006,0.004,0.003,0.002,0.002,0.0006,0.0004,0.0002}; 
TGraphErrors*  HEAO_O = new TGraphErrors(14,EEnergy,EDiffFluxO,eEEnergy,eEDiffFluxO);

// Other data are given as relative abundances (O = 1000)

// Be
// last bin is missing
Double_t  ERelAbundBe[13] = {  85.8, 105.5, 108.4, 126.1, 118.6, 120.3, 120.2, 121.2, 111.4, 107.2, 103.9, 87.3, 73.2}; 
Double_t eERelAbundBe[13] = {   3.8,   3.4,   3.2,   2.9,   2.1,   2.1,   2.1,   2.1,   2.0,   2.0,   1.8,  1.8,  2.1};
Double_t  EDiffFluxBe[13];
Double_t eEDiffFluxBe[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxBe[ii]  = ERelAbundBe[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxBe[ii] = EDiffFluxBe[ii]*sqrt(pow(eERelAbundBe[ii]/ERelAbundBe[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_Be = new TGraphErrors(13,EEnergy,EDiffFluxBe,eEEnergy,eEDiffFluxBe);

// B
Double_t  ERelAbundB[14]  = { 354.4, 362.3, 365.4, 364.9, 356.3, 321.0, 307.5, 294.8, 264.1, 251.7, 225.4,192.5,162.0,116.0};
Double_t eERelAbundB[14]  = {   7.3,   5.9,   5.7,   4.7,   3.6,   3.3,   3.2,   3.0,   2.8,   2.8,   2.6,  2.5,  3.0,  2.2};
Double_t  EDiffFluxB[14];
Double_t eEDiffFluxB[14];
for (int ii=0; ii<14; ii++) {
  EDiffFluxB[ii]  = ERelAbundB[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxB[ii] = EDiffFluxB[ii]*sqrt(pow(eERelAbundB[ii]/ERelAbundB[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_B = new TGraphErrors(14,EEnergy,EDiffFluxB,eEEnergy,eEDiffFluxB);

// C
Double_t  ERelAbundC[14]  = {1104.0,1108.7,1089.6,1092.1,1115.4,1073.3,1076.4,1084.7,1039.8,1049.0,1037.4,986.3,957.8,958.5};
Double_t eERelAbundC[14]  = {  12.7,  10.2,   9.6,   8.0,   6.2,   5.8,   5.9,   5.7,   5.4,   5.5,   5.4,  5.6,  7.1,  7.4};
Double_t  EDiffFluxC[14];
Double_t eEDiffFluxC[14];
for (int ii=0; ii<14; ii++) {
  EDiffFluxC[ii]  = ERelAbundC[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxC[ii] = EDiffFluxC[ii]*sqrt(pow(eERelAbundC[ii]/ERelAbundC[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_C = new TGraphErrors(14,EEnergy,EDiffFluxC,eEEnergy,eEDiffFluxC);

// N
Double_t  ERelAbundN[14]  = { 313.0, 314.7, 319.4, 300.8, 302.9, 291.6, 283.6, 273.5, 269.4, 252.3, 234.8,218.7,201.5,173.9};
Double_t eERelAbundN[14]  = {   6.8,   5.4,   5.2,   4.2,   3.2,   3.0,   3.0,   2.9,   2.8,   2.7,   2.6,  2.6,  3.2,  3.1};
Double_t  EDiffFluxN[14];
Double_t eEDiffFluxN[14];
for (int ii=0; ii<14; ii++) {
  EDiffFluxN[ii]  = ERelAbundN[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxN[ii] = EDiffFluxN[ii]*sqrt(pow(eERelAbundN[ii]/ERelAbundN[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_N = new TGraphErrors(14,EEnergy,EDiffFluxN,eEEnergy,eEDiffFluxN);

// F
Double_t  ERelAbundF[14]  = {  27.1,  24.0,  22.0,  21.5,  20.9,  19.8,  20.7,  19.0,  19.2,  18.1,  15.4, 15.5, 13.5, 10.1};
Double_t eERelAbundF[14]  = {   2.0,   1.5,   1.4,   1.0,   0.9,   0.8,   0.8,   0.8,   0.7,   0.7,   0.7,  0.7,  0.8,  0.8};
Double_t  EDiffFluxF[14];
Double_t eEDiffFluxF[14];
for (int ii=0; ii<14; ii++) {
  EDiffFluxF[ii]  = ERelAbundF[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxF[ii] = EDiffFluxF[ii]*sqrt(pow(eERelAbundF[ii]/ERelAbundF[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_F = new TGraphErrors(14,EEnergy,EDiffFluxF,eEEnergy,eEDiffFluxF);

// Ne
Double_t  ERelAbundNe[14] = { 151.9, 162.0, 167.6, 157.6, 157.3, 152.2, 158.4, 154.1, 155.8, 154.0, 154.5,152.2,153.9,144.2};
Double_t eERelAbundNe[14] = {   3.5,   2.9,   2.8,   2.2,   1.7,   1.6,   1.6,   1.5,   1.5,   1.5,   1.5,  1.5,  1.9,  2.1};
Double_t  EDiffFluxNe[14];
Double_t eEDiffFluxNe[14];
for (int ii=0; ii<14; ii++) {
  EDiffFluxNe[ii]  = ERelAbundNe[ii]*EDiffFluxO[ii]/1000.;
  eEDiffFluxNe[ii] = EDiffFluxNe[ii]*sqrt(pow(eERelAbundNe[ii]/ERelAbundNe[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}
TGraphErrors*  HEAO_Ne = new TGraphErrors(14,EEnergy,EDiffFluxNe,eEEnergy,eEDiffFluxNe);


// Energy GeV/n
Double_t  EEnergy2[13]    = { 0.80, 1.00, 1.25, 1.60, 2.10, 2.65, 3.35, 4.30, 5.60, 7.50,10.60,16.20,35.00};
Double_t eEEnergy2[13]    = {   0.,   0.,   0.,   0.,   0.,   0.,   0.,    0.,  0.,   0.,   0.,   0.,   0.};


// Na 
// first bin is missing
Double_t  ERelAbundNa[13] = { 34.7, 35.6, 34.7, 34.2, 33.0, 32.7, 28.5, 30.8, 29.9, 27.8, 26.1, 23.1, 19.3};
Double_t eERelAbundNa[13] = {  1.3,  1.3,  1.1,  0.8,  0.8,  0.8,  0.7,  0.7,  0.6,  0.6,  0.6,  0.7,  0.8};
Double_t  EDiffFluxNa[13];
Double_t eEDiffFluxNa[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxNa[ii]  = ERelAbundNa[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxNa[ii] = EDiffFluxNa[ii]*sqrt(pow(eERelAbundNa[ii]/ERelAbundNa[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Na = new TGraphErrors(13,EEnergy2,EDiffFluxNa,eEEnergy2,eEDiffFluxNa);

// Mg
// first bin is missing
Double_t  ERelAbundMg[13] = {208.7,206.7,206.5,203.4,198.4,202.4,194.9,202.9,204.8,192.0,197.3,199.6,190.5};
Double_t eERelAbundMg[13] = {  3.3,  3.1,  2.6,  2.0,  1.9,  1.9,  1.7,  1.8,  1.7,  1.6,  1.8,  2.2,  2.5};
Double_t  EDiffFluxMg[13];
Double_t eEDiffFluxMg[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxMg[ii]  = ERelAbundMg[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxMg[ii] = EDiffFluxMg[ii]*sqrt(pow(eERelAbundMg[ii]/ERelAbundMg[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Mg = new TGraphErrors(13,EEnergy2,EDiffFluxMg,eEEnergy2,eEDiffFluxMg);

// Al
// first bin is missing
Double_t  ERelAbundAl[13] = { 35.7, 34.6, 36.4, 34.1, 35.1, 35.8, 33.8, 33.8, 33.2, 31.2, 31.3, 30.8, 30.0};
Double_t eERelAbundAl[13] = {  1.4,  1.3,  1.1,  0.8,  0.8,  0.8,  0.7,  0.7,  0.7,  0.7,  0.7,  0.9,  1.0};
Double_t  EDiffFluxAl[13];
Double_t eEDiffFluxAl[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxAl[ii]  = ERelAbundAl[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxAl[ii] = EDiffFluxAl[ii]*sqrt(pow(eERelAbundAl[ii]/ERelAbundAl[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Al = new TGraphErrors(13,EEnergy2,EDiffFluxAl,eEEnergy2,eEDiffFluxAl);

// Si
// first bin is missing
Double_t  ERelAbundSi[13] = {152.6,158.4,155.4,154.6,151.7,156.0,158.5,160.3,162.9,161.2,163.4,175.3,175.7};
Double_t eERelAbundSi[13] = {  2.9,  2.8,  2.3,  1.7,  1.6,  1.7,  1.6,  1.6,  1.5,  1.5,  1.6,  2.1,  2.4};
Double_t  EDiffFluxSi[13];
Double_t eEDiffFluxSi[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxSi[ii]  = ERelAbundSi[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxSi[ii] = EDiffFluxSi[ii]*sqrt(pow(eERelAbundSi[ii]/ERelAbundSi[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Si = new TGraphErrors(13,EEnergy2,EDiffFluxSi,eEEnergy2,eEDiffFluxSi);

// P
// first bin is missing
Double_t  ERelAbundP[13]  = {  7.2,  7.7,  7.1,  6.9,  6.8,  6.9,  6.5,  5.9,  6.2,  6.1,  5.3,  5.3,  4.7};
Double_t eERelAbundP[13]  = {  0.6,  0.6,  0.5,  0.4,  0.4,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.4,  0.4};
Double_t  EDiffFluxP[13];
Double_t eEDiffFluxP[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxP[ii]  = ERelAbundP[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxP[ii] = EDiffFluxP[ii]*sqrt(pow(eERelAbundP[ii]/ERelAbundP[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_P = new TGraphErrors(13,EEnergy2,EDiffFluxP,eEEnergy2,eEDiffFluxP);

// S
Double_t  ERelAbundS[13]  = { 30.9, 34.6, 30.9, 31.0, 29.5, 31.8, 31.6, 30.1, 31.1, 31.0, 29.9, 31.2, 31.7};
Double_t eERelAbundS[13]  = {  1.3,  1.3,  1.0,  0.8,  0.7,  0.8,  0.7,  0.7,  0.7,  0.7,  0.7,  0.9,  1.0};
Double_t  EDiffFluxS[13];
Double_t eEDiffFluxS[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxS[ii]  = ERelAbundS[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxS[ii] = EDiffFluxS[ii]*sqrt(pow(eERelAbundS[ii]/ERelAbundS[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_S = new TGraphErrors(13,EEnergy2,EDiffFluxS,eEEnergy2,eEDiffFluxS);

// Cl
// first bin is missing
Double_t  ERelAbundCl[13] = {  6.4,  7.5,  7.2,  7.4,  6.9,  6.5,  6.4,  6.1,  5.8,  5.6,  5.2,  5.0,  3.9};
Double_t eERelAbundCl[13] = {  0.6,  0.6,  0.5,  0.4,  0.4,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.4,  0.3};
Double_t  EDiffFluxCl[13];
Double_t eEDiffFluxCl[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxCl[ii]  = ERelAbundCl[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxCl[ii] = EDiffFluxCl[ii]*sqrt(pow(eERelAbundCl[ii]/ERelAbundCl[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Cl = new TGraphErrors(13,EEnergy2,EDiffFluxCl,eEEnergy2,eEDiffFluxCl);

// Ar
// first bin is missing
Double_t  ERelAbundAr[13] = { 13.1, 14.8, 14.0, 13.2, 13.1, 13.5, 11.3, 11.0, 10.9, 10.0,  9.3,  9.7,  9.0};
Double_t eERelAbundAr[13] = {  0.9,  0.9,  0.7,  0.5,  0.5,  0.5,  0.4,  0.4,  0.4,  0.4,  0.4,  0.5,  0.5};
Double_t  EDiffFluxAr[13];
Double_t eEDiffFluxAr[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxAr[ii]  = ERelAbundAr[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxAr[ii] = EDiffFluxAr[ii]*sqrt(pow(eERelAbundAr[ii]/ERelAbundAr[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Ar = new TGraphErrors(13,EEnergy2,EDiffFluxAr,eEEnergy2,eEDiffFluxAr);

// K
// first bin is missing
Double_t  ERelAbundK[13]  = {  9.8, 11.3,  9.4,  9.8,  8.8,  8.9,  7.6,  8.2,  8.4,  7.5,  6.1,  6.6,  5.4};
Double_t eERelAbundK[13]  = {  0.8,  0.8,  0.6,  0.5,  0.4,  0.4,  0.4,  0.4,  0.4,  0.3,  0.3,  0.4,  0.4};
Double_t  EDiffFluxK[13];
Double_t eEDiffFluxK[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxK[ii]  = ERelAbundK[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxK[ii] = EDiffFluxK[ii]*sqrt(pow(eERelAbundK[ii]/ERelAbundK[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_K = new TGraphErrors(13,EEnergy2,EDiffFluxK,eEEnergy2,eEDiffFluxK);

// Ca
// first bin is missing
Double_t  ERelAbundCa[13] = { 22.0, 23.4, 23.4, 22.3, 19.8, 22.3, 19.6, 19.4, 19.2, 18.4, 18.1, 18.2, 18.7};
Double_t eERelAbundCa[13] = {  1.1,  1.1,  0.9,  0.7,  0.6,  0.7,  0.6,  0.6,  0.6,  0.5,  0.6,  0.7,  0.8};
Double_t  EDiffFluxCa[13];
Double_t eEDiffFluxCa[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxCa[ii]  = ERelAbundCa[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxCa[ii] = EDiffFluxCa[ii]*sqrt(pow(eERelAbundCa[ii]/ERelAbundCa[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Ca = new TGraphErrors(13,EEnergy2,EDiffFluxCa,eEEnergy2,eEDiffFluxCa);

// Sc
// first bin is missing
Double_t  ERelAbundSc[13] = {  4.7,  5.4,  5.4,  5.0,  3.8,  4.0,  4.1,  3.7,  3.4,  3.3,  3.1,  2.8,  2.9};
Double_t eERelAbundSc[13] = {  0.5,  0.5,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.2,  0.2,  0.2,  0.3,  0.3};
Double_t  EDiffFluxSc[13];
Double_t eEDiffFluxSc[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxSc[ii]  = ERelAbundSc[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxSc[ii] = EDiffFluxSc[ii]*sqrt(pow(eERelAbundSc[ii]/ERelAbundSc[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Sc = new TGraphErrors(13,EEnergy2,EDiffFluxSc,eEEnergy2,eEDiffFluxSc);

// Ti
// first bin is missing
Double_t  ERelAbundTi[13] = { 14.3, 16.4, 15.0, 14.6, 13.8, 13.6, 13.0, 11.9, 12.0, 10.3, 10.0, 10.0,  8.5};
Double_t eERelAbundTi[13] = {  0.9,  1.0,  0.8,  0.6,  0.5,  0.5,  0.5,  0.5,  0.4,  0.4,  0.4,  0.5,  0.5};
Double_t  EDiffFluxTi[13];
Double_t eEDiffFluxTi[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxTi[ii]  = ERelAbundTi[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxTi[ii] = EDiffFluxTi[ii]*sqrt(pow(eERelAbundTi[ii]/ERelAbundTi[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Ti = new TGraphErrors(13,EEnergy2,EDiffFluxTi,eEEnergy2,eEDiffFluxTi);

// V
// first bin is missing
Double_t  ERelAbundV[13]  = {  6.4,  7.0,  7.0,  7.7,  7.2,  6.9,  5.8,  6.7,  5.6,  5.6,  5.2,  5.3,  4.3};
Double_t eERelAbundV[13]  = {  0.6,  0.6,  0.5,  0.4,  0.4,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.4,  0.4};
Double_t  EDiffFluxV[13];
Double_t eEDiffFluxV[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxV[ii]  = ERelAbundV[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxV[ii] = EDiffFluxV[ii]*sqrt(pow(eERelAbundV[ii]/ERelAbundV[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_V = new TGraphErrors(13,EEnergy2,EDiffFluxV,eEEnergy2,eEDiffFluxV);

// Cr
// first bin is missing
Double_t  ERelAbundCr[13] = { 13.6, 15.7, 16.3, 14.0, 13.6, 15.4, 12.4, 13.9, 12.8, 12.1, 10.7, 11.8, 11.0};
Double_t eERelAbundCr[13] = {  0.9,  0.9,  0.8,  0.6,  0.5,  0.6,  0.5,  0.5,  0.5,  0.5,  0.4,  0.6,  0.6};
Double_t  EDiffFluxCr[13];
Double_t eEDiffFluxCr[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxCr[ii]  = ERelAbundCr[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxCr[ii] = EDiffFluxCr[ii]*sqrt(pow(eERelAbundCr[ii]/ERelAbundCr[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Cr = new TGraphErrors(13,EEnergy2,EDiffFluxCr,eEEnergy2,eEDiffFluxCr);

// Mn
// first bin is missing
Double_t  ERelAbundMn[13] = {  6.9,  8.4, 10.1,  9.0,  9.9,  9.5,  9.7,  8.7,  9.3,  9.0,  8.8,  9.0,  9.3};
Double_t eERelAbundMn[13] = {  0.6,  0.7,  0.6,  0.4,  0.5,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.5,  0.6};
Double_t  EDiffFluxMn[13];
Double_t eEDiffFluxMn[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxMn[ii]  = ERelAbundMn[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxMn[ii] = EDiffFluxMn[ii]*sqrt(pow(eERelAbundMn[ii]/ERelAbundMn[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Mn = new TGraphErrors(13,EEnergy2,EDiffFluxMn,eEEnergy2,eEDiffFluxMn);

// Fe
// first bin is missing
Double_t  ERelAbundFe[13] = { 86.7, 95.1,101.3,101.2, 97.0,105.7,105.1,108.9,112.6,107.2,110.0,125.6,135.5};
Double_t eERelAbundFe[13] = {  2.4,  2.4,  2.0,  1.5,  1.5,  1.5,  1.4,  1.4,  1.4,  1.4,  1.5,  1.9,  2.3};
Double_t  EDiffFluxFe[13];
Double_t eEDiffFluxFe[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxFe[ii]  = ERelAbundFe[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxFe[ii] = EDiffFluxFe[ii]*sqrt(pow(eERelAbundFe[ii]/ERelAbundFe[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Fe = new TGraphErrors(13,EEnergy2,EDiffFluxFe,eEEnergy2,eEDiffFluxFe);

// Co
// first bin is missing
Double_t  ERelAbundCo[13] = { 0.35, 0.57, 0.58, 0.48, 0.66, 0.57, 0.37, 0.60, 0.75, 0.58, 0.38, 0.61, 0.52};
Double_t eERelAbundCo[13] = { 0.10, 0.17, 0.11, 0.06, 0.10, 0.10, 0.07, 0.10, 0.10, 0.09, 0.06, 0.12, 0.11};
Double_t  EDiffFluxCo[13];
Double_t eEDiffFluxCo[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxCo[ii]  = ERelAbundCo[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxCo[ii] = EDiffFluxCo[ii]*sqrt(pow(eERelAbundCo[ii]/ERelAbundCo[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Co = new TGraphErrors(13,EEnergy2,EDiffFluxCo,eEEnergy2,eEDiffFluxCo);

// Ni
// first bin is missing
Double_t  ERelAbundNi[13] = {  3.5,  4.1,  4.9,  4.7,  4.7,  5.6,  5.9,  6.1,  5.7,  5.6,  7.0,  6.9,  6.7};
Double_t eERelAbundNi[13] = {  0.5,  0.5,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.4,  0.4,  0.5};
Double_t  EDiffFluxNi[13];
Double_t eEDiffFluxNi[13];
for (int ii=0; ii<13; ii++) {
  EDiffFluxNi[ii]  = ERelAbundNi[ii]*EDiffFluxO[ii+1]/1000.;
  eEDiffFluxNi[ii] = EDiffFluxNi[ii]*sqrt(pow(eERelAbundNi[ii]/ERelAbundNi[ii],2.) + pow(eEDiffFluxO[ii+1]/EDiffFluxO[ii+1],2.));
}
TGraphErrors*  HEAO_Ni = new TGraphErrors(13,EEnergy2,EDiffFluxNi,eEEnergy2,eEDiffFluxNi);

//////////////
/// RATIOS ///
//////////////

Double_t  ERatioBeC[13];
Double_t eERatioBeC[13];
Double_t  ERatioBC[14];
Double_t eERatioBC[14];
Double_t  ERatioCO[14];
Double_t eERatioCO[14];
Double_t  ERatioNO[14];
Double_t eERatioNO[14];
Double_t  ERatiosubFeFe[13];
Double_t eERatiosubFeFe[13];

for (int ii=0; ii<13; ii++) {
  ERatioBeC[ii] = EDiffFluxBe[ii]/EDiffFluxC[ii];
  eERatioBeC[ii]= ERatioBeC[ii]*sqrt(pow(eEDiffFluxBe[ii]/EDiffFluxBe[ii],2.) + pow(eEDiffFluxC[ii]/EDiffFluxC[ii],2.));
  // cout << EDiffFluxBe[ii] << " " << EDiffFluxC[ii] << " " << ERatioBeC[ii] << endl;
}

for (int ii=0; ii<14; ii++) {
  ERatioBC[ii]  = EDiffFluxB[ii]/EDiffFluxC[ii];
  eERatioBC[ii] = ERatioBC[ii]*sqrt(pow(eEDiffFluxB[ii]/EDiffFluxB[ii],2.) + pow(eEDiffFluxC[ii]/EDiffFluxC[ii],2.));
  ERatioCO[ii]  = EDiffFluxC[ii]/EDiffFluxO[ii];
  eERatioCO[ii] = ERatioCO[ii]*sqrt(pow(eEDiffFluxC[ii]/EDiffFluxC[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
  ERatioNO[ii]  = EDiffFluxN[ii]/EDiffFluxO[ii];
  eERatioNO[ii] = ERatioBC[ii]*sqrt(pow(eEDiffFluxN[ii]/EDiffFluxN[ii],2.) + pow(eEDiffFluxO[ii]/EDiffFluxO[ii],2.));
}

for (int ii=0; ii<13; ii++) {
  Double_t tmp1     = EDiffFluxSc[ii] + EDiffFluxV[ii] + EDiffFluxTi[ii];
  ERatiosubFeFe[ii] = tmp1/EDiffFluxFe[ii];
  Double_t tmp2     = eEDiffFluxSc[ii] + eEDiffFluxV[ii] + eEDiffFluxTi[ii];
  eERatiosubFeFe[ii]= ERatiosubFeFe[ii]*sqrt(pow(tmp2/tmp1,2.) + pow(eEDiffFluxFe[ii]/EDiffFluxFe[ii],2.));
}


TGraphErrors*  HEAO_BeC = new TGraphErrors(13,EEnergy,ERatioBeC,eEEnergy,eERatioBeC);
TGraphErrors*  HEAO_BC  = new TGraphErrors(14,EEnergy,ERatioBC, eEEnergy,eERatioBC);
TGraphErrors*  HEAO_CO  = new TGraphErrors(14,EEnergy,ERatioCO, eEEnergy,eERatioCO);
TGraphErrors*  HEAO_NO  = new TGraphErrors(14,EEnergy,ERatioNO, eEEnergy,eERatioNO);
TGraphErrors*  HEAO_subFeFe = new TGraphErrors(14,EEnergy2,ERatiosubFeFe,eEEnergy2,eERatiosubFeFe);


cout << "Loading HEAO Nuclei data ...\n";
cout << "TGraphErrors are named HEAO_xx where xx can be:\n";
cout << "Be, B, C, N, C, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni\n";
cout << "BeC, BC, CO, NO, subFeFe\n";

}
