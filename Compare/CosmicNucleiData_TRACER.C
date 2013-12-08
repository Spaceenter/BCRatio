{

/********************************************************************/

// M. Ave, ..., D. Muller, ... et al. (2008)
// The Astrophysical Journal, Volume 678, Issue 1, pp. 262-273. 
// TRACER transition radiation detector on balloon  

/********************************************************************/

// Energy GeV/n
/* Double_t  TEnergyBins[6] = {...}; ... see article table */

// Energy GeV/n, Flux in (m^2 sr s GeV/n)^-1
Double_t   TEnergyC   [8] = {0.9, 1.1, 1.5, 48, 1834, 4001, 7872, 15611};
Double_t  eTEnergyC   [8] = { 0.,  0.,  0., 0.,   0.,   0.,   0.,    0.};
Double_t   TDiffFluxC [8] = {1.71, 1.38, 8.95e-1, 7.27e-4, 4.4e-8, 5.8e-9, 9.1e-10, 2.0e-10};
Double_t ehTDiffFluxC [8] = {0.01, 0.01, 0.05e-1, 0.01e-4, 0.5e-8, 1.1e-9, 6.0e-10, 1.2e-10};
Double_t elTDiffFluxC [8] = {0.01, 0.01, 0.05e-1, 0.01e-4, 0.5e-8, 1.1e-9, 3.9e-10, 0.8e-10};
TGraphAsymmErrors* TRACER_C = new TGraphAsymmErrors(8,TEnergyC,TDiffFluxC,eTEnergyC,eTEnergyC,elTDiffFluxC,ehTDiffFluxC);

Double_t   TEnergyNe  [9] = {0.9, 1.1, 1.5, 48, 1080, 1834, 4001, 7872, 15611};
Double_t  eTEnergyNe  [9] = { 0.,  0.,  0., 0.,   0.,   0.,   0.,   0.,    0.};
Double_t   TDiffFluxNe[9] = {3.08e-1, 2.42e-1, 1.48e-1, 1.21e-4, 3.2e-8, 8.2e-9, 8.9e-10, 2.2e-10, 7.8e-11};
Double_t ehTDiffFluxNe[9] = {0.05e-1, 0.04e-1, 0.02e-1, 0.01e-4, 0.8e-8, 2.3e-9, 6.9e-10, 5.0e-10,10.1e-11};
Double_t elTDiffFluxNe[9] = {0.05e-1, 0.04e-1, 0.02e-1, 0.01e-4, 0.8e-8, 2.3e-9, 4.2e-10, 1.8e-10, 5.0e-11};
TGraphAsymmErrors* TRACER_Ne = new TGraphAsymmErrors(9,TEnergyNe,TDiffFluxNe,eTEnergyNe,eTEnergyNe,elTDiffFluxNe,ehTDiffFluxNe);

Double_t   TEnergyMg  [8] = {0.9, 1.1, 1.5, 48, 946, 1834, 4001, 10113};
Double_t  eTEnergyMg  [8] = { 0.,  0.,  0., 0.,  0.,   0.,   0.,    0.};
Double_t   TDiffFluxMg[8] = {4.02e-1, 2.96e-1, 1.89e-1, 1.60e-4, 3.6e-8, 6.3e-9, 2.4e-10, 6.7e-11};
Double_t ehTDiffFluxMg[8] = {0.06e-1, 0.05e-1, 0.02e-1, 0.01e-4, 0.7e-8, 2.9e-9, 5.6e-10,15.4e-11};
Double_t elTDiffFluxMg[8] = {0.06e-1, 0.05e-1, 0.02e-1, 0.01e-4, 0.7e-8, 2.1e-9, 2.0e-10, 5.5e-11};
TGraphAsymmErrors* TRACER_Mg = new TGraphAsymmErrors(8,TEnergyMg,TDiffFluxMg,eTEnergyMg,eTEnergyMg,elTDiffFluxMg,ehTDiffFluxMg);

Double_t   TEnergySi  [8] = {0.9, 1.1, 1.5, 28, 138, 627, 946, 2088};
Double_t  eTEnergySi  [8] = { 0.,  0.,  0., 0.,  0.,  0.,  0.,   0.};
Double_t   TDiffFluxSi[8] = {3.07e-1, 2.27e-1, 1.39e-1, 5.23e-4, 7.2e-6, 1.3e-7, 4.0e-8, 2.2e-9};
Double_t ehTDiffFluxSi[8] = {0.06e-1, 0.04e-1, 0.02e-1, 0.03e-4, 0.2e-6, 0.3e-7, 0.8e-8, 1.3e-9};
Double_t elTDiffFluxSi[8] = {0.06e-1, 0.04e-1, 0.02e-1, 0.03e-4, 0.2e-6, 0.3e-7, 0.8e-8, 0.9e-9};
TGraphAsymmErrors* TRACER_Si = new TGraphAsymmErrors(8,TEnergySi,TDiffFluxSi,eTEnergySi,eTEnergySi,elTDiffFluxSi,ehTDiffFluxSi);

Double_t   TEnergyS  [9] = {0.9, 1.1, 1.5, 2.0, 21,70, 223, 644, 1257};
Double_t  eTEnergyS  [9] = { 0.,  0.,  0.,  0., 0.,0.,  0.,  0.,   0.};
Double_t   TDiffFluxS[9] = {7.5e-2, 6.3e-2, 3.3e-2, 2.0e-2, 2.1e-4, 1.1e-5, 5.2e-7, 2.1e-8, 2.6e-9};
Double_t ehTDiffFluxS[9] = {0.6e-2, 0.5e-2, 0.3e-2, 0.2e-2, 0.1e-4, 0.1e-5, 0.9e-7, 1.6e-8, 2.5e-9};
Double_t elTDiffFluxS[9] = {0.6e-2, 0.5e-2, 0.3e-2, 0.2e-2, 0.1e-4, 0.1e-5, 0.9e-7, 1.0e-8, 1.4e-9};
TGraphAsymmErrors* TRACER_S = new TGraphAsymmErrors(9,TEnergyS,TDiffFluxS,eTEnergyS,eTEnergyS,elTDiffFluxS,ehTDiffFluxS);

Double_t   TEnergyAr  [9] = {0.9, 1.1, 1.5, 2.0, 21,70, 223, 558, 1105};
Double_t  eTEnergyAr  [9] = { 0.,  0.,  0.,  0., 0.,0.,  0.,  0.,   0.};
Double_t   TDiffFluxAr[9] = {4.0e-2, 2.9e-2, 1.5e-2, 8.9e-3, 7.1e-5, 3.8e-6, 1.5e-7, 1.9e-8, 3.2e-9};
Double_t ehTDiffFluxAr[9] = {0.5e-2, 0.4e-2, 0.2e-2, 1.4e-3, 0.3e-5, 0.4e-6, 0.5e-7, 2.5e-8, 7.3e-9};
Double_t elTDiffFluxAr[9] = {0.5e-2, 0.4e-2, 0.2e-2, 1.4e-3, 0.3e-5, 0.4e-6, 0.5e-7, 1.3e-8, 2.6e-9};
TGraphAsymmErrors* TRACER_Ar = new TGraphAsymmErrors(9,TEnergyAr,TDiffFluxAr,eTEnergyAr,eTEnergyAr,elTDiffFluxAr,ehTDiffFluxAr);

Double_t   TEnergyCa  [9] = {0.9, 1.1, 1.5, 2.0, 21,70, 223, 558, 1105};
Double_t  eTEnergyCa  [9] = { 0.,  0.,  0.,  0., 0.,0.,  0.,  0.,   0.};
Double_t   TDiffFluxCa[9] = {5.5e-2, 3.7e-2, 2.2e-2, 1.3e-2, 9.9e-5, 4.1e-6, 2.8e-7, 1.8e-8, 1.9e-9};
Double_t ehTDiffFluxCa[9] = {0.6e-2, 0.4e-2, 0.2e-2, 0.2e-2, 0.4e-5, 0.4e-6, 0.6e-7, 1.7e-8, 4.3e-9};
Double_t elTDiffFluxCa[9] = {0.6e-2, 0.4e-2, 0.2e-2, 0.2e-2, 0.4e-5, 0.4e-6, 0.6e-7, 1.0e-8, 1.6e-9};
TGraphAsymmErrors* TRACER_Ca = new TGraphAsymmErrors(9,TEnergyCa,TDiffFluxCa,eTEnergyCa,eTEnergyCa,elTDiffFluxCa,ehTDiffFluxCa);

Double_t   TEnergyFe  [11] = {0.9, 1.1, 1.5, 2.0, 16, 30, 60, 183, 566, 733, 1421};
Double_t  eTEnergyFe  [11] = { 0.,  0.,  0.,  0., 0., 0., 0.,  0.,  0.,  0.,   0.};
Double_t   TDiffFluxFe[11] = {1.21e-1, 1.65e-1, 1.02e-1, 6.2e-2, 1.22e-3, 2.80e-4, 4.7e-5, 1.9e-6, 1.3e-7, 5.2e-8, 7.4e-9};
Double_t ehTDiffFluxFe[11] = {0.05e-1, 0.04e-1, 0.02e-1, 0.2e-2, 0.02e-3, 0.04e-4, 0.1e-5, 0.1e-6, 0.6e-9, 1.5e-8, 2.2e-9};
Double_t elTDiffFluxFe[11] = {0.05e-1, 0.04e-1, 0.02e-1, 0.2e-2, 0.02e-3, 0.04e-4, 0.1e-5, 0.1e-6, 0.4e-9, 1.5e-8, 2.2e-9};
TGraphAsymmErrors* TRACER_Fe = new TGraphAsymmErrors(11,TEnergyFe,TDiffFluxFe,eTEnergyFe,eTEnergyFe,elTDiffFluxFe,ehTDiffFluxFe);

/********************************************************************/

// A. Obermeier et al.
// The Astrophysical Journal, 742:14 (11pp), 2011 November 20
// ENERGY SPECTRA OF PRIMARY AND SECONDARY COSMIC-RAY NUCLEI MEASURED WITH TRACER

/********************************************************************/

Double_t   OEnergyBC[6] = {0.9,1.1,1.5,23.8,173,2070};
Double_t  eOEnergyBC[6] = {0.,0.,0.,0.,0.,0.};
Double_t   ORatioBC[6]  = {0.36,0.33,0.32,0.180,0.082,0.12}; 
Double_t ehORatioBC[6]  = {0.022,0.022,0.022,0.009,0.016,0.15};
Double_t elORatioBC[6]  = {0.022,0.022,0.022,0.098,0.014,0.09};
TGraphAsymmErrors* TRACER_BC = new TGraphAsymmErrors(9,OEnergyBC,ORatioBC,eOEnergyBC,eOEnergyBC,elORatioBC,ehORatioBC);

cout << "Loading TRACER Nuclei data ...\n";
cout << "TGraphErrors are named TRACER_xx where xx can be:\n";
cout << "C, Ne, Mg, Si, S, Ar, Ca, Fe\n";
cout << "BC" << endl;
}
