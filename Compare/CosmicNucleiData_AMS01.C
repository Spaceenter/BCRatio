{

/********************************************************************/

// AMS01 Unpublished (2009)

/********************************************************************/

Double_t  NTEnergyBC[10] = { 0.446, 0.725, 1.179, 1.915, 3.113, 5.059, 8.223, 13.364, 21.719, 35.298};
Double_t eNTEnergyBC[10] = {    0.,    0.,    0.,    0.,    0.,    0.,    0.,     0.,     0.,     0.};
Double_t  NTRatioBC[10]  = {0.3237,0.3375,0.3407,0.3039,0.2674,0.2376,0.2044, 0.1801, 0.1390, 0.1300};
Double_t eNTRatioBC[10]  = {0.0188,0.0133,0.0131,0.0133,0.0152,0.0174,0.0200, 0.0245, 0.0329, 0.0509};
TGraphErrors* AMS01_BC = new TGraphErrors(10,NTEnergyBC,NTRatioBC,eNTEnergyBC,eNTRatioBC);

cout << "Loading AMS01 Nuclei data ...\n";
cout << "TGraphErrors are named AMS01_xx where xx can be:\n";
cout << "BC\n";

}
