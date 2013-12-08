#ifndef ELEMENT_H
#define ELEMENT_H


#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>


#define NCS 10
#define DEDXTABLE 132
#define NA 293
#define NZ 118


using namespace std;


class Element {

 public:

  ////////////////////////////////////
  // Members
  ////////////////////////////////////

  // Name of the element
  char   Name[256];
  // Mass number
  int    A;
  // Atomic number
  int    Z;
  // Radiation lenght [g/cm2]
  double RadLen; 
  // Interaction lenght [g/cm2]
  double IntLen; 
  // Collision lenght [g/cm2] (ColCrossSec > IntCrossSec)
  double ColLen;

  ////////////////////////////////////
  // Stopping power and range tables 
  ////////////////////////////////////

  // Stopping power and Range tables (for protons), Kinetic Energy binning [MeV]
  double Energy[DEDXTABLE];
  // Stopping power and Range tables (for protons), dEdx [MeV cm2/g]
  double dEdx[DEDXTABLE];
  // Stopping power and Range tables (for protons), range [g/cm2]
  double Range[DEDXTABLE];  
 
  ////////////////////////////////////
  // Cross section part
  ////////////////////////////////////

  enum CSType { 
    kInteraction   = 0,
    kCollision     = 1, 
    kGeometric     = 2,
    kBradtPeters   = 3,
    kSihver        = 4,
    kKox           = 5,
    kShen          = 6,
    kTripathi      = 7,
    kTripathiUniv  = 8,
    kGlauber       = 9
  };

  // Pi 
  static double PI;
  // Atomic Mass Unit, Mass(C12)/12 [MeV]
  static double AMU;

  ////////////////////////////////////
  // Mass and Radius tables
  ////////////////////////////////////

  //! Mass table [AMU]
  static double MassTable[NA][NZ];
  //! Radius table [fm]
  static double RadiusTable[NA][NZ];

 public:

  //! c-tor
  Element(const char* name, int a, int z, double rl, double il, double cl) : A(a), Z(z), RadLen(rl), IntLen(il), ColLen(cl) { strcpy(Name,name); }
  //! d-tor
  ~Element() {}
  //! Get mass number (
  int    GetA() { return A; }
  //! Get atomic number 
  int    GetZ() { return Z; }
  //! Get mass [GeV/c2] 
  double GetMass() { return GetMass(GetZ(),GetA()); }                  
  //! Get radius [fm]
  double GetRadius() { return GetRadius(GetZ(),GetA()); } 
  //! Print info
  void   Info(int more=0);

  //! Get radiation lenght (g/cm2)
  double GetRadLen() { return RadLen; } 
  //! Get the interaction lenght (g/cm2)
  double GetIntLen() { return IntLen; } 
  //! Get collision cross section (g/cm2)
  double GetColLen() { return ColLen; } // g/cm2

  ///////////////////////////////////////
  // Mass and Radius evaluation 
  ///////////////////////////////////////

  //! Read mass table file 
  static void   ReadMassTableFile(char* filename);
  //! Get approximate mass [GeV/c2]
  static double GetApproxMass(int A) { return A*AMU; }
  //! Get mass [GeV/c2] 
  static double GetTableMass(int Z, int A);
  //That'the one to be used
  //! Get mass best evaluator [GeV/c2] 
  static double GetMass(int Z, int A); 

  //! Radius evaluation [fm]
  static double GetClassicRadius(int A) { return 1.2*pow(A,1./3.); }
  //! Liquid drop model (Meyer) radius [fm]
  static double GetMeyerRadius(int A) { return sqrt(3./5.)*pow(A,1./3.)*(1.15 + 1.8*pow(A,-2./3.) - 1.2*pow(A,-4./3.)); }
  //! Read radius table file
  static void   ReadRadiusTableFile(char* filename);
  //! Get radius best evaluator [fm]
  static double GetTableRadius(int Z, int A);
  //! Get radius best evaluator [fm]
  //Ok, the one to be used
  static double GetRadius(int Z, int A);

  ///////////////////////////////////////
  // Cross Section
  ///////////////////////////////////////

  //! Return cross section name by code 
  static const char* GetCrossSectionTypeName(int type);

  //! Generic method for the total reaction cross section (s_R = s_T - s_e) retrieval (E is kinetic energy [MeV/n] of projectile) [barn]
  double GetCrossSection(int type, int Ap, int Zp, double E);
  //! Generic method for mean free path retrieval (E is kinetic energy [MeV/n] of projectile) [cm]
  double GetMeanFreePath(int type, int Ap, int Zp, double E);
  //! Generic method for survival probability retrieval (x is width (cm); rho is density (g/cm3); E is kinetic energy [MeV/n])
  double GetSurvivalProb(double x, double rho, int type, int Ap, int Zp, double E);

  //! Get total cross section (s_i = s_T - s_e - s_qe), valid only for protons/neutrons [barn] 
  double GetIntCrossSection() { return 10*GetMass()/(GetIntLen()*6.022); }
  //! Get collision cross section (s_c = s_T), valid only for protons/neutrons [barn] 
  double GetColCrossSection() { return 10*GetMass()/(GetColLen()*6.022); }
  //! Get Sihver total reaction cross section, independent from energy [barn]
  double GetSihverCrossSection(int Ap);        
  //! Get geometric total reaction cross section, independent from energy [barn]
  double GetGeometricCrossSection(int Ap); 
  //! Get Bradt-Peters total reaction cross section, independent from energy [barn]
  double GetBradtPetersCrossSection(int Ap);    
  //! Get Kox total reaction cross section [barn] 
  double GetKoxCrossSection(int Ap, int Zp, double E); 
  //! Get Shen total reaction cross section [barn] 
  double GetShenCrossSection(int Ap, int Zp, double E); 
  //! Get Tripathi total reaction cross section, accurate for heavy elements [barn] 
  double GetTripathiCrossSection(int Ap, int Zp, double E);
  //! Get Tripathi total reaction cross section, universal formula [barn]
  double GetTripathiUnivCrossSection(int Ap, int Zp, double E); // TO-BE-FIX: problem with the Rc parameters, the Coulomb multipliers
  //! Get Glauber Gribov cross section
  double GetGlauberCrossSection(int Ap, int Zp, double E);
  
  //! Get kinetic energy in center-of-mass system [MeV]
  double GetKEcm(int Ap, int Zp, double E);
  //! Get the Kox-Shen c-parameter
  double GetC(double E);
 
  ///////////////////////////////////////
  // Energy loss
  ///////////////////////////////////////

  //! Read the stopping power and range tables (from PSTAR)
  void   ReadStoppingPowerFile(char* filename);
  //! Get the table value
  double GetEnergy(int i) { return Energy[i]; }
  //! Get the table value
  double GetdEdx(int i)   { return dEdx[i];   }
  //! Get the table value
  double GetRange(int i)  { return Range[i];  }
  //! Lower edge  
  int    GetLowerIndex(double energy);
  //! Upper edge 
  int    GetUpperIndex(double energy);
  //! Get stopping power for a given energy point (from log10-interpolation) 
  double GetdEdx(double energy ); //MeV
  //! Get the range (g/cm^2) for a given kinetic energy energy
  double GetRange(double energy, int type = 0) { return (type==0) ? GetRangeFromInterpolation(energy) : GetRangeFromdEdxIntegral(energy); }
  //! Get the range (g/cm^2) for a given kinetic energy (range table log10-interpolation)  
  double GetRangeFromInterpolation(double energy);
  //! Get the range (g/cm^2) for a given kinetic energy (dE/dx integral)  
  double GetRangeFromdEdxIntegral(double energy);
  //! Get kinetic energy for a given range (range-energy relation inversion)
  double GetKineticEnergyForRange(double range, int type = 0);
  //! Get MIP (defined as the minimum on the dE/dx table)
  double GetMIP();

};

#endif
