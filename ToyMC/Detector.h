#ifndef DETECTOR_H
#define DETECTOR_H


#include "Material.h"
#include "Element.h"


#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>


#define VERBOSE 0


using namespace std;

class Detector {

 public:

  //! Name of the detector
  char Name[256];
  //! Material list: the detector is given by a set of different materials
  vector<Material*> mat_list;
  //! Thickness list (cm): each material of the detector, which has its own density, has also his thickness
  vector<double> thickness_list; 

 public:

  //! Ctor
  Detector(const char* name) { strcpy(Name,name); }
  //! Dtor
  ~Detector() {} 

  //! Get number of materials
  int GetNMaterial() { return (int) mat_list.size(); }
  //! Get material pointer
  Material* GetMaterial(int imat) { return ( (imat>=0)&&(imat<GetNMaterial()) ) ? mat_list.at(imat) : 0; }

  //! Add a material to the detector
  void   AddMaterial(Material* mat, double thickness) { mat_list.push_back(mat); thickness_list.push_back(thickness); }
  //! Some detector information
  void   Info(int more=0);

  //ATTENZIONE PER ALCUNI VALORI DI THETA SPUTA FUORI PROB MAGGIORI DI 1!!!!!!!!!Da sistemare mettendo limitazioni su theta!!!
  //! Get survival probability. Here we introduce the collision angle.
  double GetSurvivalProb(int type, double theta = 0/* rad */, int A_Projectile=0,int Z_Projectile=0, double E=0 /* kinetic energy [MeV/n] */);
  //! Get total grammage for a given inclination
  double GetGrammage(double theta = 0/* rad */);
  
  //! Detector thickness in radiation lenght units (for a given inclination)
  double GetRadLenFraction(double theta = 0);
  //! Detector thickness in interaction lenght units (for a given inclination)
  double GetIntLenFraction(double theta = 0);
  //! Detector thickness in collision lenght units (for a given inclination)
  double GetColLenFraction(double theta = 0);

  //return the B over C fraction after passage through the detector, for an initial B and C composition
  double GetBoverC(double B_0, double C_0, int type, double theta, double E);
  double GetMeasuredBoverC(int type, double theta, double B_0, double C_0, double E);

  //! Get the kinetic energy at the top of the instrument
  double GetToiEnergy(double k, double theta = 0/* rad */);
  //! Get the proton minimum energy to traverse downgoind the entire detector (with angle theta)
  double GetMinEnergy(double theta = 0/* rad */);
  //! Get the MIP energy loss (only a simplified evaluation)
  double GetMipEnergyLoss(double theta = 0/* rad */); 
  //! Get the energy loss from a downgoing particle (from material 0 to n) with energy and inclination
  double GetEnergyLoss(double energy, double theta = 0) { return 0; }
  
};

#endif

