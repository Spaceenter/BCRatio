#ifndef MATERIAL_H
#define MATERIAL_H


#include "Element.h"


#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>


using namespace std;


class Material {

 public:

  //! Material name
  char   Name[256];
  //! Density  
  double Density; // g/cm3
  //! Element list
  vector<Element*> elem_list;
  //! Chemical formula (or percentage of number of elements)
  vector<double>   weight_list; 

 public:

  //! Ctor 
  Material(const char* name, double rho) : Density(rho) { strcpy(Name,name); }
  //! Dtor 
  ~Material() {} 

  //! Get number of elemenst
  int GetNElement() { return (int) elem_list.size(); }
  //! Get material pointer
  Element* GetElement(int ielem) { return ( (ielem>=0)&&(ielem<GetNElement()) ) ? elem_list.at(ielem) : 0; }
  //! Adding an element to the material. 
  void   AddElement(Element* elem, double weight) { elem_list.push_back(elem); weight_list.push_back(weight); }
  //! Get some info about the material
  void   Info(int more=0);

  //! Get compound atomic number (meaningful in case of chemical formula) 
  double GetA();
  //! Get compound mass (meaningful in case of chemical formula)
  double GetMass();
  //! Get averaged mass 
  double GetAverageMass();
  //! Get the material density 
  double GetDensity() { return Density; }

  //! Get material radiation lenght (weighted over elements masses)  
  double GetRadLen();
  //! Get material interaction lenght (weighted over elements masses)  
  double GetIntLen();
  //! Get material collision lenght (weighted over elements masses)  
  double GetColLen();
  //! The material lenght in units of radiation lenght for a given thikness x  
  double GetRadLenFraction(double x /*cm*/);
  //! The material lenght in units of interaction lenght for a given thikness x  
  double GetIntLenFraction(double x /*cm*/);
  //! The material lenght in units of collision lenght for a given thikness x  
  double GetColLenFraction(double x /*cm*/);

  //! Get the average cross section. 
  double GetCrossSection(int type, int Ap=0, int Zp=0, double E=0);
  //! Get the mean free path
  double GetMeanFreePath(int type, int Ap=0, int Zp=0, double E=0);
  //! Get survival probability: requires the thickness x of material, in order to find the total grammage of it
  double GetSurvivalProb(double x /*cm*/, int type, int Ap=0, int Zp=0, double E=0);

  //! Get stopping power for a given energy point [MeV cm2/g]
  double GetdEdx(double energy /* MeV */);
  //! Get the range (g/cm^2) for a given kinetic energy energy 
  // I'm using the integral definition, I don't know how to do in other ways ... [g/cm2]. 
  double GetRange(double energy /* MeV */);
  //! Get kinetic energy for a given range (range-energy relation inversion) [MeV]. 
  double GetKineticEnergyForRange(double range /* g/cm2 */);
  //! Get MIP (defined as the minimum on the dE/dx table).
  double GetMIP();

};

#endif
