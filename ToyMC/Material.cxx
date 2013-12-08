#include "Material.h"
#include "Element.h"
#include "TMath.h"


void  Material::Info(int more) {
  cout << "--------------------" << endl;
  cout << "Material:  " << Name << endl;
  cout << "Density:   " << GetDensity() << " g/cm3 " << endl;
  cout << "RadLen:    " << GetRadLen() << " g/cm2 " << endl;
  cout << "IntLen:    " << GetIntLen() << " g/cm2 " << endl;
  cout << "ColLen:    " << GetColLen() << " g/cm2 " << endl;
  cout << "NElements: " << (int)elem_list.size() << endl;
  for (int i=0; i<(int)elem_list.size(); i++) cout << elem_list.at(i)->Name << " " << weight_list.at(i) << endl;
  if (more>0) {
    for (int ii=0; ii<8; ii++) {
      cout << Element::GetCrossSectionTypeName(ii) << " Cross Section [barn] (Cex. 12 @ 1GeV/n) = " << GetCrossSection(ii,12,6,1000.) << endl;
    }
  }
  cout << "MIP value: " << GetMIP() << " MeV cm2/g" << endl; 
  cout << "--------------------" << endl;
}


double Material::GetA() {
  double m=0;
  for (int i=0; i<(int)elem_list.size(); i++) 
    m+=weight_list.at(i)*elem_list.at(i)->GetA();
  return m;
}


double Material::GetMass() {
  double m=0;
  for (int i=0; i<(int)elem_list.size();i++) 
    m+=weight_list.at(i)*elem_list.at(i)->GetMass();
  return m;
}


double Material::GetAverageMass() {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size();i++) {
    m+=weight_list.at(i)*elem_list.at(i)->GetMass();
    s+=weight_list.at(i);
  }
  return m/s;
}


double Material::GetRadLen() {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)/elem_list.at(i)->GetRadLen();
    s+=weight_list.at(i);
  }
  return s/m;
}


// This is the one to be used for computing survival probabilities, because it doesn't take into account eleastic scattering.
double Material::GetIntLen() {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)/elem_list.at(i)->GetIntLen();
    s+=weight_list.at(i);
  }
  return s/m;
}


// Collision lenght does include eleastic scattering.
double Material::GetColLen() {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)/elem_list.at(i)->GetColLen();
    s+=weight_list.at(i);
  }
  return s/m;
}


// Returns the grammage of the material of thickness x in units of rad length
double Material::GetRadLenFraction(double x) {
  return x*GetDensity()/GetRadLen();
}


// Returns the grammage of the material of thickness x in units of int length. 
// The survival probability is then easily computed as the exp of this quantity with minus sign
double Material::GetIntLenFraction(double x) {
  return x*GetDensity()/GetIntLen();
}


// Returns the grammage of the material of thickness x in units of coll length
double Material::GetColLenFraction(double x) { 
  return x*GetDensity()/GetColLen();
}


double Material::GetMeanFreePath(int type, int A_Projectile, int Z_Projectile, double E) {
  /*
    Mass:            A = A_{H_{2}O} = 2A_{H} + A_{O} = \sum_{i} w_{i} A_{i}
    Numeric density: n = \frac{N_{A} \rho}{A}
                     n_{i} = n_{H} = 2 n_{H_{2}O} = w_{i} n
                     n_{i} = \frac{w_{i} N_{A} \rho}{\sum_{i} w_{i} A_{i}}
    Mean free path:  \lambda [g/cm2] = \frac{\rho}{\sum_{i} n_{i} \sigma_{i}} 
                     \lambda = \frac{\sum_{i} w_{i} A_{i} \rho}{\sum_{i} w_{i} N_{A} \sigma_{i}} 
                     \lambda = \frac{1}{N_{A}} \times \frac{\sum_{i} w_{i} A_{i}}{\sum_{i} w_{i} \sigma_{i}}
  */ 
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)*elem_list.at(i)->GetA();
    s+=weight_list.at(i)*elem_list.at(i)->GetCrossSection(type,A_Projectile,Z_Projectile,E);
  }
  return m/(s*0.602214129);
}


// Get the cross section of the effective elemental species simulating the material. 
double Material::GetCrossSection(int type, int A_Projectile, int Z_Projectile, double E) {
  double m=0;
  for (int i=0; i<(int)elem_list.size(); i++) 
    m+=elem_list.at(i)->GetCrossSection(type,A_Projectile,Z_Projectile,E)*weight_list.at(i);
  return m;
}


double Material::GetSurvivalProb(double x, int type, int A_Projectile, int Z_Projectile, double E) {
  return exp(-x*GetDensity()/GetMeanFreePath(type,A_Projectile,Z_Projectile,E));
}


double Material::GetRange(double energy) {
  // per avere la tabella di energia devo utilizzare un elemento come riferimento ...
  Element* ref = elem_list.at(0); // used as reference
  int lastbin = min(ref->GetLowerIndex(energy),DEDXTABLE);
  if (lastbin==0) return 0;
  // first bin
  int    index = 0;
  double de = ref->GetEnergy(index);
  double m1 = 0; 
  double m2 = 0;
  double s  = 0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m1 += weight_list.at(i)*elem_list.at(i)->GetdEdx(index);
    s  += weight_list.at(i);
  }
  double invdedx = s/m1;
  double sum     = invdedx*de;
  // other bins
  for (index=1; index<=lastbin; index++) {
    de = ref->GetEnergy(index) - ref->GetEnergy(index-1);
    m1 = 0.;
    m2 = 0.;
    s =  0.;
    for (int i=0; i<(int)elem_list.size(); i++) {
      m1 += weight_list.at(i)*elem_list.at(i)->GetdEdx(index);
      m2 += weight_list.at(i)*elem_list.at(i)->GetdEdx(index-1); 
      s  += weight_list.at(i);
    }
    invdedx = s*(1/m1 + 1/m2)/2;
    sum += invdedx*de;
  }
  // last bin
  index = lastbin;
  de = energy - ref->GetEnergy(index);
  m1 = 0.;
  m2 = 0.;
  s  = 0.;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m1 += weight_list.at(i)*elem_list.at(i)->GetdEdx(index);
    m2 += weight_list.at(i)*elem_list.at(i)->GetdEdx(index-1);
    s  += weight_list.at(i);
  }
  invdedx = s*(1/m1 + 1/m2)/2;
  sum += invdedx*de;
  return sum;
}


double Material::GetdEdx(double energy) {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)*elem_list.at(i)->GetdEdx(energy);
    s+=weight_list.at(i);
  }
  return m/s;
}


double Material::GetKineticEnergyForRange(double range) {
  Element* ref = elem_list.at(0); // used as reference
  for (int i=0; i<DEDXTABLE; i++) {
    double energy1 = ref->GetEnergy(i);
    double energy2 = ref->GetEnergy(i+1);
    double range1  = GetRange(energy1);
    double range2  = GetRange(energy2);
    if ( (range<GetRange(energy1))||(range>GetRange(energy2)) ) continue;
    double width = range2 - range1;
    double wu    = (range2 - range)/width;
    double wl    = (range - range1)/width;
    double eu    = energy2;
    double el    = energy1;
    return wl*eu + wu*el;
  }
  return 0.;
}


double Material::GetMIP() {
  double m=0;
  double s=0;
  for (int i=0; i<(int)elem_list.size(); i++) {
    m+=weight_list.at(i)*elem_list.at(i)->GetMIP();
    s+=weight_list.at(i);
  }
  return m/s;
}
