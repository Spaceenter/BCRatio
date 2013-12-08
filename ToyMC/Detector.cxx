#include "Detector.h"
#include "Element.h"
#include "TMath.h"
#include <vector>


void  Detector::Info(int more) {
  cout << "--------------------" << endl;
  cout << "Detector: " << Name << endl;
  cout << "NMaterials: " << (int) mat_list.size() << endl;
  cout << "Total grammage:  " << GetGrammage(0) << " g/cm2 (theta=0)" << endl;
  cout << "MIP Energy Loss: " << GetMipEnergyLoss(0) << " MeV (theta=0)" << endl;
  cout << "Min. Energy:     " << GetMinEnergy(0) << " MeV (theta=0)" << endl;
  cout << "Survival Probability (interaction lenght --> CStot - CSqe) for a vertical proton: " << GetSurvivalProb(0,0.,1,1,10000) << endl;
  cout << "Survival Probability (collision lenght --> CStot) for a vertical proton:          " << GetSurvivalProb(1,0.,1,1,10000) << endl;
  cout << "x/X0:   " << GetRadLenFraction(0) << " (theta=0)" << endl;
  cout << "x/Xint: " << GetIntLenFraction(0) << " (theta=0)" << endl;
  cout << "x/Xcol: " << GetColLenFraction(0) << " (theta=0)" << endl;
  if (more>1) {
    printf("n.                           Name Thickness Grammage    x/X0  x/Xint  x/Xcol\n");
    printf("                                       (cm) (g/cm^2)     (%%)     (%%)     (%%)\n");
    for (int i=0; i<(int)mat_list.size(); i++)
      printf("%2d %30s   %7.3f  %7.3f %7.3f %7.3f %7.3f\n",
             i,mat_list.at(i)->Name,thickness_list.at(i),thickness_list.at(i)*mat_list.at(i)->Density,
             100.*thickness_list.at(i)*mat_list.at(i)->Density/mat_list.at(i)->GetRadLen(),
             100.*thickness_list.at(i)*mat_list.at(i)->Density/mat_list.at(i)->GetIntLen(),
	     100.*thickness_list.at(i)*mat_list.at(i)->Density/mat_list.at(i)->GetColLen());
  }
  if (more>0) {
    for (int ii=0; ii<NCS; ii++) {
      cout << Element::GetCrossSectionTypeName(ii) << " survival probability (ex. C12 @ 1GeV/n) = " << GetSurvivalProb(ii,0.,12,6,1000.) << endl;
    }
  }
  if (more>2) {
    for (int i=0; i<(int)mat_list.size(); i++) { 
      mat_list.at(i)->Info();
    }
  }
  cout << "--------------------" << endl;
}


double  Detector::GetSurvivalProb(int type, double theta /* rad */, int A_Projectile, int Z_Projectile, double E /* kinetic energy [MeV/n] */) {
  double prob = 1;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta); // thickness [cm]
    prob *= mat_list.at(i)->GetSurvivalProb(x,type,A_Projectile,Z_Projectile,E);
  }
  return prob;
}


double Detector::GetMeasuredBoverC(int type, double theta, double B_0, double C_0, double E) {   
  double crossCB=0.07;
  double lambdaC, lambdaB, lambdaCB;
  double B_x, C_x;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta);
    double grammage= x*mat_list.at(i)->GetDensity();
    lambdaC  = mat_list.at(i)->GetMeanFreePath(type,12,6,E);
    lambdaB  = mat_list.at(i)->GetMeanFreePath(type,11,5,E);
    lambdaCB = 10*mat_list.at(i)->GetAverageMass()/(6.022*crossCB);
    B_x = exp(-grammage/lambdaB)*B_0+exp(-grammage/lambdaB)*C_0*lambdaC*lambdaB/(lambdaCB*(lambdaC-lambdaB))*(exp(grammage*(lambdaC-lambdaB)/(lambdaC*lambdaB))-1);
    C_x = exp(-grammage/lambdaC)*C_0;
    C_0 = C_x;
    B_0 = B_x; 
  }
  return B_x/C_x;
}


double Detector::GetGrammage(double theta) {
  double grammage = 0.;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta); // thickness [cm]
    grammage += x*mat_list.at(i)->GetDensity(); // grammage [g/cm2]
  }
  return grammage;
}


double Detector::GetRadLenFraction(double theta) {
  double sum = 0;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta); // thickness [cm]
    sum += mat_list.at(i)->GetRadLenFraction(x);
  }
  return sum;
}


double Detector::GetIntLenFraction(double theta) {
  double sum = 0;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta); // thickness [cm]
    sum += mat_list.at(i)->GetIntLenFraction(x);
  }
  return sum;
}


double Detector::GetColLenFraction(double theta) {
  double sum = 0;
  for (int i=0; i<(int)mat_list.size(); i++) {
    double x = thickness_list.at(i)/cos(theta); // thickness [cm]
    sum += mat_list.at(i)->GetColLenFraction(x);
  }
  return sum;
}


//Qui dentro andrebbe implementata tutta la robaccia della conversione C in B e ridecadimento di B. La cosa andrebbe iterata strato per strato! 
//Che energia dare al boro così prodotto? Si potrebbe approssimare lasciando la solita energia del carbonio in entrata, almeno come prima approssimazione dovrebbe andare, visto la debole dipendenza in energia delle sezioni d'urto per energie piuttosto elevante (>GeV/n sicuramente OK)
double Detector::GetBoverC(double B_0, double C_0, int type, double theta, double E) {   
  return B_0*GetSurvivalProb(type,theta,11,5,E)/(C_0*GetSurvivalProb(type,theta,12,6,E));
}


double Detector::GetMinEnergy(double theta /* rad */) {
  return GetToiEnergy(0,theta);
}


double Detector::GetMipEnergyLoss(double theta /* rad */) {
  double mip = 0.;
  for (int i=0; i<(int)mat_list.size(); i++) {
    mip += thickness_list.at(i)*mat_list.at(i)->Density/cos(theta)*mat_list.at(i)->GetMIP();
  }
  return mip;
}


double Detector::GetToiEnergy(double k, double theta /* rad */) {
  double xtop = 0.;
  double xbot = 0.;
  double etop = k;
  double ebot = k;
  for (int i=(int)mat_list.size()-1; i>=0; i--) {
    xbot = xtop;
    xtop = xbot + thickness_list.at(i)*mat_list.at(i)->Density/cos(theta);
    ebot = etop;
    etop = mat_list.at(i)->GetKineticEnergyForRange(xtop-xbot+mat_list.at(i)->GetRange(ebot)); 
    // etop = ebot + mat_list.at(i)->GetKineticEnergyForRange(xtop) - mat_list.at(i)->GetKineticEnergyForRange(xbot);
    // if (VERBOSE>0) 
    // printf("%2d  %20s    xbot=%8.4f g/cm2    xtop=%8.4f g/cm2    dx=%8.4f g/cm2    Ebot=%8.4f MeV    Etop=%8.4f MeV   dE=%8.4f MeV \n",i,mat_list.at(i)->Name,xbot,xtop,xtop-xbot,ebot,etop,etop-ebot);
  }
  return etop;
}




