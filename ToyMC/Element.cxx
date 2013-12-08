#include "Element.h"


// const 
double Element::PI = 3.14159265;
double Element::AMU = 0.931494061;
double Element::MassTable[NA][NZ] = {{0}};
double Element::RadiusTable[NA][NZ] = {{0}};


void Element::Info(int more) {
  cout << "--------------------" << endl;
  cout << "Element: " << Name << endl;
  cout << "A: " << GetA() << endl;
  cout << "Z: " << GetZ() << endl;
  cout << "Mass: " << GetMass() << " GeV/c2" << endl; 
  cout << "Radius: " << GetRadius() << " fm" << endl;
  cout << "RadLen: " << GetRadLen() << " g/cm2" << endl;
  cout << "IntLen: " << GetIntLen() << " g/cm2" << endl;
  cout << "ColLen: " << GetColLen() << " g/cm2" << endl;
  if (more>0) {
    for (int ii=0; ii<10; ii++) {
      cout << GetCrossSectionTypeName(ii) << " Cross Section [barn] (ex. C12 @ 1GeV/n) = " << GetCrossSection(ii,12,6,1000.) << endl;
    }
  }
  cout << "--------------------" << endl;
}


///////////////////////////////////
// Mass and Radius Table
///////////////////////////////////


double Element::GetTableMass(int Z, int A) {
  double mass = 0;
  if ( (Z>0)&&(Z<=NZ)&&(A>0)&&(A<=NA) ) mass = MassTable[A-1][Z-1];
  return mass;
}


// Get the tabulated masses as function of Z and A, if not tabulated use the approximate formula.
double Element::GetMass(int Z, int A) {
  double mass = GetTableMass(Z,A);
  if (mass<=0) mass = GetApproxMass(A); 
  return mass;
}


// Read masses from tabulated values.
void Element::ReadMassTableFile(char* filename) {
  FILE* file = fopen(filename,"r");
  if (file==0) {
    printf("Error opening file %s\n",filename);
    return;
  }
  int i=0;
  while (!(feof(file))) {
    int Z,A;
    float M;
    fscanf(file,"%d%d%f",&Z,&A,&M);
    if ( (Z>0)&&(Z<=NZ)&&(A>0)&&(A<=NA) ) MassTable[A-1][Z-1] = M*AMU;
    i++;
  }
}


// Get the tabulated radii as function of Z and A, if not tabulated use the approximate formula.
double Element::GetTableRadius(int Z, int A) {
  double radius = 0;
  if ( (Z>0)&&(Z<=NZ)&&(A>0)&&(A<=NA) ) radius = RadiusTable[A-1][Z-1];
  return radius;
}


// Get radius from Meyer approximate formula if not tabulated, otherwise use the (most precise) tabulated ones.
double Element::GetRadius(int Z, int A) {
  double radius = GetTableRadius(Z,A);
  // if no element in the database use Meyer formula
  if (radius<=0) radius = GetMeyerRadius(A);
  return radius;
}


// Read the tabulated radii.
void Element::ReadRadiusTableFile(char* filename) {
  FILE* file = fopen(filename,"r");
  if (file==0) {
    printf("Error opening file %s\n",filename);
    return;
  }
  int i=0;
  while (!(feof(file))) {
    int Z,A;
    float R;
    fscanf(file,"%d%d%f",&Z,&A,&R);
    if ( (Z>0)&&(Z<=NZ)&&(A>0)&&(A<=NA) ) RadiusTable[A-1][Z-1] = R;
    i++;
  }
}


///////////////////////////////////
// Cross Section
///////////////////////////////////


// Define the names of the different cross section models.
const char* Element::GetCrossSectionTypeName(int type) {
  const char* names[10] = {
    "Interaction  ",
    "Collision    ",
    "Geometric    ",
    "BrandPeters  ",
    "Sihver       ",
    "Kox          ",
    "Shen         ",
    "Tripathi     ",
    "TripathiUniv ",
    "Glauber      "
  };
  return names[type];
}


// Geometrical cross section (energy independent and without Coulomb barrier). 
double Element::GetGeometricCrossSection(int Ap) { 
  double At = GetA();
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  double r0 = 1.35; // fm
  double CS = PI*r0*r0*pow(cubicAt+cubicAp,2.);
  return CS*0.01; // fm^2 -> barn
}


// Add to the geometrical cross section the Coulomb barrier, but in an energy independent way. 
// The result is a smaller cross section due to Coulomb repulsion of nuclei. 
// Recall that the Coulomb repulsion is important for low energies, tipically below a ten Mev per nucleon.
double Element::GetBradtPetersCrossSection(int Ap) {  
  double At = GetA();
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);  
  double r0 = 1.35; // fm
  double b = 0.83;
  double CS = PI*r0*r0*pow(cubicAt+cubicAp-2*b,2.);
  return CS*0.01; // fm^2 -> barn
}


// Independent on energy: good for E>100Mev/n. B_0 transparency parameter, describing the nucleus nucleus overlap
double Element::GetSihverCrossSection(int Ap) { 
  double At = GetA();
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  double r0 = 1.36; // fm
  double b0 = 1.581 - 0.876*(1./cubicAt + 1./cubicAp); // transparency
  double CS = PI*r0*r0*pow((cubicAt + cubicAp - b0*(1./cubicAt + 1./cubicAp)),2.);
  return CS*0.01; // fm^2 -> barn
}


// Get the center of mass energy, giving the projectile A and Z and its kinetic energy in MeV per nucleon
double Element::GetKEcm(int Ap, int Zp, double E) {
  // double At = GetA();
  // double Zt = GetZ();
  double Mt = GetMass();
  double Mp = GetMass(Zp,Ap);
  double KElab = E*Ap; // Ekin/n -> Ekin
  double Elab = KElab + Mp; // Ekin -> Etot
  double Plab = sqrt(Elab*Elab - Mp*Mp); // Etot -> P
  double Ecm = sqrt(Mp*Mp + Mt*Mt + 2*Elab*Mt); // Elab -> Ecm
  double Pcm = Plab*Mt/Ecm; // Ecm -> Pcm
  double KEcm = sqrt(Pcm*Pcm + Mp*Mp) - Mp; // Ekincm
  return KEcm;
}


// Return an energy-dependent parameter required for Kox and Shen cross section parametrizations
double Element::GetC(double E) {
  double x = log10(E);
  return (x<1.5) ? (-10/pow(1.5,5.) + 2.0)*pow(x/1.5,3.) : -10/pow(x,5.) + 2.0;
}


// Energy dependent Kox cross section
double Element::GetKoxCrossSection(int Ap, int Zp, double E) {
  // source/processes/hadronic/cross_sections/src/G4IonsKoxCrossSection.cc 
  double At = GetA();
  double Zt = GetZ();
  double rc = 1.3; // fm
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  // Coulomb barrier
  double KEcm = GetKEcm(Ap,Zp,E);  
  double Bc = Zt*Zp/(rc*(cubicAp + cubicAt));
  // Coulomb barrier could prevent the reaction
  if (KEcm<=Bc) return 0;
  // Radius
  double r0 = 1.1; // fm
  double Rvol = r0*(cubicAp + cubicAt); // fm
  double a = 1.85;
  double c = GetC(E);
  double D;
  D = 5.*(At - 2*Zt)*Zp/(Ap*At); // fm
  // quite an unjustified solution in order to symmetrize the cross section between projectile and target.
  // if(At>=Ap) D = 5.*(At - 2*Zt)*Zp/(Ap*At);;
  // else       D = 5.*(Ap - 2*Zp)*Za/(Ap*At);
  double Rsurf = r0*((a*cubicAp*cubicAt)/(cubicAp + cubicAt) - c) + D; // fm
  double R = Rvol + Rsurf;
  // CS
  double CS = PI*R*R*(1 - Bc/KEcm);
  return CS*0.01; // fm^2 -> barn
}


// Shen cross section. 
// the differences between Shen and Kox is below 30 MeV/n, where Shen works a little bit better!
// c was not given in Shen and is taken to be equal to Kox
double Element::GetShenCrossSection(int Ap, int Zp, double E /* kinetic energy [MeV/n] */) {
  // source/processes/hadronic/cross_sections/src/G4IonsShenCrossSection.cc 
  // I found an error in Shen and Cox in the GEANT 4 documentation! A factor 5 in R_2 and a different definition of the D parameter in Kox.
  double At = GetA();
  double Zt = GetZ();
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  // Coulomb barrier 
  double Rt = 1.12*cubicAt - 0.94/cubicAt; // fm
  double Rp = 1.12*cubicAp - 0.94/cubicAp; // fm
  double r = Rt + Rp + 3.2; // fm
  double b = 1.0; // MeV/fm
  double KEcm = GetKEcm(Ap,Zp,E);
  double Bc = 1.44*Zt*Zp/r - b*Rt*Rp/(Rt+Rp); // MeV  
  // Coulomb barrier could prevent the reaction
  if (KEcm<=Bc) return 0;
  // Radius
  double alpha = 1; // fm
  double beta = 0.176; // MeV^1/3 fm
  double r0 = 1.1; // fm
  double c = GetC(E);
  double R1 = r0*(cubicAt + cubicAp  + 1.85*cubicAt*cubicAp/(cubicAt + cubicAp) - c);
  double R2;
  R2 = alpha*(At - 2*Zt)*Zp/(Ap*At);
  // quite an unjustified solution in order to symmetrize the cross section between projectile and target.
  // if(At>=Ap) R2 = alpha*(At - 2*Zt)*Zp/(Ap*At);
  // else       R2 = alpha*(Ap - 2*Zp)*Zt/(Ap*At);
  double R3 = beta/pow(KEcm,1./3.)*cubicAt*cubicAp/(cubicAt + cubicAp);
  double R = R1 + R2 + R3;
  double CS = PI*R*R*(1 - Bc/KEcm);
  return CS*0.01; // fm^2 -> barn
}


// OK, but this is not good for light nuclei. 
// Moreover, on Geant they speak about a validity range in enrgy. 
double Element::GetTripathiCrossSection(int Ap, int Zp, double E) {
  // Nuclear Instruments and Methods in Physics Research B 117 (1996) 347-349
  // "Universal parameterization of absorption cross sections"
  // R. K. Tripathi, Francis A. Cucinotta, John W. Wilson
  double At = GetA();
  double Zt = GetZ();
  double Ac = 12;
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  // double cubicAc = pow(Ac,1./3.);
  // Coulomb barrier
  double KEcm = GetKEcm(Ap,Zp,E);
  double r_rms_p = GetRadius(Zp,Ap);
  double r_rms_t = GetRadius();
  double r_rms_c = GetRadius(6,12); 
  double rp = 1.29*r_rms_p;
  double rt = 1.29*r_rms_t;
  double rc = 1.29*r_rms_c; 
  double R = rp + rt + 1.2*(cubicAt+cubicAp)/pow(KEcm,1./3.);
  double Bc = 1.44*Zt*Zp/R;
  // Coulomb barrier could prevent the reaction
  if (KEcm<=Bc) return 0;
  double S = cubicAt*cubicAp/(cubicAt+cubicAp);
  double rho_p = Ap/(PI*rp*rp*rp*4./3.);
  double rho_t = At/(PI*rt*rt*rt*4./3.);
  double rho_c = Ac/(PI*rc*rc*rc*4./3.);
  double D = 1.75*(rho_t + rho_p)/(rho_c + rho_c);
  double C = D*(1-exp(-E/40)) - 0.292*exp(-E/792.)*cos(0.229*pow(E,0.453));
  double delta;
  delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(At - 2*Zt)*Zp/(At*Ap); 
  // quite an unjustified solution in order to symmetrize the cross section between projectile and target.
  // if(At>=Ap) delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(At - 2*Zt)*Zp/(At*Ap); 
  // else       delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(Ap - 2*Zp)*Zt/(At*Ap); 
  double r0 = 1.1; // fm
  double CS = PI*r0*r0*pow(cubicAt+cubicAp+delta,2)*(1 -Bc/KEcm);
  return CS*0.01; // fm^2 -> barn
}


// Includes in the description also light nuclei, both as target and as projectile. 
// Is it really reversible between target and projectile? To be studied. 
// For heavy nuclei it does return the familiar Tripathi cross section.
double Element::GetTripathiUnivCrossSection(int Ap, int Zp, double E) {
  // Accurate universal parameterization of absorption cross sections III: light systems
  // R.K. Tripathi, F.A. Cucinotta, J.W. Wilson
  // Nuclear Instruments and Methods in Physics Research B 155 (1999) 349±356  
  double At = GetA();
  double Zt = GetZ();
  double Ac = 12;
  // double cubicAc = pow(Ac,1./3.);
  double cubicAt = pow(At,1./3.);
  double cubicAp = pow(Ap,1./3.);
  // Coulomb barrier
  double KEcm = GetKEcm(Ap,Zp,E);
  double r_rms_p = GetRadius(Zp,Ap);
  double r_rms_t = GetRadius();
  double r_rms_c = GetRadius(6,12);
  double rp = 1.29*r_rms_p;
  double rt = 1.29*r_rms_t;
  double rc = 1.29*r_rms_c;
  double R = rp + rt + 1.2*(cubicAt+cubicAp)/pow(KEcm,1./3.);
  double Bc = 1.44*Zt*Zp/R;
  // Coulomb barrier could prevent the reaction
  if (KEcm<=Bc) return 0;
  double S = cubicAt*cubicAp/(cubicAt+cubicAp);
  // Set Rc to one (arbitrary choice). If tabulated it will be modified later.
  double Rc=1.;
  double rho_p = Ap/(PI*rp*rp*rp*4./3.);
  double rho_t = At/(PI*rt*rt*rt*4./3.);
  double rho_c = Ac/(PI*rc*rc*rc*4./3.);
  // default high-Z nuclei
  double D = 1.75*(rho_t + rho_p)/(rho_c + rho_c);
  double T1 = 40;  
  // up to know it was the usual Thripati
  // the only problem is that some Rc parameters are not listed. Whenever this is the case we put Rc=1.
  // the symmetrization of the cross sections is done by hand for each case.
  // p(n) + X systems
  if ( (Ap==1) && (At>1) )  { 
    D = 1.85 + 0.16/(1+exp(500-E)/200);
    if      (Zp==0) T1 = 18;
    else if (Zp==1) T1 = 23; 
    if      ( (At==2) && (Zt==1) ) Rc = 13.5;
    else if ( (At==3) && (Zt==2) ) Rc = 21;
    else if ( (At==4) && (Zt==2) ) Rc = 27;
    else if (Zt==3)                Rc = 2.2;
  }
  // X + p(n) systems
  else if ( (At==1) && (Ap>1) )  { 
    D = 1.85 + 0.16/(1+exp(500-E)/200);
    if      (Zt==0) T1 = 18;
    else if (Zt==1) T1 = 23; 
    if      ( (Ap==2) && (Zp==1) ) Rc = 13.5;
    else if ( (Ap==3) && (Zp==2) ) Rc = 21;
    else if ( (Ap==4) && (Zp==2) ) Rc = 27;
    else if (Zp==3)                Rc = 2.2;
  }
  // d + X systems
  else if ( (Ap==2) && (Zp==1) && (At>1) ){
    D = 1.65 + 0.1/(1+exp(500-E)/200);
    T1 = 23;
    if      ( (Zt==1) && (At==2) ) Rc = 13.5;
    else if ( (Zt==2) && (At==4) ) Rc = 13.5;
    else if (Zt==6)                Rc = 6.;
  }  
  // X + d systems 
  else if ( (At==2) && (Zt==1) && (Ap>1)){ 
    D = 1.65 + 0.1/(1+exp(500-E)/200); 
    T1 = 23; 
    if      ( (Zp==1) && (Ap==2) ) Rc = 13.5;
    else if ( (Zp==2) && (Ap==4) ) Rc = 13.5;
    else if (Zp==6)                Rc = 6.;
  } 
  // 3He + X systems
  else if ( (Ap==3) && (Zp==2) ) {
    D = 1.55; 
    T1 = 40;
  } 
  // X + 3He systems
  else if ( (At==3) && (Zt==2) ){
    D = 1.55; 
    T1 = 40;
  }
  // 4He + X systems
  else if ( (Ap==4) && (Zp==2) ) {
    if ( (At==4) && (Zt==2) ) {
      double G = 75;                         
      T1 = 40;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zt<=4){
      double G = 300;                         
      T1 = 40;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zt<=7){
      double G = 300;                         
      T1 = 25;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zt<=13){
      double G = 300;                         
      T1 = 25;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zt<=26){
      double G= 300;
      T1 = 40; 
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else {  
      static int maxerr = 0;
      if (maxerr<100) { 
        cout << "Element::GetTripathiUnivCrossSection-W 4He+X dataset non available." << endl;	
        maxerr++;
      }
      return 0.;
    }
  }
  // X + 4He systems
  else if ( (At==4) && (Zt==2) )  {
    if ( (Ap==4) && (Zp==2) ){
      double G = 75;                         
      T1 = 40;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zp<=4){
      double G = 300;                         
      T1 = 40;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zp<=7){
      double G = 300;                         
      T1 = 25;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zp<=13){
      double G = 300;                         
      T1 = 25;
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else if (Zp<=26){
      double G= 300;
      T1 = 40; 
      D = 2.77 - 8.0e-3*At + 1.8e-5*At*At - 0.8/(1 + exp((250.0 - E)/G));
    }
    else { 
      cout << "Element::GetTripathiUnivCrossSection-W X+4He dataset non available." << endl;
      return 0.;
    }
  }
  double C = D*(1-exp(-E/T1)) - 0.292*exp(-E/792.)*cos(0.229*pow(E,0.453));
  double delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(At - 2*Zt)*Zp/(At*Ap);
  // quite an unjustified solution in order to symmetrize the cross section between projectile and target.
  if (At>=Ap) delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(At - 2*Zt)*Zp/(At*Ap); 
  else        delta = 1.85*S + 0.16*S/pow(KEcm,1./3.) - C + 0.91*(Ap - 2*Zp)*Zt/(At*Ap); 
  double X1 = 2.83 - (0.031*At) + (0.00017*At*At);
  if ( ((Zp==0)&&(Zt==2)&&(At==4))||((Zt==0)&&(Zp==2)&&(Ap==4)) ) X1=5.2;
  double SL = 1.2 + 1.6*(1-exp(-E/15));
  double Xm = 1 - X1*exp(-E/(X1*SL));
  double r0 = 1.1; // fm
  double CS = PI*r0*r0*pow(cubicAt+cubicAp+delta,2)*(1 - Bc/KEcm);
  // usual Tripathi
  if((Zp>2)&&(Zt>2)) return CS*0.01; //fm^2->barn
  // light modified Tripathi
  return 0.01*PI*r0*r0*pow(cubicAt+cubicAp+delta,2)*(1 -Rc*Bc/KEcm)*Xm; // fm^2 -> barn
}

// Glauber-Gribov cross section
// treats the nucleus-nucleus cross section as the convolution of the nucleon-nucleon cross sections, assuming a gaussian nucleon distribution within the nuclei. 
// Therefore it requires the energy-dependent nucleon nucleon (total) cross sections. 
// This parametrization has the advantage to be manifestly symmetric with respect to target-projectile interchange.
// http://indico.cern.ch/getFile.py/access?contribId=4&sessionId=0&resId=0&materialId=slides&confId=186631
double Element::GetGlauberCrossSection(int Ap, int Zp, double E) {
  // define the target A and Z
  int Zt = GetZ();
  int At = GetA();
  // define neutron numbers
  int Nt = At-Zt;
  int Np = Ap-Zp;
  // center of mass kinetic energy energy
  double KEcm = GetKEcm(Ap,Zp,E);
  // nucleon-nucleon total cross sections in mbarn
  double sigma_pp;
  double sigma_nn;
  double sigma_pn;
  // Glauber-Gribov model paramteres
  double r=1.2; //this is the value that shoud be introduced somehow: they claim it must be in the range 1.1-1.28 fm and that it has to be taken from the fit.
  double Rt, Rp;
  // special case of H target or pojectile (??)
  // Can the Glauber-Gribov can be used for H target? If yes, how?
  if (At==1) Rt=1.1;
  if (Ap==1) Rp=1.1;
  if ((At>1)&&(At<=50)) Rt=r*pow(At,1./3.)*(1-pow(At,-2./3.));
  if ((Ap>1)&&(Ap<=50)) Rp=r*pow(Ap,1./3.)*(1-pow(Ap,-2./3.));
  if (Ap>50) Rt=1.*pow(At,0.27);
  if (At>50) Rp=1.*pow(Ap,0.27);
  double R_min = Rt+Rp;
  double Bc = 1.44*Zt*Zp/R_min;  //is it really so? CHECK!
  // Coulomb barrier could prevent the reaction
  if (KEcm<=Bc) return 0;
  // Nucleon nucleon energy-dependent cross section parametrization in mbarn
  // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20080014212_2008013699.pdf	
  // mass and momentum of the proton with lab kin energy E
  double m = 938.3;
  double p = sqrt(E*(E+2*m));
  if (E<0.1) {
    sigma_pp = exp(6.51*exp(-pow(E/134,0.7)));
    sigma_nn = sigma_pp;
    sigma_pn = 26000*exp(-pow(E/0.282,0.3));
  }
  if ( (E>=0.1)&&(E<25) ) {
    sigma_pp = exp(6.51*exp(-pow(E/134,0.7)));
    sigma_nn = sigma_pp;
    sigma_pn = 38 + 12500*exp(-1.187*pow((E-0.1),0.35));
  }
  if ( (E>=25)&&(p<500) ) {
    sigma_pp = (1+5/E)*(40 + 109*cos(0.199*sqrt(E))*exp(-0.451*pow((E-25),0.258)));
    sigma_nn = sigma_pp;
    sigma_pn = 38 + 12500*exp(-1.187*pow((E-0.1),0.35));
  }
  if ( (p>=500)&&(p<1800) ) {
    sigma_pp = (1 + 5/E)*(40+109*cos(0.199*sqrt(E))*exp(-0.451*pow((E-25),0.258)));
    sigma_nn = sigma_pp; 
    sigma_pn = 40 + 10*cos(0.00369*p-0.943)*exp(-0.00895741*pow(p,0.8)+2);
  }
  if ( (p>=1800)&&(p<2000) ) {
    sigma_pp = 158.547/(pow((E*(E+2*m)),0.08));
    sigma_nn = sigma_pp;
    sigma_pn = 40 + 10*cos(0.00369*p-0.943)*exp(-0.00895741*pow(p,0.8)+2);
  }
  if ( (p>=2000)&&(p<4700) ) {
    sigma_pp = 158.547/(pow((E*(E+2*m)),0.08));
    sigma_nn = sigma_pp;
    sigma_pn = 35.80 + 0.308*log(2*m*(2*m + E)/(5380*5380)) + 40.15*pow((1000000/(2*m*(2*m + E))),0.458) - 30.00*pow((1000000/(2*m*(2*m + E))),0.545);
  }
  if (p>=4700) {
    sigma_pp = 35.45 + 0.308*log(2*m*(2*m + E)/(5380*5380)) + 42.53*pow((1000000/(2*m*(2*m + E))),0.458) - 33.34*pow((1000000/(2*m*(2*m + E))),0.545);
    sigma_nn = sigma_pp;
    sigma_pn = 35.80 + 0.308*log(2*m*(2*m + E)/(5380*5380)) + 40.15*pow((1000000/(2*m*(2*m + E))),0.458) - 30.00*pow((1000000/(2*m*(2*m + E))),0.545);
  }
  // mbarn->barn
  sigma_pp = sigma_pp/1000;
  sigma_nn = sigma_nn/1000;
  sigma_pn = sigma_pn/1000;
  //special case p+p scattering
  if ( (At==1)&&(Ap==1) ) return sigma_pp;
  double area = PI*(Rp*Rp+Rt*Rt);
  double AAsigma = Zp*Zt*sigma_pp+(Zp*Nt+Zt*Np)*sigma_pn+Nt*Np*sigma_nn;
  // factor 100 for the conversion fm^2->barn
  double cross = area*(1 - Bc/KEcm)*log(1+100*AAsigma/(area));	
  // factor 0.01 for the conversion fm^2->barn	
  return cross*0.01;
}


// General method for getting the cross section, specifying its type by means of the type parameter.
double Element::GetCrossSection(int type, int Ap, int Zp, double E /* kinetic energy [MeV/n] */) {
  switch (type) {
  case 0:
    return GetIntCrossSection();
    break;
  case 1:
    return GetColCrossSection();
    break;
  case 2:
    return GetGeometricCrossSection(Ap);
    break;
  case 3:
    return GetBradtPetersCrossSection(Ap);
    break;
  case 4:
    return GetSihverCrossSection(Ap);
    break;
  case 5:
    return GetKoxCrossSection(Ap,Zp,E);
    break;
  case 6:
    return GetShenCrossSection(Ap,Zp,E);
    break;    
  case 7:
    return GetTripathiCrossSection(Ap,Zp,E);
    break;    
  case 8:
    return GetTripathiUnivCrossSection(Ap,Zp,E);
    break; 
  case 9:
    return GetGlauberCrossSection(Ap,Zp,E);
    break;     
  }
  return 0.;
}


double Element::GetMeanFreePath(int type, int Ap, int Zp, double E) {
  return GetA()/(GetCrossSection(type,Ap,Zp,E)*0.602214129);
}


// returns the survival probability giving the thickness and density of the target element.
double Element::GetSurvivalProb(double x, double rho, int type, int Ap, int Zp, double E) {
  return exp(-x*rho/GetMeanFreePath(type,Ap,Zp,E));
}


///////////////////////////////////
// Energy Loss
///////////////////////////////////


void Element::ReadStoppingPowerFile(char* filename) {
  FILE* file = fopen(filename,"r");
  if (file==0) {
    printf("Error opening file %s\n",filename);
    return;
  }
  int i=0;
  while (!(feof(file))) {
    float a,b,c;
    fscanf(file,"%E %E %E",&a,&b,&c);
    Energy[i] = a;
    dEdx[i]   = b;
    Range[i]  = c;
    i++;
    if (i>=DEDXTABLE) return;
  }
}


int Element::GetLowerIndex(double energy) {
  if (energy<0) return -1;
  for (int i=0; i<DEDXTABLE; i++) 
    if (log10(energy)<=log10(Energy[i])) 
      return i;
  return -1;
}


int Element::GetUpperIndex(double energy) {
  if (energy<0) return -1;
  for (int i=DEDXTABLE-1; i>=0; i--)
    if (log10(energy)>=log10(Energy[i]))
      return i;
  return -1;
}


double Element::GetMIP() {
  double mip = 10000.;
  for (int i=0; i<DEDXTABLE; i++)
    if (mip>dEdx[i]) mip = dEdx[i];
  return mip;
}


double Element::GetRangeFromInterpolation(double energy) {
  int iu = GetUpperIndex(energy);
  int il = GetLowerIndex(energy);
  if (iu==il) return Range[iu];
  double width = log10(Energy[iu])-log10(Energy[il]);
  double wu    = (log10(Energy[iu])-log10(energy))/width;
  double wl    = (log10(energy)-log10(Energy[il]))/width;
  double ru    = Range[iu];
  double rl    = Range[il];
  return wl*ru + wu*rl;
}


double Element::GetRangeFromdEdxIntegral(double energy) {
  int lastbin = min(GetLowerIndex(energy),DEDXTABLE);
  if (lastbin==0) return 0;
  // first bin 
  double de = Energy[0];
  double invdedx = 1./dEdx[0];
  double s = invdedx*de;
  // other bins
  for (int i=1; i<=lastbin; i++) {
    de = Energy[i] - Energy[i-1];
    invdedx = (1./dEdx[i] + 1./dEdx[i-1])/2.; // regola del trapezio
    s += invdedx*de;
  }
  // last bin
  de = energy - Energy[lastbin];
  invdedx = (1./GetdEdx(energy) + 1./dEdx[lastbin])/2.;
  s += invdedx*de;
  return s;
}


double Element::GetdEdx(double energy) {
  int iu = GetUpperIndex(energy);
  int il = GetLowerIndex(energy);
  if (iu==il) return dEdx[iu];
  double width = log10(Energy[iu])-log10(Energy[il]);
  double wu    = (log10(Energy[iu])-log10(energy))/width;
  double wl    = (log10(energy)-log10(Energy[il]))/width;
  double eu    = dEdx[iu];
  double el    = dEdx[il];
  return wl*eu + wu*el;
}


double Element::GetKineticEnergyForRange(double range, int type) {
  for (int i=0; i<DEDXTABLE; i++) {
    double energy1 = Energy[i];
    double energy2 = Energy[i+1];
    double range1  = GetRange(energy1,type);
    double range2  = GetRange(energy2,type);
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
  
   
    
