#include "Element.h"
#include "Material.h"
#include "Detector.h"


static Double_t amu  = 0.93146; // GeV/c2
static Double_t mass[28] = { 1.0079,  4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797,
  22.9897, 24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078,
  44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.693};
static Double_t z[28] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};


static Int_t element_counter = 0;
Element* TGeoElement2Element(TGeoElement* element) {
  Element* elem = new Element(Form("%s_ELEM%05d",element->GetName(),element_counter),element->A(),element->Z(),10,10,10);
  element_counter++;
  return elem;
}


static Int_t material_counter = 0;
Material* TGeoMaterial2Material(TGeoMaterial* gmat) {
  if (gmat->InheritsFrom("TGeoMixture")) { // TGeoMixture
    TGeoMixture* gmix = (TGeoMixture*) gmat;
    Double_t* z = gmix->GetZmixt();
    Double_t* a = gmix->GetAmixt();
    Double_t* w = gmix->GetWmixt();
    Material* material = new Material(Form("%s_MAT%05d",gmix->GetName(),material_counter),gmix->GetDensity());
    for (Int_t ielem=0; ielem<gmix->GetNelements(); ielem++) {
      Element* element = TGeoElement2Element(gmix->GetElement(ielem));
      material->AddElement(element,w[ielem]);
    }
  }
  else { // TGeoMaterial
    Material* material = new Material(Form("%s_MAT%05d",gmat->GetName(),material_counter),gmat->GetDensity());
    Element* element = TGeoElement2Element(gmat->GetElement());
    material->AddElement(element,1);
  }
  material_counter++;
  return material;
}


void CleanUp(Detector* detector) {
  if (!detector) return;
  for (Int_t imat=0; imat<detector->GetNMaterial(); imat++) {
    Material* material = detector->GetMaterial(imat);
    if (!material) return; 
    for (Int_t ielem=0; ielem<material->GetNElement(); ielem++) {
      Element* element = material->GetElement(ielem);
      if (!element) return;
      delete element;
    }
    delete material;
  }
}


Detector* CreateDetectorOnTrajectory(Double_t point[3], Double_t dir[3], Double_t z_limit, Bool_t verbose) {
  // if needed reload geometry 
  if (gGeoManager==0) TGeoManager::Import("data/ams02_dec2011.root");
  // initialization of trajectory and global vars
  gGeoManager->InitTrack(point,dir);
  Int_t index = 0;
  Double_t tot_grammage = 0.;
  Double_t tot_radlen = 0.;
  Detector* detector = new Detector("detector");
  Double_t this_point[3] = {point[0],point[1],point[2]};
  Double_t this_dir[3] = {dir[0],dir[1],dir[2]}; 
  // propagate 
  while ( (!gGeoManager->IsOutside())&&(this_point[2]>z_limit) ) {
    // current
    TGeoNode *cnode = gGeoManager->GetCurrentNode();
    TGeoVolume *cvol = gGeoManager->GetCurrentVolume();
    TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();
    for (Int_t i=0; i<3; i++) {
      this_point[i] = gGeoManager->GetCurrentPoint()[i];
      this_dir[i] = gGeoManager->GetCurrentDirection()[i];
    }
    // calculate step up to next boundary
    gGeoManager->FindNextBoundary();
    Double_t snext = gGeoManager->GetStep();
    Double_t this_grammage = snext*cmat->GetDensity();
    Double_t this_radlen = snext/cmat->GetRadLen();
    tot_grammage += this_grammage;
    tot_radlen += this_radlen;
    // step up to z_limit 
    Double_t last_point[3] = {
      this_point[0] + this_dir[0]/this_dir[2]*(z_limit - this_point[2]), 
      this_point[1] + this_dir[1]/this_dir[2]*(z_limit - this_point[2]),
      z_limit
    };
    Double_t z_end_point = this_point[2] - snext*sqrt(1/(1 + pow(this_dir[0]/this_dir[2],2) + pow(this_dir[1]/this_dir[2],2)));
    Double_t last_step = sqrt(pow(last_point[0]-this_point[0],2)+pow(last_point[1]-this_point[1],2)+pow(last_point[2]-this_point[2],2));
    Double_t step = (z_end_point<z_limit) ? last_step: snext; //(z_final<z_limit) ? last_step : snext;
    if (verbose)
      printf("%3d %10.2f %10.4f %10.4f %20s %20s %20s %8.3f %8.3f %8.3f %8.3g %8.3f %8.3f %8.3g %8.3g %3d %3d %3d \n",
        index,this_point[0],this_point[1],this_point[2],
        cnode->GetName(),cvol->GetName(),cmat->GetName(),
        gGeoManager->GetStep(),last_step,
        cmat->GetDensity(),cmat->GetRadLen(),this_grammage,tot_grammage,100*this_radlen,100*tot_radlen,
        gGeoManager->IsEntering(),gGeoManager->IsOnBoundary(),gGeoManager->IsOutside()
      );
    // do not add vacuum slabs
    if (cmat->GetDensity()>1e-10) detector->AddMaterial(TGeoMaterial2Material(cmat),step);
    // move on 
    TGeoNode *newNode = gGeoManager->Step();
    for (Int_t i=0; i<3; i++) {
      this_point[i] = gGeoManager->GetCurrentPoint()[i];
      this_dir[i] = gGeoManager->GetCurrentDirection()[i];
    }
    index++;
  }
  if (verbose) detector->Info(2);
  return detector;
}


// check the detector construction for a given trajectory 
void OneTrajectory(Double_t point[3], Double_t dir[3]) {
  // Double_t point[3] = {2,2,200};
  // Double_t dir[3] = {0,0.19,-0.9};
  //   70  Under TRD 
  //   60  Under UTOF 
  //  -75  Under LTOF 
  // -138  Over ECAL 
  Double_t z_limit=-138;
  Double_t grammage;
  Double_t radlen;
  Double_t survival[10];
  Detector* detector = CreateDetectorOnTrajectory(point,dir,z_limit,true);
  if (!detector) continue;
  CleanUp(detector);
  if (detector) delete detector;
}


// coordinates of the edges of each tracker plane
static Double_t tracker_layers_z[9] = {158.920,53.060,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882};
static Double_t tracker_planes_edges[9][4] = {
  {-62.14,  -47.40,   62.14,   47.40},
  {-62.14,  -40.10,   62.14,   40.10},
  {-49.70,  -43.75,   49.70,   43.75},
  {-49.72,  -43.75,   49.72,   43.75},
  {-49.71,  -36.45,   49.70,   36.45},
  {-49.72,  -36.45,   49.72,   36.45},
  {-49.72,  -43.75,   49.71,   43.75},
  {-49.72,  -43.75,   49.71,   43.75},
  {-45.62,  -29.48,   45.55,   29.53}
};
// check if a trajectory passes through tracker
Int_t GetPatternInsideTracker(Double_t point[3], Double_t dir[3]) {
  Int_t pattern = 0;
  for (int ilayer=0; ilayer<9; ilayer++) {
    Double_t tanxz = dir[0]/dir[2];
    Double_t tanyz = dir[1]/dir[2];
    Double_t x = point[0] + tanxz*(tracker_layers_z[ilayer] - point[2]);
    Double_t y = point[1] + tanyz*(tracker_layers_z[ilayer] - point[2]);
    Bool_t isinlayer = false;
    if ( (x>tracker_planes_edges[ilayer][0])&&
         (x<tracker_planes_edges[ilayer][2])&&
         (y>tracker_planes_edges[ilayer][1])&&
         (y<tracker_planes_edges[ilayer][3]) ) {
      if ((ilayer+1)==9) isinlayer = true;
      else {
        if ( (sqrt(x*x+y*y)<tracker_planes_edges[ilayer][2]) ) isinlayer = true;
      }
    }
    if (isinlayer) pattern |= (1<<ilayer);
  }
  return pattern;
}


// circle is a simple approximation
static float trd_layers_z[2] = {158,79}; // ad occhio
static float trd_radius[2] = {110,77}; // ad occhio
// check if a trajectory passes through TRD
int GetPatternInsideTRD(double point[3], double dir[3]) {
  int pattern = 0;
  // direction
  double tanxz = dir[0]/dir[2];
  double tanyz = dir[1]/dir[2];
  float x,y;
  // bottom
  x = point[0] + tanxz*(trd_layers_z[1] - point[2]);
  y = point[1] + tanyz*(trd_layers_z[1] - point[2]);
  bool passTrdBot = (x*x+y*y)<(trd_radius[1]*trd_radius[1]); 
  // top 
  x = point[0] + tanxz*(trd_layers_z[0] - point[2]);
  y = point[1] + tanyz*(trd_layers_z[0] - point[2]);
  bool passTrdTop = (x*x+y*y)<(trd_radius[1]*trd_radius[1]);
  // return value
  if (passTrdBot) pattern |= 0x1;
  if (passTrdTop) pattern |= 0x2; 
  return pattern;
}


// random extraction of a Trajectory with check
bool GenerateTrajectory(Double_t point[3], Double_t dir[3]) {
  // upper face of a cube of 380 cm side centrated in the AMS02 origin
  point[0] = -195.*2 + 2*2*195.*gRandom->Rndm();
  point[1] = -195.*2 + 2*2*195.*gRandom->Rndm();
  point[2] = 195.;
  // angular isotropic distribution (sin * cos)
  Double_t theta = acos(sqrt(gRandom->Rndm())); 
  Double_t phi = 2.*TMath::Pi()*gRandom->Rndm();
  dir[0] = sin(theta)*cos(phi);
  dir[1] = sin(theta)*sin(phi);
  dir[2] = -cos(theta);
  return true;
} 


Double_t GenerateEnergy(Double_t Emin = 10, Double_t Emax = 1000000) {
  return 0;  
}


void DisplayManyTrackExample() {
  DrawGeometry();
  Double_t point[3];
  Double_t dir[3];
  for (Int_t ievent=0; ievent<1000000; ievent++) {
    if (!GenerateTrajectory(point,dir)) continue;
    // if ((GetPatternInsideTracker(point,dir)&0xff)!=0xff) continue;
    if (GetPatternInsideTRD(point,dir)!=0) continue;
    if ((GetPatternInsideTracker(point,dir)&0xfe)!=0xfe) continue;
    DrawOneTrack(point,dir);
  }
}


void HeightStudyMonteCarlo(Int_t nevent, Int_t Z, Int_t A, Double_t E) {
  // prepare output
  TFile* file = TFile::Open(Form("ntuple/height_inside_A%02d_Z%02d.root",A,Z),"recreate");
  TTree* ntuple = new TTree("Ntuple","Ntuple"); 
  Int_t event; ntuple->Branch("event",&event,"event/I");
  ntuple->Branch("A",&A,"A/I");
  ntuple->Branch("Z",&Z,"Z/I");
  ntuple->Branch("E",&E,"E/D");
  Double_t z_limit; ntuple->Branch("z_limit",&z_limit,"z_limit/D"); 
  Double_t point[3]; ntuple->Branch("point[3]",point,"point[3]/D");  
  Double_t dir[3]; ntuple->Branch("dir[3]",dir,"dir[3]/D");   
  Double_t grammage; ntuple->Branch("grammage",&grammage,"grammage/D");
  Double_t survival[10]; ntuple->Branch("survival[10]",survival,"survival[10]/D");
  for (event=0; event<nevent; event++) {
    // progress
    if (event%Int_t(nevent/10)==0) cout << "Generated " << event << " trajectories in L1+Inner of " << nevent << endl; 
    // extraction of a trajectory
    if (!GenerateTrajectory(point,dir)) continue; 
    // check if is inside L1+Inner
    if ((GetPatternInsideTracker(point,dir)&0xff)!=0xff) continue;
    // if (GetPatternInsideTRD(point,dir)!=0) continue;
    // if ((GetPatternInsideTracker(point,dir)&0xfe)!=0xfe) continue;
    // create materials up to a certain level 
    for (int iz=0; iz<180+138; iz++) {
      z_limit = 180 - iz;
      // create material stack 
      Detector* detector = CreateDetectorOnTrajectory(point,dir,z_limit,false);
      if (!detector) continue; 
      // calculated needed things 
      grammage = detector->GetGrammage();
      for (Int_t i= 0; i<10; i++) survival[i] = detector->GetSurvivalProb(i,0,A,Z,E);  
      ntuple->Fill();
      // clean
      CleanUp(detector);
      if (detector) delete detector;
    }
  }
  cout << "File " << file->GetName() << " generated." << endl; 
  ntuple->Write();
  file->Close();
}


void ProductionHeightStudyMonteCarlo() {
  // 10 GeV/n
  Int_t nevent = 10000; //100000;
  HeightStudyMonteCarlo(nevent,8,16,10000);
  HeightStudyMonteCarlo(nevent,6,12,10000);
  HeightStudyMonteCarlo(nevent,5,11,10000);
  HeightStudyMonteCarlo(nevent,5,10,10000);
}


void SurvivalStudyMonteCarlo(Int_t nevent, Int_t Z, Int_t A) {
  // prepare output
  TFile* file = TFile::Open(Form("ntuple/survival_A%02d_Z%02d.root",A,Z),"recreate");
  TTree* ntuple = new TTree("Ntuple","Ntuple");
  Int_t event; ntuple->Branch("event",&event,"event/I");
  ntuple->Branch("A",&A,"A/I");
  ntuple->Branch("Z",&Z,"Z/I");
  Double_t E; ntuple->Branch("E",&E,"E/D");
  Double_t z_limit = -75; ntuple->Branch("z_limit",&z_limit,"z_limit/D");
  Double_t point[3]; ntuple->Branch("point[3]",point,"point[3]/D");
  Double_t dir[3]; ntuple->Branch("dir[3]",dir,"dir[3]/D");
  Double_t grammage; ntuple->Branch("grammage",&grammage,"grammage/D");
  Double_t survival[10]; ntuple->Branch("survival[10]",survival,"survival[10]/D");
  for (event=0; event<nevent; event++) {
    // progress
    if (event%Int_t(nevent/10)==0) cout << "Generated " << event << " trajectories in L1+Inner of " << nevent << endl;
    // extraction of a trajectory
    if (!GenerateTrajectory(point,dir)) continue;
    // check if is inside Inner
    if ((GetPatternInsideTracker(point,dir)&0xfe)!=0xfe) continue;
    // create material stack 
    Detector* detector = CreateDetectorOnTrajectory(point,dir,z_limit,false);
    if (!detector) continue;
    // calculated needed things 
    grammage = detector->GetGrammage();
    for (int ie=0; ie<101; ie++) {
      E = 100.*pow(10.,(log10(100000.)-log10(100.))*ie/100.);
      for (Int_t i=0; i<10; i++) survival[i] = detector->GetSurvivalProb(i,0,A,Z,E);
      ntuple->Fill(); 
    }
    // clean
    CleanUp(detector);
    if (detector) delete detector;
  }
  cout << "File " << file->GetName() << " generated." << endl;
  ntuple->Write();
  file->Close();
}


void ProductionSurvivalStudyMonteCarlo() { 
  Int_t nevent = 10000; //100000;
  SurvivalStudyMonteCarlo(nevent,8,16);
  SurvivalStudyMonteCarlo(nevent,6,12);
  SurvivalStudyMonteCarlo(nevent,5,11);
  SurvivalStudyMonteCarlo(nevent,5,10);
}


void AllChargesStudyMonteCarlo() {
  // prepare output
  TFile* file = TFile::Open("ntuple/survival_allcharges_L1_upper.root","recreate");
  TTree* ntuple = new TTree("Ntuple","Ntuple");
  Int_t event; ntuple->Branch("event",&event,"event/I");
  Int_t A; ntuple->Branch("A",&A,"A/I");
  Int_t Z; ntuple->Branch("Z",&Z,"Z/I");
  Double_t E; ntuple->Branch("E",&E,"E/D");
  Double_t z_limit = 50; ntuple->Branch("z_limit",&z_limit,"z_limit/D");
  Double_t point[3]; ntuple->Branch("point[3]",point,"point[3]/D");
  Double_t dir[3]; ntuple->Branch("dir[3]",dir,"dir[3]/D");
  Double_t grammage; ntuple->Branch("grammage",&grammage,"grammage/D");
  Double_t survival[10]; ntuple->Branch("survival[10]",survival,"survival[10]/D");
  for (int iz=0; iz<28; iz++) {
    Z = z[iz];
    A = int(mass[iz]+0.5);
    E = 10000; // 10 GeV/n
    cout << A << " " << Z << endl;
    Int_t nevent = 1000000;
    for (event=0; event<nevent; event++) {
      // progress
      if (event%Int_t(nevent/10)==0) cout << "Generated " << event << " trajectories in L1+Inner of " << nevent << endl;
      // extraction of a trajectory
      if (!GenerateTrajectory(point,dir)) continue;
      // check if is inside L1+Inner
      if ((GetPatternInsideTracker(point,dir)&0xff)!=0xff) continue;
      // check if is inside Inner
      // if ((GetPatternInsideTracker(point,dir)&0xff)!=0xfe) continue;
      // create material stack 
      Detector* detector = CreateDetectorOnTrajectory(point,dir,z_limit,false);
      if (!detector) continue;
      // calculated needed things 
      grammage = detector->GetGrammage();
      for (Int_t i=0; i<10; i++) survival[i] = detector->GetSurvivalProb(i,0,A,Z,E);
      ntuple->Fill();
      // clean
      CleanUp(detector);
      if (detector) delete detector;
    }
  }
  cout << "File " << file->GetName() << " generated." << endl;
  ntuple->Write();
  file->Close();
}

void Thesis_Theory() //survival probablity v.s. energy
{
  TFile* file = TFile::Open("ntuple/toymc_sur_prob.root","recreate");
  TTree* ntuple = new TTree("ntuple","ntuple"); 

  int event;
  int pid;
  double E;
  double z_limit = -75;
  double survival[10];

  ntuple->Branch("event",&event,"event/I");
  ntuple->Branch("pid",&pid,"pid/I");
  ntuple->Branch("E",&E,"E/D");
  ntuple->Branch("z_limit",&z_limit,"z_limit/D"); 
  ntuple->Branch("survival[10]",survival,"survival[10]/D");

  int A, Z;
  double point[3] = {5, 5, 195}; //TOF good path length
  double dir[3] = {0, 0, -1};

  for(pid=0; pid<3; pid++)
  {
    if(pid==0) {A = 10; Z = 5;}
    if(pid==1) {A = 11; Z = 5;}
    if(pid==2) {A = 12; Z = 6;}

    Detector *detector = CreateDetectorOnTrajectory(point, dir, z_limit, false);
    if(!detector) cout<<"Error in creating detector!"<<endl;
    for (int ie=0; ie<101; ie++)
    {
      E = 5.*pow(10.,(log10(100000.)-log10(5.))*ie/100.);
      for (int i=0; i<10; i++) survival[i] = detector->GetSurvivalProb(i,0,A,Z,E);
      ntuple->Fill();
    }
    CleanUp(detector);
    if(detector) delete detector;
  }

  ntuple->Write();
  file->Close();
}

void Thesis_Material() //material v.s. depth
{
  TFile* file = TFile::Open("ntuple/toymc_material.root","recreate");
  TTree* ntuple = new TTree("ntuple","ntuple"); 

  int event;
  double z_limit;
  double grammage;

  ntuple->Branch("event",&event,"event/I");
  ntuple->Branch("z_limit",&z_limit,"z_limit/D"); 
  ntuple->Branch("grammage",&grammage,"grammage/D");

  double point[3] = {5, 5, 195}; //TOF good path length
  double dir[3] = {0, 0, -1};

  for (z_limit=180; z_limit>-138; z_limit--) 
  {
    Detector *detector = CreateDetectorOnTrajectory(point, dir, z_limit, false);
    if(!detector) cout<<"Error in creating detector!"<<endl;
    grammage = detector->GetGrammage();
    ntuple->Fill();
    CleanUp(detector);
    if(detector) delete detector;
  }

  ntuple->Write();
  file->Close();
}
