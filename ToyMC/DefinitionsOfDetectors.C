{
	
  int depth = 7;

  if (depth==0) {
    cout << "Select detector depth: " << endl;
    cout << " 1 - All AMS (no ECAL)" << endl;
    cout << " 2 - AMS - LUSS (no LToF+RICH+ECAL) " << endl;
    cout << " 3 - Upper AMS (1N-Tk + TRD + UToF + 1-Tk)" << endl; 
    cout << " 4 - Upper AMS no TRD" << endl;
    cout << " 5 - All AMS, no TRD, no RICH" << endl;
    cout << " 6 - AMS on Test Beam: 20m Air + 2 Scintillators + Upper AMS (1N-Tk + TRD + UToF + 1-Tk)" << endl;
    cout << " 7 - AMS (no ECAL, no RICH) " <<endl;
  }
   
  cout << "Defining AMS02 Detector ..." << endl;

  Detector* AMS02 = new Detector("AMS02_PM");


  if (depth==6) {
    AMS02->AddMaterial(Air,        200.);
    AMS02->AddMaterial(TofScint,     1.);
    AMS02->AddMaterial(TofScint,     1.); 
    cout << "Defining Air and Scintillators ..." << endl;
  }

  // Tk 1N
  AMS02->AddMaterial(TkSkin,         0.070);
  AMS02->AddMaterial(TkHoneycombExt, 3.86);
  AMS02->AddMaterial(TkSkin,         0.070);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSi,           0.03);
  cout << "Defining 1N-Tk ..." << endl;

  if ( (depth!=4)&&(depth!=5) ) {  
    // TRD
    AMS02->AddMaterial(Fleece,         2*20);
    AMS02->AddMaterial(TrdXe,          0.6*0.8*20);
    AMS02->AddMaterial(TrdCO2,         0.6*0.2*20);
    cout << "Defining TRD ..." << endl;
  }

  // TOF Upper
  AMS02->AddMaterial(TofHoneycomb,   10.);
  AMS02->AddMaterial(TofScint,       1.);
  AMS02->AddMaterial(TofScint,       1.);
  cout << "Defining Upper TOF ..." << endl;

  // Tk 1
  AMS02->AddMaterial(TkSi,           0.03);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSkin,         0.070);
  AMS02->AddMaterial(TkHoneycombExt, 3.86);
  AMS02->AddMaterial(TkSkin,         0.070);
  cout << "Defining 1-Tk ..." << endl;

  if ( (depth==3)||(depth==4)||(depth==6) ) return AMS02;

  // Tk 2 
  AMS02->AddMaterial(TkSi,           0.03);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkHoneycombInt, 1.24);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSi,           0.03);
  cout << "Defining 2-Tk ..." << endl;

  // Tk 3
  AMS02->AddMaterial(TkSi,           0.03);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkHoneycombInt, 1.24);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSi,           0.03);
  cout << "Defining 3-Tk ..." << endl;

  // Tk 4 
  AMS02->AddMaterial(TkSi,           0.03);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkHoneycombInt, 1.24);
  AMS02->AddMaterial(TkSkin,         0.022);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSi,           0.03);
  cout << "Defining 4-Tk ..." << endl;
  
  // Tk 5   //inert piece of material
  AMS02->AddMaterial(TkHoneycombInt, 1.24);
  cout << "Defining 5-Tk ..." << endl;

  if (depth==2) return AMS02;

  // TOF Lower
  AMS02->AddMaterial(TofScint,       1.);
  AMS02->AddMaterial(TofScint,       1.);
  AMS02->AddMaterial(TofHoneycomb,   10.);
  cout << "Defining Lower ToF ..." << endl;

  // RICH
  if ((depth!=5) && (depth!=7)) {
    AMS02->AddMaterial(RichRadiatorNaF,0.5);
    cout << "Defining RICH NaF ..." << endl;
  }
  if (depth!=7)
  // Tk 6
  AMS02->AddMaterial(TkSi,           0.03);
  AMS02->AddMaterial(TkFoam,         0.5);
  AMS02->AddMaterial(TkSkin,         0.070);
  AMS02->AddMaterial(TkHoneycombExt, 3.86);
  AMS02->AddMaterial(TkSkin,         0.070);
  cout << "Defining 6-Tk ..." << endl;
}

