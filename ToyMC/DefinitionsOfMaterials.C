//I modified it giving the thickness and the number of moduli of each material
//Actually this is not necessary, because one can define the class Detector which is given by a combination of Materials and their thickness

{
  //Define the fleece TRD component. There are 20 of them.
  Material* Fleece = new Material("Fleece",0.06); // C2 H4
  Fleece->AddElement(C,2);
  Fleece->AddElement(H,4);
  
  Material* TrdXe = new Material("TrdXe",0.006); 
  TrdXe->AddElement(Xe,1);
 
  //Substituted O with CO2!
  Material* TrdCO2 = new Material("TrdO",0.0015); 
  TrdCO2->AddElement(O,2);
  TrdCO2->AddElement(C,1);

  Material* TofHoneycomb= new Material("Tof_Honeycomb",0.085);
  TofHoneycomb->AddElement(Al,1);
    
  Material* TofScint = new Material("Polyvyniltoluene_up",1.03); // CH2 CH C6 H4 CH3 = C9 H10
  TofScint->AddElement(C, 9);
  TofScint->AddElement(H,10);
   
  Material* TkSi = new Material("Silicon", 2.33); 
  TkSi->AddElement(Si,1);
 
  Material* TkKapton = new Material("Kapton",1.47); // C22 H10 N2 05
  TkKapton->AddElement(C,22);
  TkKapton->AddElement(H,10);
  TkKapton->AddElement(N,2);
  TkKapton->AddElement(O,5);

  Material* TkFoam = new Material("TkFoam",0.1); 
  TkFoam->AddElement(C,1);

  Material* TkSkin = new Material("TkCarbonSkin",2.210); // da bill
  TkSkin->AddElement(C,1);

  Material* TkHoneycombInt = new Material("TkHoneycombInt",0.016); // da bill
  TkHoneycombInt->AddElement(Al,1);
  Material* TkHoneycombExt = new Material("TkHoneycombExt",0.032); // da bill
  TkHoneycombExt->AddElement(Al,1);
    
  Material* RichRadiatorNaF = new Material("RichNaF",2.558);
  RichRadiatorNaF->AddElement(Na,1);
  RichRadiatorNaF->AddElement(F,1);
    
  Material* Air = new Material("Air",1.2e-3); // specific gravity
  Air->AddElement(C, 0.000124);
  Air->AddElement(N, 0.755267);
  Air->AddElement(O, 0.231781);
  Air->AddElement(Ar,0.012827);

  Material* MLI = new Material("MLI",0);
  MLI->AddElement(C,5);
  MLI->AddElement(H,4);
  MLI->AddElement(O,2);

  Material* LEP = new Material("LEP",0);
  LEP->AddElement(C,1);

  Material* LeadFibers = new Material("LeadFibers",6.9);
  LeadFibers->AddElement(Pb,1);
}

