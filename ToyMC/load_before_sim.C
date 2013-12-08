{
  // classes used to estimate material properties
  gROOT->ProcessLine(".L Element.cxx+");
  gROOT->ProcessLine(".L Material.cxx+");
  gROOT->ProcessLine(".L Detector.cxx+");

  // load AMS geometry
  TGeoManager::Import("data/ams02_dec2011.root");
}
