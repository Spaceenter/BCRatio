void SetVisibilitiesAndColors(TGeoVolume* volume) { 
  TString volume_name(volume->GetName());
  TString material_name(volume->GetMaterial()->GetName());
  bool visibility=kTRUE;
  if ( (material_name.Contains("VACUUM"))|| (volume_name.Contains("ECFW"))||(volume_name.Contains("PMTB")) ) {
    volume->SetVisibility(kFALSE);
  }
  else if (material_name.Contains("HONEY")) {
    volume->SetLineColor(kGray);
    volume->SetTransparency(10);
    volume->SetVisibility(visibility);
  }
  else if (material_name.Contains("MAGNET")) {
    volume->SetLineColor(kRed+1);
    volume->SetTransparency(10);
    volume->SetVisibility(visibility);
  }
  else if ( 
    (material_name.Contains("TOF"))||
    (material_name.Contains("MYLAR"))||
    (material_name.Contains("SCINT")) ) {
    volume->SetLineColor(kBlue);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else if ( (material_name.Contains("SILICON")) ) {
    volume->SetLineColor(kViolet);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else if ( (material_name.Contains("CARBON"))||(material_name.Contains("FOAM")) ) {
    volume->SetLineColor(kGray+1);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else if (material_name.Contains("TRD")) {
    volume->SetLineColor(kGreen+1);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else if (material_name.Contains("EC")) {
    volume->SetLineColor(kOrange+1);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else if (material_name.Contains("RICH")) {
    volume->SetLineColor(kWhite);
    volume->SetTransparency(5);
    volume->SetVisibility(visibility);
  }
  else {
    volume->SetLineColor(kWhite);
    volume->SetTransparency(1);
    volume->SetVisibility(); //kFALSE);
    // cout << volume_name << " " << material_name << endl;
  }
  for (int ii=0; ii<volume->GetNdaughters(); ii++) {
    TGeoVolume* dau = volume->GetNode(ii)->GetVolume();
    SetVisibilitiesAndColors(dau);
  }
}


void DrawGeometry() {
  TEveManager::Create();
  if (gGeoManager==0) TGeoManager::Import("data/ams02_dec2011.root");
  TGeoVolume* top = gGeoManager->GetTopVolume();
  TGeoNode* node1 = top->FindNode("FMOT_00002_2");
  TEveGeoTopNode* inn = new TEveGeoTopNode(gGeoManager, node1);
  gEve->AddGlobalElement(inn);
  SetVisibilitiesAndColors(top);
  gEve->Redraw3D(kTRUE);
  TGLViewer *v = gEve->GetDefaultGLViewer();
  v->DoDraw();
}


void DrawOneTrack(Double_t point[3],Double_t dir[3]) {
  TEveTrackList *list = new TEveTrackList();
  TEveTrackPropagator* prop = list->GetPropagator();
  prop->SetMaxZ(200.);
  TEveTrack *track = 0;
  prop->SetMagField(0.);
  TEveRecTrackD *rc = new TEveRecTrackD();
  rc->fV.Set(point[0],point[1],point[2]);
  rc->fP.Set(dir[0],dir[1],dir[2]);
  TEveTrack* track = new TEveTrack(rc,prop);
  gEve->AddElement(list);
  list->AddElement(track);
  track->SetLineColor(kMagenta+1);
  track->MakeTrack();
  gEve->Redraw3D(kTRUE);
}

