{
  char name[100];
  float nb[31], nc[31], rbc[31];

  TFile *fout = new TFile("OutSysErr.root","RECREATE");

  TTree *stree = new TTree("stree","stree");

  stree->Branch("nb",&nb,"nb[31]/F");
  stree->Branch("nc",&nc,"nc[31]/F");
  stree->Branch("rbc",&rbc,"rbc[31]/F");

  TH1D *hebrig, *hecrig, *hebcrig;
  for(int i=0; i<500; i++)
  {
    sprintf(name,"output%d.root",i);
    cout<<"Reading file: "<<name<<endl;
    TFile fin(name,"READ");
    if(fin.IsZombie()) continue;
    hebrig = (TH1D*)fin.Get("hebrig");
    hecrig = (TH1D*)fin.Get("hecrig");
    hebcrig = (TH1D*)fin.Get("hebcrig");
    for(int j=1; j<=31; j++)
    {
      nb[j] = hebrig->GetBinContent(j);
      nc[j] = hecrig->GetBinContent(j);
      rbc[j] = hebcrig->GetBinContent(j);
    }
    stree->Fill();
  }

  fout->Write();
  fout->Close();
}
