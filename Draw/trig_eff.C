{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("trig_eff.root");

  TH1D *ha = (TH1D*)f.Get("hetrigr");
  ha->SetLineWidth(2);
  ha->SetLineColor(2);

  TCanvas *c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetLogx();
  c->SetGridx();
  c->SetGridy();
  ha->Draw();
  ha->GetXaxis()->SetTitle("Rigidity (GV)");
  ha->GetYaxis()->SetTitle("Trigger Efficiency Ratio Between B and C");
}
