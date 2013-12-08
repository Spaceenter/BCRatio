{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("result.root");

  //eff. ratio for two methods

  TH1D *h1r = (TH1D*)f.Get("her1");
  h1r->SetLineWidth(2);
  h1r->SetLineColor(2);
  h1r->GetFunction("fer1")->SetBit(TF1::kNotDraw);
  
  TH1D *h2r = (TH1D*)f.Get("her2");
  h2r->SetLineWidth(2);
  h2r->SetLineStyle(2);

  TCanvas *cr = new TCanvas("cr","cr",1000,600);
  cr->cd();
  cr->SetLogx();
  cr->SetGridx();
  cr->SetGridy();
  h1r->Draw();
  h2r->Draw("same");
  h1r->GetXaxis()->SetTitle("Rigidity (GV)");
  h1r->GetYaxis()->SetTitle("Selection Efficiency Ratio Between B and C");
  TLegend *lr = new TLegend(0.6,0.6,0.9,0.9);
  lr->AddEntry(h1r,"Method One");
  lr->AddEntry(h2r,"Method Two");
  lr->Draw("same");

  //eff. for B and C using method 1

  TH1D *hb = (TH1D*)f.Get("he_1");
  hb->SetLineWidth(2);
  hb->SetLineColor(2);

  TH1D *hc = (TH1D*)f.Get("he_3");
  hc->SetLineWidth(2);
  hc->SetLineStyle(2);

  TCanvas *ce = new TCanvas("ce","ce",1000,600);
  ce->cd();
  ce->SetLogx();
  ce->SetGridx();
  ce->SetGridy();
  hb->Draw();
  hc->Draw("same");
  hb->GetXaxis()->SetTitle("Rigidity (GV)");
  hb->GetYaxis()->SetTitle("Selection Efficiency");
  TLegend *le = new TLegend(0.6,0.6,0.9,0.9);
  le->AddEntry(hb,"Boron");
  le->AddEntry(hc,"Carbon");
  le->Draw("same");

  //L1 XY efficiency

  TH1D *herp1 = (TH1D*)f.Get("herp1");
  herp1->SetLineWidth(2);
  herp1->SetLineColor(2);
  herp1->GetFunction("ferp1")->SetBit(TF1::kNotDraw);
  
  TCanvas *cr1 = new TCanvas("cr1","cr1",1000,600);
  cr1->cd();
  cr1->SetLogx();
  cr1->SetGridx();
  cr1->SetGridy();
  herp1->Draw();
  herp1->GetXaxis()->SetTitle("Rigidity (GV)");
  herp1->GetYaxis()->SetTitle("L1 Hit Pickup Efficiency Ratio Between B and C");
}
