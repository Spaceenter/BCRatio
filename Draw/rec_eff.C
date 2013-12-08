{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  //TOF 

  TFile fm("int_mc_dpmjet.root");

  TH1D *hb = (TH1D*)fm.Get("hetof_1");
  hb->SetLineWidth(2);
  hb->SetLineColor(2);
  TH1D *hc = (TH1D*)fm.Get("hetof_2");
  hc->SetLineWidth(2);
  hc->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  hb->Draw();
  hb->GetXaxis()->SetTitle("Log10(Rigidity)");
  hb->GetYaxis()->SetTitle("TOF Reconstruction Efficiency");
  hb->GetXaxis()->SetTitleSize(0.047);
  hb->GetYaxis()->SetTitleSize(0.047);
  hb->GetYaxis()->SetRangeUser(0.95,1.05);
  hc->Draw("same");
  TLegend *l1 = new TLegend(0.6,0.6,0.9,0.9);
  l1->AddEntry(hb,"Boron");
  l1->AddEntry(hc,"Carbon");
  l1->Draw("same");

  //Tracker

  TFile fd("sur_prob.root");

  TH1D *hbt = (TH1D*)fd.Get("hetr_1");
  hbt->SetLineWidth(2);
  hbt->SetLineColor(2);
  TH1D *hct = (TH1D*)fd.Get("hetr_3");
  hct->SetLineWidth(2);
  hct->SetLineStyle(2);

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  hbt->Draw();
  hbt->GetXaxis()->SetTitle("Log10(Rigidity)");
  hbt->GetYaxis()->SetTitle("Tracker Reconstruction Efficiency");
  hbt->GetXaxis()->SetTitleSize(0.047);
  hbt->GetYaxis()->SetTitleSize(0.047);
  hbt->GetYaxis()->SetRangeUser(0.8,1.2);
  hct->Draw("same");
  TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
  l2->AddEntry(hbt,"Boron");
  l2->AddEntry(hct,"Carbon");
  l2->Draw("same");
}
