{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("sur_prob.root");

  TH1D *hs = (TH1D*)f.Get("hs0_0_10");
  hs->SetLineWidth(2);
  hs->SetLineColor(2);

  TH1D *ht = (TH1D*)f.Get("ht_0_10");
  ht->SetLineWidth(2);
  ht->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogy();
  hs->Draw();
  ht->Draw("same");
  hs->GetXaxis()->SetTitle("Tracker Layer 1 Charge");
  hs->GetYaxis()->SetTitle("Entries");
  hmax = hs->GetMaximum();
  TLine *l1 = new TLine(4.9,0,4.9,hmax);
  TLine *l2 = new TLine(5.1,0,5.1,hmax);
  TLine *l3 = new TLine(5.9,0,5.9,hmax);
  TLine *l4 = new TLine(6.1,0,6.1,hmax);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l4->SetLineWidth(2);
  l1->SetLineColor(4);
  l2->SetLineColor(4);
  l3->SetLineColor(3);
  l4->SetLineColor(3);
  l1->Draw("same");
  l2->Draw("same");
  l3->Draw("same");
  l4->Draw("same");
  TLegend *l = new TLegend(0.6,0.6,0.9,0.9);
  l->AddEntry(hs,"All Particles");
  l->AddEntry(ht,"Z = 1 or 2");
  l->Draw("same");

  TH1D *hq = (TH1D*)f.Get("hqin6_0_10");
  hq->SetLineWidth(2);
  hq->SetLineColor(2);

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  c2->SetLogy();
  hq->Draw();
  hq->GetXaxis()->SetTitle("Inner Tracker Charge");
  hq->GetYaxis()->SetTitle("Entries");

  TH1D *hsp = (TH1D*)f.Get("hmspr56_0");
  hsp->SetLineWidth(2);
  hsp->SetLineColor(2);

  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->cd();
  c3->SetGridx();
  c3->SetGridy();
  hsp->Draw();
  hsp->GetXaxis()->SetTitle("Log_{10}(Rigidity)");
  hsp->GetYaxis()->SetTitle("Survival Probability Ratio Between B and C");
  hsp->GetXaxis()->SetRangeUser(0.5,1.39);
  hsp->GetYaxis()->SetRangeUser(0.8,1.2);
}
