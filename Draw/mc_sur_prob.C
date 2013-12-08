{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("mc_sur_prob.root");

  RooPlot *rp = (RooPlot*)f.Get("qframe_2_10");

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetLogy();
  rp->Draw();
  rp->GetXaxis()->SetTitle("MC Inner Tracker Charge Estimator");
  rp->GetYaxis()->SetTitle("Entries");
  rp->GetXaxis()->SetRangeUser(0,6.5);
  TPaveText *t1 = new TPaveText(2.3,300,5,900);
  t1->AddText("^{12}C: 10.0GV - 12.6GV");
  t1->Draw("same");

  TH1D *h1 = (TH1D*)f.Get("hsurpr10");
  h1->SetLineWidth(2);
  h1->SetLineColor(2);
  TH1D *h2 = (TH1D*)f.Get("hsurpr11");
  h2->SetLineWidth(2);
  h2->SetLineStyle(2);

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  h1->Draw();
  h2->Draw("same");
  h1->GetXaxis()->SetTitle("Log_{10}(Rigidity)");
  h1->GetYaxis()->SetTitle("Survival Probability Ratio");
  h1->GetXaxis()->SetRangeUser(0.8,1.89);
  h1->GetYaxis()->SetRangeUser(0.8,1.2);
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(h1,"^{10}B/^{12}C");
  leg->AddEntry(h2,"^{11}B/^{12}C");
  leg->Draw("same");
}
