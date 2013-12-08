{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("purity.root");
  int HMax;

  //B: 4.066 - 5.326 GV

  RooPlot *h1 = (RooPlot*)f.Get("qframeB_5");
  h1->GetXaxis()->SetTitleSize(0.047);
  h1->GetYaxis()->SetTitleSize(0.047);
  h1->SetTitleSize(0.047);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetLogy();
//  c1->SetGridx();
//  c1->SetGridy();
  h1->Draw();
  h1->GetXaxis()->SetRangeUser(3.5,9.5);
  h1->GetXaxis()->SetTitle("Tracker Layer 1 Charge");
  h1->GetYaxis()->SetTitle("Entries");
  HMax=h1->GetMaximum();
  TLine *ll1a = new TLine(4,0,4,HMax);
  ll1a->SetLineWidth(2);
  ll1a->SetLineColor(6);
  ll1a->Draw("same");
  TLine *ll1b = new TLine(5.5,0,5.5,HMax);
  ll1b->SetLineWidth(2);
  ll1b->SetLineColor(6);
  ll1b->Draw("same");

  //B: 139.226 - 210.249 GV

  RooPlot *h2 = (RooPlot*)f.Get("qframeB_17");
  h2->GetXaxis()->SetTitleSize(0.047);
  h2->GetYaxis()->SetTitleSize(0.047);
  h2->SetTitleSize(0.047);

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetLogy();
//  c2->SetGridx();
//  c2->SetGridy();
  h2->Draw();
  h2->GetXaxis()->SetRangeUser(3.5,9.5);
  h2->GetXaxis()->SetTitle("Tracker Layer 1 Charge");
  h2->GetYaxis()->SetTitle("Entries");
  HMax=h2->GetMaximum();
  TLine *ll2a = new TLine(4,0,4,HMax);
  ll2a->SetLineWidth(2);
  ll2a->SetLineColor(6);
  ll2a->Draw("same");
  TLine *ll2b = new TLine(5.5,0,5.5,HMax);
  ll2b->SetLineWidth(2);
  ll2b->SetLineColor(6);
  ll2b->Draw("same");

  //C: 4.066 - 5.326 GV

  RooPlot *h3 = (RooPlot*)f.Get("qframeC_5");
  h3->GetXaxis()->SetTitleSize(0.047);
  h3->GetYaxis()->SetTitleSize(0.047);
  h3->SetTitleSize(0.047);

  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->cd();
  c3->SetLogy();
//  c3->SetGridx();
//  c3->SetGridy();
  h3->Draw();
  h3->GetXaxis()->SetRangeUser(3.5,9.5);
  h3->GetXaxis()->SetTitle("Tracker Layer 1 Charge");
  h3->GetYaxis()->SetTitle("Entries");
  HMax=h3->GetMaximum();
  TLine *ll3a = new TLine(5,0,5,HMax);
  ll3a->SetLineWidth(2);
  ll3a->SetLineColor(6);
  ll3a->Draw("same");
  TLine *ll3b = new TLine(6.8,0,6.8,HMax);
  ll3b->SetLineWidth(2);
  ll3b->SetLineColor(6);
  ll3b->Draw("same");

  //C: 139.226 - 210.249 GV

  RooPlot *h4 = (RooPlot*)f.Get("qframeC_17");
  h4->GetXaxis()->SetTitleSize(0.047);
  h4->GetYaxis()->SetTitleSize(0.047);
  h4->SetTitleSize(0.047);

  TCanvas *c4 = new TCanvas("c4","c4",1000,600);
  c4->cd();
  c4->SetLogy();
//  c4->SetGridx();
//  c4->SetGridy();
  h4->Draw();
  h4->GetXaxis()->SetRangeUser(3.5,9.5);
  h4->GetXaxis()->SetTitle("Tracker Layer 1 Charge");
  h4->GetYaxis()->SetTitle("Entries");
  HMax=h4->GetMaximum();
  TLine *ll4a = new TLine(5,0,5,HMax);
  ll4a->SetLineWidth(2);
  ll4a->SetLineColor(6);
  ll4a->Draw("same");
  TLine *ll4b = new TLine(6.8,0,6.8,HMax);
  ll4b->SetLineWidth(2);
  ll4b->SetLineColor(6);
  ll4b->Draw("same");
}
