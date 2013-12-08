{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("charge_id.root");

  //Z=3 to Z=8

  TH1D *ha = (TH1D*)f.Get("hqtrinall");
  ha->Scale(1.0/ha->Integral());
  ha->SetLineWidth(2);

  TCanvas *c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  ha->Draw();
  ha->GetXaxis()->SetTitle("Inner Tracker Charge");
  ha->GetYaxis()->SetTitle("Normalized Entries");

  TF1 *f3 = new TF1("f3","gaus",0,10);
  TF1 *f4 = new TF1("f4","gaus",0,10);
  TF1 *f5 = new TF1("f5","gaus",0,10);
  TF1 *f6 = new TF1("f6","gaus",0,10);
  TF1 *f7 = new TF1("f7","gaus",0,10);
  TF1 *f8 = new TF1("f8","gaus",0,10);

  ha->Fit("f3","","",2.8,3.3);
  ha->Fit("f4","","",3.8,4.3);
  ha->Fit("f5","","",4.7,5.3);
  ha->Fit("f6","","",5.7,6.3);
  ha->Fit("f7","","",6.7,7.3);
  ha->Fit("f8","","",7.7,8.3);

  f3->SetLineWidth(2);
  f4->SetLineWidth(2);
  f5->SetLineWidth(2);
  f6->SetLineWidth(2);
  f7->SetLineWidth(2);
  f8->SetLineWidth(2);

  f3->SetLineColor(9);
  f4->SetLineColor(3);
  f5->SetLineColor(4);
  f6->SetLineColor(2);
  f7->SetLineColor(6);
  f8->SetLineColor(7);

  f3->Draw("same");
  f4->Draw("same");
  f5->Draw("same");
  f7->Draw("same");
  f8->Draw("same");
  f6->Draw("same");

  TLine *line0 = new TLine(2.5,0,8.5,0);
  line0->SetLineWidth(2);
  line0->Draw("same");

  //Z=5 and Z=6
 
  TH1D *h6 = (TH1D*)f.Get("hqtrin6");
  double scale = h6->Integral();
  h6->Scale(1.0/scale);
  h6->SetLineWidth(2);
  h6->SetLineStyle(2);

  TH1D *h5 = (TH1D*)f.Get("hqtrin5");
  h5->Scale(1.0/scale);
  h5->SetLineWidth(2);
  h5->SetLineColor(2);

  TCanvas *cr = new TCanvas("cr","cr",1000,600);
  cr->cd();
  cr->SetLogy();
  cr->SetGridx();
  cr->SetGridy();
  h6->Draw();
  h5->Draw("same");
  h6->GetXaxis()->SetTitle("Inner Tracker Charge");
  h6->GetYaxis()->SetTitle("Normalized Entries");
  TLegend *lr = new TLegend(0.6,0.6,0.9,0.9);
  lr->AddEntry(h5,"Boron");
  lr->AddEntry(h6,"Carbon");
  lr->Draw("same");

  TLine *ll1 = new TLine(4.6,0,4.6,0.1);
  TLine *ll2 = new TLine(5.4,0,5.4,0.1);
  TLine *ll3 = new TLine(5.6,0,5.6,0.1);
  TLine *ll4 = new TLine(6.4,0,6.4,0.1);
  ll1->SetLineWidth(2);
  ll2->SetLineWidth(2);
  ll3->SetLineWidth(2);
  ll4->SetLineWidth(2);
  ll1->SetLineColor(4);
  ll2->SetLineColor(4);
  ll3->SetLineColor(3);
  ll4->SetLineColor(3);
  ll1->Draw("same");
  ll2->Draw("same");
  ll3->Draw("same");
  ll4->Draw("same");

  TLine *ll1 = new TLine(4.6,0,4.6,0.1);
  ll1->SetLineWidth(2);
  ll1->SetLineColor(4);
  ll1->Draw("same");

  line0->Draw("same");
}
