{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  TFile f("matrix_1.05.root");

  TF1 *ff = new TF1("ff","[0]*(1/[1]*TMath::Exp(-(x)*(x)/2/[1]/[1]) + [3]/[2]*TMath::Exp(-(x)*(x)/2/[2]/[2]))",-5,5);
  ff->SetParameter(0,0.01);
  ff->SetParameter(1,0.3);
  ff->SetParameter(2,0.7);
  ff->SetParameter(3,0.3);
  ff->SetLineWidth(2);
  ff->SetLineColor(2);

  TH1D *ha = (TH1D*)f.Get("hres_256");
  ha->SetLineWidth(2);
  ha->Fit("ff","","",-2.5,2.5);

  TCanvas *c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  ha->Draw();
  ha->GetXaxis()->SetTitle("(1/R_{m}-1/R_{t})/(1/R_{t})");
  ha->GetYaxis()->SetTitle("Normalized Entries");
  ff->Draw("same");
  TLine *ll = new TLine(-5,0,5,0);
  ll->SetLineWidth(2);
  ll->Draw("same");

  TH2D *hh = (TH2D*)f.Get("hmatrix");

  TCanvas *cc = new TCanvas("cc","cc",1000,600);
  cc->cd();
  cc->SetGridx();
  cc->SetGridy();
  hh->Draw("colz");
  hh->GetXaxis()->SetTitle("1/R_{t} (1/GV)");
  hh->GetYaxis()->SetTitle("1/R_{m} (1/GV)");
  hh->GetZaxis()->SetRangeUser(0.005,0.52);
}
