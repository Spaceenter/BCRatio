{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("ExpT.root");

  TH1D *h11 = (TH1D*)f.Get("etb11");
  h11->Scale(1.0/h11->GetBinContent(25));
  h11->SetLineWidth(2);
  h11->SetLineColor(2);

  TH1D *h12 = (TH1D*)f.Get("etc12");
  h12->Scale(1.0/h12->GetBinContent(25));
  h12->SetLineWidth(2);
  h12->SetLineStyle(2);

  TCanvas *c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetGridx();
  c->SetGridy();
  h11->Draw();
  h12->Draw("same");
  h11->GetXaxis()->SetTitle("Log_{10}(E_{k}/A)");
  h11->GetYaxis()->SetTitle("Normalized Exposure Time");
  TLegend *l = new TLegend(0.6,0.6,0.9,0.9);
  l->AddEntry(h11,"^{11}B");
  l->AddEntry(h12,"^{10}B ^{12}C");
  l->Draw("same");
}
