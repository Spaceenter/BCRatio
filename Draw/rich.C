{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f6("rich_60.root");
  TFile f7("rich_70.root");

  TH1D *h6 = (TH1D*)f6.Get("hebcrig");
  h6->SetLineWidth(3);
  h6->SetLineColor(2);
  h6->SetLineStyle(2);
  for(int i=0; i<h6->GetNbinsX(); i++) //Add SurProb
  {
    h6->SetBinContent(i+1,h6->GetBinContent(i+1)/1.055);
  }

  TH1D *h7 = (TH1D*)f7.Get("hebcrig");
  h7->SetLineWidth(3);
  h7->SetLineColor(4);
  h7->SetLineStyle(3);
  for(int i=0; i<h7->GetNbinsX(); i++) //Add SurProb
  {
    h7->SetBinContent(i+1,h7->GetBinContent(i+1)/1.055);
  }

  TH1D *hr = (TH1D*)f7.Get("hebcbetar");
  hr->SetLineWidth(3);
  for(int i=0; i<hr->GetNbinsX(); i++) //Add SurProb
  {
    hr->SetBinContent(i+1,hr->GetBinContent(i+1)/1.055);
  }

  TCanvas *c = new TCanvas("c","c",1000,600);
  c->cd();
  c->SetLogx();
  c->SetGridx();
  c->SetGridy();
  h6->Draw();
  h7->Draw("same");
  hr->Draw("same");
  h6->GetXaxis()->SetTitle("E_{k}/A (GeV/n)");
  h6->GetYaxis()->SetTitle("Boron to Carbon Ratio");
  TLegend *l = new TLegend(0.6,0.6,0.9,0.9);
  l->AddEntry(h6,"^{11}B/(^{11}B+^{10}B)=0.6");
  l->AddEntry(h7,"^{11}B/(^{11}B+^{10}B)=0.7");
  l->AddEntry(hr,"RICH");
  l->Draw("same");
}
