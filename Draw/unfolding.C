{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile f("unfolding.root");

  TH1D *hunb = (TH1D*)f.Get("hunb");
  hunb->SetLineWidth(2);
  hunb->SetLineStyle(2);
  TH1D *hunreb = (TH1D*)f.Get("hunreb");
  hunreb->SetLineWidth(2);
  hunreb->SetLineColor(2);
  TH1D *hunc = (TH1D*)f.Get("hunc");
  hunc->SetLineWidth(2);
  hunc->SetLineStyle(2);
  TH1D *hunrec = (TH1D*)f.Get("hunrec");
  hunrec->SetLineWidth(2);
  hunrec->SetLineColor(2);

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  hunc->Draw("E");
  hunb->Draw("Esame");
  hunrec->Draw("same");
  hunreb->Draw("same");
  hunc->GetXaxis()->SetTitle("1/Rigidity (1/GV)");
  hunc->GetYaxis()->SetTitle("Raw Counts of Boron and Carbon");
  hunc->GetXaxis()->SetRangeUser(0.001,0.25);
  hunc->GetYaxis()->SetRangeUser(80,500000);
  TLegend *l1 = new TLegend(0.6,0.6,0.9,0.9);
  l1->AddEntry("hunc","Before Unfolding");
  l1->AddEntry("hunrec","After Unfolding");
  l1->Draw("same");

  TH1D *hunbc = (TH1D*)f.Get("hunbc");
  hunbc->SetLineWidth(2);
  hunbc->SetLineStyle(2);
  TH1D *hunrebc = (TH1D*)f.Get("hunrebc");
  hunrebc->SetLineWidth(2);
  hunrebc->SetLineColor(2);

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetLogx();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  hunbc->Draw();
  hunrebc->Draw("same");
  hunbc->GetXaxis()->SetTitle("1/Rigidity (1/GV)");
  hunbc->GetYaxis()->SetTitle("Raw Boron to Carbon Ratio");
  hunbc->GetXaxis()->SetRangeUser(0.001,0.25);
  hunbc->GetYaxis()->SetRangeUser(0.065,0.345);
  TLegend *l2 = new TLegend(0.6,0.6,0.9,0.9);
  l2->AddEntry("hunbc","Before Unfolding");
  l2->AddEntry("hunrebc","After Unfolding");
  l2->Draw("same");
}
