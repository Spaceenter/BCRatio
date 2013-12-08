{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  TFile fm("toymc_material.root");
  TCanvas *ct = new TCanvas("ct","ct",1000,600);

  ct->cd();
  ntuple->Draw("grammage:z_limit");
  TGraph *grm = new TGraph(ntuple->GetSelectedRows(),ntuple->GetV2(),ntuple->GetV1());

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetGridx();
  c1->SetGridy();
  TH1F *frame1 = c1->DrawFrame(-95,0,175,20);
  frame1->GetXaxis()->SetTitle("Z (cm)");
  frame1->GetYaxis()->SetTitle("Amount of Material Above Z (g/cm^{2})");
  grm->SetMarkerStyle(21);
  grm->SetMarkerSize(0.5);
  grm->SetMarkerColor(kRed);
  grm->Draw("P same");

  ct->cd();
  TFile fs("toymc_sur_prob.root");
  ntuple->Draw("survival[6]:E/1000","pid==0 && E/1000>0.007");
  TGraph *grs0 = new TGraph(ntuple->GetSelectedRows(),ntuple->GetV2(),ntuple->GetV1());
  ntuple->Draw("survival[6]:E/1000","pid==1 && E/1000>0.007");
  TGraph *grs1 = new TGraph(ntuple->GetSelectedRows(),ntuple->GetV2(),ntuple->GetV1());
  ntuple->Draw("survival[6]:E/1000","pid==2 && E/1000>0.007");
  TGraph *grs2 = new TGraph(ntuple->GetSelectedRows(),ntuple->GetV2(),ntuple->GetV1());

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetLogx();
  c2->SetGridx();
  c2->SetGridy();
  TH1F *frame2 = c2->DrawFrame(0.007,0.33,100,0.55);
  frame2->GetXaxis()->SetTitle("E_{k}/A (GeV/n)");
  frame2->GetYaxis()->SetTitle("Survival Probability at Z = -75cm");
  grs0->SetMarkerStyle(23);
  grs0->SetMarkerSize();
  grs0->SetMarkerColor(kRed);
  grs0->Draw("P same");
  grs0->GetXaxis()->SetTitle("E_{k}/A (GeV/n)");
  grs0->GetYaxis()->SetTitle("Survival Probability at Z = -75cm");
  grs0->GetYaxis()->SetRangeUser(0.33,0.55);
  grs1->SetMarkerStyle(22);
  grs1->SetMarkerSize(1);
  grs1->SetMarkerColor(kBlue);
  grs1->Draw("P same");
  grs2->SetMarkerStyle(20);
  grs2->SetMarkerSize(1);
  grs2->SetMarkerColor(kGreen);
  grs2->Draw("P same");
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(grs0,"^{10}B","P");
  leg->AddEntry(grs1,"^{11}B","P");
  leg->AddEntry(grs2,"^{12}C","P");
  leg->Draw("same");
}
