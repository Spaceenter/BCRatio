{
  TFile f("result.root");

  TH1D *hbco = (TH1D*)f.Get("hbco");
  hbco->SetLineColor(2);
  TH1D *hbcf = (TH1D*)f.Get("hbcf");
  hbcf->SetLineStyle(2);
  TH1D *hbct = (TH1D*)f.Get("hbct");
  hbct->SetLineColor(4);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  c1->SetLogx();
  c1->SetLogy();
  hbco->Draw();
  hbcf->Draw("same");
  hbco->GetXaxis()->SetRangeUser(0.001,0.7);
  hbco->GetYaxis()->SetRangeUser(0.07,0.5);

  TCanvas *c2 = new TCanvas();
  c2->cd();
  c2->SetLogx();
  c2->SetLogy();
  hbct->Draw();
  hbcf->Draw("same");
  hbct->GetXaxis()->SetRangeUser(0.001,0.7);
  hbct->GetYaxis()->SetRangeUser(0.07,0.5);

  TCanvas *c3 = new TCanvas();
  c3->cd();
  c3->SetLogx();
  c3->SetLogy();
  hbct->Draw();
  hbco->Draw("same");
  hbct->GetXaxis()->SetRangeUser(0.001,0.7);
  hbct->GetYaxis()->SetRangeUser(0.07,0.5);

  TH1D *hbt = (TH1D*)f.Get("hbt");
  hbt->SetLineColor(2);
  TH1D *hct = (TH1D*)f.Get("hct");
  hct->SetLineColor(2);
  TH1D *hbf = (TH1D*)f.Get("hbf");
  hbf->SetLineColor(4);
  TH1D *hcf = (TH1D*)f.Get("hcf");
  hcf->SetLineColor(4);

  TFile fo("result.root");
  double scale = 50;
  TH1D *hbo = (TH1D*)fo.Get("hunb");
  hbo->Scale(scale);
  TH1D *hco = (TH1D*)fo.Get("hunc");
  hco->Scale(scale);

  TCanvas *c4 = new TCanvas();
  c4->cd();
  c4->SetLogx();
  c4->SetLogy();
  hbt->Draw();
  hbf->Draw("same");
  hct->Draw("same");
  hcf->Draw("same");
  hbo->Draw("same");
  hco->Draw("same");
  hbt->GetXaxis()->SetRangeUser(0.001,0.7);
  hbt->GetYaxis()->SetRangeUser(3000,25000000);

  TH1D *hsigma = (TH1D*)f.Get("hsigma");

  TCanvas *c5 = new TCanvas();
  c5->cd();
  c5->SetLogx();
  hsigma->SetMarkerStyle(22);
  hsigma->Draw("P");
  hsigma->GetXaxis()->SetRangeUser(0.001,0.1);

  //print bin centers
  for(int i=1; i<=hsigma->GetNbinsX(); i++)
  {
    double center = 1.0/(hsigma->GetBinLowEdge(i)+0.5*hsigma->GetBinWidth(i));
    if(center<0 || center>1500) continue;
    cout<<Form("%4.3f",center)<<endl;
  }

  //upper limit
  int NBin = 42;
  double XBin[43] = {
    -0.794281, -0.582072, -0.427533, -0.322789, -0.245942, -0.187758, -0.145243, -0.11217, -0.0866927, -0.0668673, -0.050984, -0.0386518, -0.0288675, -0.0210975, -0.0150987, -0.0105141, -0.00718257, -0.00475627, -0.00299555, -0.00173286, -0.000745535, 0, 0.000745535, 0.00173286, 0.00299555, 0.00475627, 0.00718257, 0.0105141, 0.0150987, 0.0210975, 0.0288675, 0.0386518, 0.050984, 0.0668673, 0.0866927, 0.11217, 0.145243, 0.187758, 0.245942, 0.322789, 0.427533, 0.582072, 0.794281
  };
  double BinCenter[9] = {22.313, 29.621, 40.028, 55.254, 78.086, 113.016, 167.520, 258.004, 422.975};
  double Value[9] = {0, 0, 0, 0.000806485, 0.028071, 0.0737699, 0.133113, 0.228461, 0.470343};
  double Error[9] = {0.00232418, 0.00285193, 0.00542002, 0.0613907, 0.0224958, 0.0348759, 0.0579517, 0.109124, 0.245548};
  TH1D *hlimit = new TH1D("hlimit","hlimit",NBin,XBin);
  for(int i=0; i<9; i++)
  {
    int iBin = hlimit->FindBin(1.0/BinCenter[i]);
    hlimit->SetBinContent(iBin,Value[i]);
    hlimit->SetBinError(iBin,2*Error[i]);
  }
  TCanvas *c6 = new TCanvas();
  c6->cd();
  c6->SetLogx();
  hlimit->Draw("P");
  hlimit->SetLineColor(2);
  hlimit->SetLineWidth(2);
  hlimit->GetXaxis()->SetRangeUser(0.001,0.1);
}
