{
  //Settings
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  //Get histograms
  TFile f("study_cuts.root");
  TH1D *hsbeta[2], *hsrig[2];
  TH1D *hsg[10], *hsb[10];
  char name[100];
  for(int i=0; i<2; i++)
  {
    sprintf(name,"hsbeta_%d",i);
    hsbeta[i] = (TH1D*)f.Get(name);
    hsbeta[i]->SetLineWidth(2);
    if(i==1) hsbeta[i]->SetLineColor(2);
    else hsbeta[i]->SetLineStyle(2);

    sprintf(name,"hsrig_%d",i);
    hsrig[i] = (TH1D*)f.Get(name);
    hsrig[i]->SetLineWidth(2);
    if(i==1) hsrig[i]->SetLineColor(2);
    else hsrig[i]->SetLineStyle(2);
  }
  for(int i=0; i<10; i++)
  {
    sprintf(name,"hsg_%d",i);
    hsg[i] = (TH1D*)f.Get(name);
    hsg[i]->Scale(1/hsg[i]->Integral());
    hsg[i]->SetLineWidth(2);
    hsg[i]->SetLineColor(2);
    sprintf(name,"hsb_%d",i);
    hsb[i] = (TH1D*)f.Get(name);
    hsb[i]->Scale(1/hsb[i]->Integral());
    hsb[i]->SetLineWidth(2);
    hsb[i]->SetLineStyle(2);
  }

  ////////////////////
  //Check
  ////////////////////
  
  double scaleb = hsbeta[0]->Integral();
  hsbeta[0]->Scale(1.0/scaleb);
  hsbeta[1]->Scale(1.0/scaleb);
  TCanvas *cc0 = new TCanvas("cc0","cc0",1000,600);
  cc0->cd();
  cc0->SetLogy();
  cc0->SetGridx();
  cc0->SetGridy();
  hsbeta[0]->Draw();
  hsbeta[1]->Draw("same");
  hsbeta[0]->GetXaxis()->SetTitle("Beta");
  hsbeta[0]->GetYaxis()->SetTitle("Normalized Entries");
//  hsbeta[0]->GetXaxis()->SetTitleSize(0.047);
//  hsbeta[0]->GetYaxis()->SetTitleSize(0.047);
  TLegend *lc0 = new TLegend(0.5,0.6,0.9,0.9);
  lc0->AddEntry(hsbeta[0],"Before TOF Cuts");
  lc0->AddEntry(hsbeta[1],"After TOF Cuts");
  lc0->Draw("same");
  TLine *line1 = new TLine(0,0,1.5,0);
  line1->SetLineWidth(2);
  line1->Draw("same");
  
  double scaler = hsrig[0]->Integral();
  hsrig[0]->Scale(1.0/scaler);
  hsrig[1]->Scale(1.0/scaler);
  TCanvas *cc1 = new TCanvas("cc1","cc1",1000,600);
  cc1->cd();
  cc1->SetLogy();
  cc1->SetGridx();
  cc1->SetGridy();
  hsrig[0]->Draw();
  hsrig[1]->Draw("same");
  hsrig[0]->GetXaxis()->SetTitle("1/Rigidity (1/GV)");
  hsrig[0]->GetYaxis()->SetTitle("Normalized Entries");
//  hsrig[0]->GetXaxis()->SetTitleSize(0.047);
//  hsrig[0]->GetYaxis()->SetTitleSize(0.047);
  TLegend *lc1 = new TLegend(0.5,0.6,0.9,0.9);
  lc1->AddEntry(hsrig[0],"Before Tracker Cuts");
  lc1->AddEntry(hsrig[1],"After Tracker Cuts");
  lc1->Draw("same");
  TLine *line2 = new TLine(-1,0,1,0);
  line2->SetLineWidth(2);
  line2->Draw("same");

  ////////////////////
  //Draw TOF
  ////////////////////
 
  TCanvas *c0 = new TCanvas("c0","c0",1000,600);
  c0->cd();
  c0->SetLogy();
  c0->SetGridx();
  c0->SetGridy();
  hsg[0]->Draw();
  hsb[0]->Draw("same");
  hsg[0]->GetXaxis()->SetTitle("TOF Time Chi2/dof");
  hsg[0]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[0]->GetXaxis()->SetTitleSize(0.047);
  hsg[0]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l0 = new TLegend(0.5,0.6,0.9,0.9);
  l0->AddEntry(hsg[0],"Good TOF Sample");
  l0->AddEntry(hsb[0],"Bad TOF Sample");
  l0->Draw("same");

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  hsg[1]->Draw();
  hsb[1]->Draw("same");
  hsg[1]->GetXaxis()->SetTitle("TOF Coordinate Chi2/dof");
  hsg[1]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[1]->GetXaxis()->SetTitleSize(0.047);
  hsg[1]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l1 = new TLegend(0.5,0.6,0.9,0.9);
  l1->AddEntry(hsg[1],"Good TOF Sample");
  l1->AddEntry(hsb[1],"Bad TOF Sample");
  l1->Draw("same");
 
  TCanvas *c2 = new TCanvas("c2","c2",1000,600);
  c2->cd();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  hsg[2]->Draw();
  hsb[2]->Draw("same");
  hsg[2]->GetXaxis()->SetTitle("TOF Integer Charge Probability");
  hsg[2]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[2]->GetXaxis()->SetTitleSize(0.047);
  hsg[2]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l2 = new TLegend(0.5,0.6,0.9,0.9);
  l2->AddEntry(hsg[2],"Good TOF Sample");
  l2->AddEntry(hsb[2],"Bad TOF Sample");
  l2->Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3",1000,600);
  c3->cd();
  c3->SetLogy();
  c3->SetGridx();
  c3->SetGridy();
  hsg[3]->Draw();
  hsb[3]->Draw("same");
  hsg[3]->GetXaxis()->SetTitle("TOF Cleanliness");
  hsg[3]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[3]->GetXaxis()->SetTitleSize(0.047);
  hsg[3]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l3 = new TLegend(0.5,0.6,0.9,0.9);
  l3->AddEntry(hsg[3],"Good TOF Sample");
  l3->AddEntry(hsb[3],"Bad TOF Sample");
  l3->Draw("same");
 
  ////////////////////
  //Draw Tracker
  ////////////////////
  
  TCanvas *c4 = new TCanvas("c4","c4",1000,600);
  c4->cd();
  c4->SetLogy();
  c4->SetGridx();
  c4->SetGridy();
  hsg[4]->Draw();
  hsb[4]->Draw("same");
  hsg[4]->GetXaxis()->SetTitle("Tracker Rigidity Chi2/dof");
  hsg[4]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[4]->GetXaxis()->SetTitleSize(0.047);
  hsg[4]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l4 = new TLegend(0.5,0.6,0.9,0.9);
  l4->AddEntry(hsg[4],"Good Tracker Sample");
  l4->AddEntry(hsb[4],"Bad Tracker Sample");
  l4->Draw("same");

  TCanvas *c5 = new TCanvas("c5","c5",1000,600);
  c5->cd();
  c5->SetLogy();
  c5->SetGridx();
  c5->SetGridy();
  hsg[5]->Draw();
  hsb[5]->Draw("same");
  hsg[5]->GetXaxis()->SetTitle("Beta Rigidity Agreement");
  hsg[5]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[5]->GetXaxis()->SetTitleSize(0.047);
  hsg[5]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l5 = new TLegend(0.5,0.6,0.9,0.9);
  l5->AddEntry(hsg[5],"Good Tracker Sample");
  l5->AddEntry(hsb[5],"Bad Tracker Sample");
  l5->Draw("same");

  TCanvas *c6 = new TCanvas("c6","c6",1000,600);
  c6->cd();
  c6->SetLogy();
  c6->SetGridx();
  c6->SetGridy();
  hsg[6]->Draw();
  hsb[6]->Draw("same");
  hsg[6]->GetXaxis()->SetTitle("Inner Tracker Cleanliness");
  hsg[6]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[6]->GetXaxis()->SetTitleSize(0.047);
  hsg[6]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l6 = new TLegend(0.5,0.6,0.9,0.9);
  l6->AddEntry(hsg[6],"Good Tracker Sample");
  l6->AddEntry(hsb[6],"Bad Tracker Sample");
  l6->Draw("same");

  TCanvas *c7 = new TCanvas("c7","c7",1000,600);
  c7->cd();
  c7->SetLogy();
  c7->SetGridx();
  c7->SetGridy();
  hsg[7]->Draw();
  hsb[7]->Draw("same");
  hsg[7]->GetXaxis()->SetTitle("Inner Tracker Energy Deposition Discrepancy");
  hsg[7]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[7]->GetXaxis()->SetTitleSize(0.047);
  hsg[7]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l7 = new TLegend(0.5,0.6,0.9,0.9);
  l7->AddEntry(hsg[7],"Good Tracker Sample");
  l7->AddEntry(hsb[7],"Bad Tracker Sample");
  l7->Draw("same");

  TCanvas *c8 = new TCanvas("c8","c8",1000,600);
  c8->cd();
  c8->SetGridx();
  c8->SetGridy();
  hsg[8]->Draw();
  hsb[8]->Draw("same");
  hsg[8]->GetXaxis()->SetTitle("Inner Tracker Track Pattern");
  hsg[8]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[8]->GetXaxis()->SetTitleSize(0.047);
  hsg[8]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l8 = new TLegend(0.5,0.6,0.9,0.9);
  l8->AddEntry(hsg[8],"Good Tracker Sample");
  l8->AddEntry(hsb[8],"Bad Tracker Sample");
  l8->Draw("same");

  ////////////////////
  //Draw Fragmentation
  ////////////////////

  TCanvas *c9 = new TCanvas("c9","c9",1000,600);
  c9->cd();
  c9->SetGridx();
  c9->SetGridy();
  hsg[9]->Draw();
  hsb[9]->Draw("same");
  hsg[9]->GetXaxis()->SetTitle("Number of Energetic Secondary Tracks");
  hsg[9]->GetYaxis()->SetTitle("Normalized Entries");
  hsg[9]->GetXaxis()->SetTitleSize(0.047);
  hsg[9]->GetYaxis()->SetTitleSize(0.047);
  TLegend *l9 = new TLegend(0.5,0.6,0.9,0.9);
  l9->AddEntry(hsg[9],"Non-Fragmentation Sample");
  l9->AddEntry(hsb[9],"Fragmentation Sample");
  l9->Draw("same");
}
