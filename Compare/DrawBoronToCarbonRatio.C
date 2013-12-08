void DrawBoronToCarbonRatio() {

  gROOT->ProcessLine(".x CosmicNucleiData_HEAO.C");
  gROOT->ProcessLine(".x CosmicNucleiData_TRACER.C");
  gROOT->ProcessLine(".x CosmicNucleiData_CREAMII.C");
  gROOT->ProcessLine(".x CosmicNucleiData_AMS01.C");
  gROOT->ProcessLine(".x CosmicNucleiData_ATIC2.C");
  gROOT->ProcessLine(".x CosmicNucleiData_AMS02.C");
  gROOT->ProcessLine(".x l9CosmicNucleiData_AMS02.C");
//  gROOT->ProcessLine(".x try.C");

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->SetLogy();
  canvas->SetLogx();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetLogy(kFALSE);

  TH1F* frame = canvas->DrawFrame(1e-1,0.,10.e3,0.40);
  frame->SetXTitle("Kinetic Energy (GeV/n)");
  frame->SetYTitle("Boron-to-Carbon Ratio");
  frame->GetYaxis()->SetTitleOffset(1.20);
  HEAO_BC->SetMarkerSize(1.5);
  HEAO_BC->SetMarkerColor(kMagenta+1);
  HEAO_BC->SetMarkerStyle(24);
  HEAO_BC->Draw("P SAME");
  CREAMII_BC->SetMarkerSize(1.5);
  CREAMII_BC->SetMarkerColor(kAzure+1);
  CREAMII_BC->SetMarkerStyle(25);
  CREAMII_BC->Draw("P SAME");
  AMS01_BC->SetMarkerSize(1.5);
  AMS01_BC->SetMarkerColor(kGreen+1);
  AMS01_BC->SetMarkerStyle(26);
  AMS01_BC->Draw("P SAME");
  ATIC2_BC->SetMarkerSize(1.5);
  ATIC2_BC->SetMarkerColor(kGray+1);
  ATIC2_BC->SetMarkerStyle(27);
  ATIC2_BC->Draw("P SAME");
  TRACER_BC->SetMarkerSize(1.5);
  TRACER_BC->SetMarkerColor(kRed+1);
  TRACER_BC->SetMarkerStyle(28);
  TRACER_BC->Draw("P SAME");
  AMS02_BC->SetMarkerSize(1.5);
  AMS02_BC->SetMarkerColor(kBlue+1);
  AMS02_BC->SetMarkerStyle(21);
  AMS02_BC->Draw("P SAME");
  AMS02_BC_l9->SetMarkerSize(1.5);
  AMS02_BC_l9->SetMarkerColor(kRed+1);
  AMS02_BC_l9->SetMarkerStyle(21);
  AMS02_BC_l9->Draw("P SAME");
  /*
  AMS02_BC_try->SetMarkerSize(1.5);
  AMS02_BC_try->SetMarkerColor(kOrange+1);
  AMS02_BC_try->SetMarkerStyle(22);
  AMS02_BC_try->Draw("P SAME");
  */
  TLegend* legend = new TLegend(0.55,0.55,0.88,0.88,"","brNDC");
  legend->SetTextFont(52);
  legend->SetBorderSize(0);
  legend->AddEntry(HEAO_BC,"HEAO (A&A 1990)","P");
  legend->AddEntry(ATIC2_BC,"ATIC-2 (ICRC 2007)","P");
  legend->AddEntry(CREAMII_BC,"CREAM (Astropart. Phys. 2008)","P");
  legend->AddEntry(AMS01_BC,"AMS01 (ApJ 2009)","P");
  legend->AddEntry(TRACER_BC,"TRACER (ApJ 2011)","P");
  legend->AddEntry(AMS02_BC,"AMS02 L1+Inner","P");
  legend->AddEntry(AMS02_BC_l9,"AMS02 L1+Inner+L9","P");
 // legend->AddEntry(AMS02_BC_try,"AMS02 L1+Inner+L9 try","P");
  legend->Draw("SAME");
  canvas->Update();

}
