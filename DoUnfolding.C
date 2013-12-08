{
  //Set environment variables
  char libname[255];
  sprintf(libname,"$AMSWD/lib/%s/ntuple_slc4_PG.so",gSystem->GetBuildArch());
  gSystem->Load(libname);
  gInterpreter->AddIncludePath("$AMSSRC/include");
  TString curIncludePath(gSystem->GetIncludePath());
  gSystem->SetIncludePath("-I/$AMSSRC/include"+curIncludePath);

  BayesianUnfolder unfolder;

  TFile fm("Draw/matrix_1.05.root");
//  TFile fm("Draw/matrix_1.13.root");
//  TFile fm("Draw/matrix_0.97.root");
  TH2D *hmatrix = (TH2D*)fm.Get("hmatrix");

  TFile fs("Draw/result.root");
  TH1D *hunb = (TH1D*)fs.Get("hunb");
  TH1D *hunc = (TH1D*)fs.Get("hunc");

  TH1D hreb = *hunb;
  hreb.Reset();
  TH1D hrec = *hunc;
  hrec.Reset();

  //Temp modifications
  for(int i=0; i<hunb->GetNbinsX(); i++)
  {
    if(hunb->GetBinCenter(i+1)<-0.0015) hunb->SetBinContent(i+1,0);
  }
  for(int i=0; i<hunc->GetNbinsX(); i++)
  {
    if(hunc->GetBinCenter(i+1)<-0.0015) hunc->SetBinContent(i+1,0);
  }

  //Unfolding
  unfolder.computeAll(*hmatrix, *hunb, hreb, 100, true, true, 10);
  unfolder.computeAll(*hmatrix, *hunc, hrec, 100, true, true, 10);

  //Change the format
  TH1D *hunreb= (TH1D*)hreb.Clone("hunreb");
  TH1D *hunrec= (TH1D*)hrec.Clone("hunrec");

  //Calculate the ratio
  TH1D *hunrebc = (TH1D*)hunreb->Clone("hunrebc");
  hunrebc->Divide(hunrec);
  TH1D *hunbc = (TH1D*)hunb->Clone("hunbc");
  hunbc->Divide(hunc);
  for(int i=1; i<=hunbc->GetNbinsX(); i++)
  {
    if(hunb->GetBinContent(i)==0) continue;
    double a = hunb->GetBinContent(i);
    double ea = TMath::Sqrt(a);
    double b = hunc->GetBinContent(i);
    double eb = TMath::Sqrt(b);
    hunbc->SetBinError(i,a/b*sqrt(ea*ea/a/a+eb*eb/b/b));
  }

  //Store results
  TFile *fr = new TFile("unfolding.root","RECREATE");
  fr->cd();
  hunb->Write();
  hunc->Write();
  hunbc->Write();
  hunreb->Write();
  hunrec->Write();
  hunrebc->Write();
  fr->Write();
  fr->Close();

  //Print out the correction factors
  cout<<"Correction factors of unfolding:"<<endl;
  TH1D *hcorr = (TH1D*)hunrebc->Clone("hcorr");
  hcorr->Divide(hunbc);
  int count;
  cout<<endl<<"Bin edges:"<<endl;
  count = 0;
  for(int i=hcorr->GetNbinsX(); i>0; i--)
  {
    double BinRightEdge = hcorr->GetBinLowEdge(i) + hcorr->GetBinWidth(i);
    if(1/BinRightEdge > 4000) break;
    cout<<1/BinRightEdge<<", ";
    count++;
  }
  cout<<endl<<"count = "<<count<<endl;
  cout<<endl<<"Correction factors:"<<endl;
  count = 0;
  for(int i=hcorr->GetNbinsX(); i>0; i--)
  {
    double BinLeftEdge = hcorr->GetBinLowEdge(i);
    if(1/BinLeftEdge > 4000) break;
    cout<<hcorr->GetBinContent(i)<<", ";
    count++;
  }
  cout<<endl<<"count = "<<count<<endl;
  cout<<endl;
}
