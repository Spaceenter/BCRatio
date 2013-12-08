{
  //Environment
  gSystem->Load("libRooFit");
  using namespace RooFit;

  //Read histogram
  TFile ff("Draw/int_mc_dpmjet.root");
  TH1D *htot[4], *hrec[4][18], *htmp[4][18];
  char name[100];
  for(int i=0;i<4;i++)
  {
    sprintf(name,"htot_%d",i);
    htot[i] = (TH1D*)ff.Get(name);
    for(int j=0;j<18;j++)
    {
      sprintf(name,"hrec_%d_%d",i,j+2);
      hrec[i][j] = (TH1D*)ff.Get(name);
      sprintf(name,"htmp_%d_%d",i,j+2);
      htmp[i][j] = (TH1D*)ff.Get(name);
    }
  }

  //RooFit
  RooRealVar *Charge = new RooRealVar("Charge","Charge",0,10);
  RooRealVar *Frac[4][18];
  RooDataHist *Spec[4][18], *Temp[4][18];
  RooHistPdf *PdfH[4][18];
  RooAddPdf *Pdf[4][18];
  RooPlot *qframe[4][18];
  double rup, rlow;
  int xup, xlow;
  double NSur[4][18], Prob[4][18];
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<18;j++)
    {
      //Create objects
      sprintf(name,"Frac_%d_%d",i,j+2);
      Frac[i][j] = new RooRealVar(name,name,0,100000);
      sprintf(name,"Spec_%d_%d",i,j+2);
      Spec[i][j] = new RooDataHist(name,name,RooArgSet(*Charge),hrec[i][j]);
      sprintf(name,"Temp_%d_%d",i,j+2);
      Temp[i][j] = new RooDataHist(name,name,RooArgSet(*Charge),htmp[i][j]);
      sprintf(name,"PdfH_%d_%d",i,j+2);
      PdfH[i][j] = new RooHistPdf(name,name,RooArgSet(*Charge),*Temp[i][j],0);
      sprintf(name,"Pdf_%d_%d",i,j+2);
      Pdf[i][j] = new RooAddPdf(name,name,RooArgList(*PdfH[i][j]),RooArgList(*Frac[i][j]));
      sprintf(name,"qframe_%d_%d",i,j+2);
      qframe[i][j] = Charge->frame();
      qframe[i][j]->SetNameTitle(name,name);

      //Range
      if(i==0 || i==1)
      {
	rlow=4.4;
	rup=5.1;
	xlow=89;
	xup=102;
      }
      else if(i==2)
      {
	rlow=5.4;
	rup=6.1;
	xlow=109;
	xup=122;
      }
      else
      {
	rlow=7.4;
	rup=8.1;
	xlow=149;
	xup=162;
      }

      //Fit
      Pdf[i][j]->fitTo(*Spec[i][j],Extended(),Range(rlow,rup));

      //Plot 
      Spec[i][j]->plotOn(qframe[i][j]);
      Pdf[i][j]->plotOn(qframe[i][j],Range(0,10));

      /*
      //Survival probability
      Prob[i][j]=Frac[i][j]->getVal()/htmp[i][j]->Integral(xlow,xup)*htmp[i][j]->Integral();
      Prob[i][j]=Prob[i][j]/Spec[i][j]->sum(kFALSE);
      */

      //Survived events
      NSur[i][j]=Frac[i][j]->getVal()/htmp[i][j]->Integral(xlow,xup)*htmp[i][j]->Integral();
    }
  }

  //Save results into histograms
  TH1D *hnsur[4], *hsurp[4];
  double a, b, ea, eb;
  for(int i=0;i<4;i++)
  {
    sprintf(name,"hnsur_%d",i);
    hnsur[i] = new TH1D(name,name,18,0.2,2);
    sprintf(name,"hsurp_%d",i);
    hsurp[i] = new TH1D(name,name,18,0.2,2);
    for(int j=0;j<18;j++) 
    {
      hnsur[i]->SetBinContent(j+1,NSur[i][j]);
      a=hnsur[i]->GetBinContent(j+1);
      b=htot[i]->GetBinContent(j+1);
      ea=TMath::Sqrt(a);
      eb=TMath::Sqrt(b);
      hsurp[i]->SetBinContent(j+1,a/b);
      hsurp[i]->SetBinError(j+1,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));
    }
  }

  //Calculate survival probability ratios
  TH1D *hsurpr10 = (TH1D*)hsurp[0]->Clone("hsurpr10");
  TH1D *hsurpr11 = (TH1D*)hsurp[1]->Clone("hsurpr11");
  for(int r=1; r<=18; r++)
  {
    a=hsurp[0]->GetBinContent(r);
    b=hsurp[2]->GetBinContent(r);
    ea=hsurp[0]->GetBinError(r);
    eb=hsurp[2]->GetBinError(r);
    hsurpr10->SetBinContent(r,a/b);
    hsurpr10->SetBinError(r,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));

    a=hsurp[1]->GetBinContent(r);
    b=hsurp[2]->GetBinContent(r);
    ea=hsurp[1]->GetBinError(r);
    eb=hsurp[2]->GetBinError(r);
    hsurpr11->SetBinContent(r,a/b);
    hsurpr11->SetBinError(r,a/b*TMath::Sqrt(ea*ea/a/a+eb*eb/b/b));
  }

  //Save results into root file
  TFile *f = new TFile("mc_sur_prob.root","RECREATE");
  f->cd();
  hsurpr10->Write();
  hsurpr11->Write();
  for(int r=0; r<18; r++) qframe[2][r]->Write();
  f->Write();
  f->Close();
}
