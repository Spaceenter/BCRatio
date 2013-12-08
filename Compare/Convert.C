{
  TFile f("rr.root");
  double SurProb = 0.055;
  double RelSysErr = 0.03;

  //**********************************
  //L1 + Inner
  //**********************************

  cout<<"L1 + Inner :"<<endl;

  TH1D *hebc = (TH1D*)f.Get("hebcrig");
  cout<<"NBin = "<<hebc->GetNbinsX()<<endl;

  for(int i=1; i<=hebc->GetNbinsX(); i++)
  {
    double x = hebc->GetBinCenter(i);
    cout<<x<<", ";
  }
  cout<<endl;

  for(int i=0;i<hebc->GetNbinsX();i++)
  {
    cout<<"0, ";
  }
  cout<<endl;

  for(int i=1; i<=hebc->GetNbinsX(); i++)
  {
    cout<<hebc->GetBinContent(i)/(1+SurProb)<<", ";
  }
  cout<<endl;

  for(int i=1; i<=hebc->GetNbinsX(); i++)
  {
    double staterr = hebc->GetBinError(i);
    double syserr = hebc->GetBinContent(i)/(1+SurProb)*RelSysErr; 
    cout<<TMath::Sqrt(staterr*staterr+syserr*syserr)<<", ";
  }
  cout<<endl<<endl;

  //**********************************
  //L1 + Inner + L9
  //**********************************

  cout<<"L1 + Inner + L9:"<<endl;

  TH1D *hebcl9 = (TH1D*)f.Get("hebcrigl9");
  cout<<"NBin = "<<hebcl9->GetNbinsX()<<endl;

  for(int i=1; i<=hebcl9->GetNbinsX(); i++)
  {
    double x = hebcl9->GetBinCenter(i);
    cout<<x<<", ";
  }
  cout<<endl;

  for(int i=0;i<hebcl9->GetNbinsX();i++)
  {
    cout<<"0, ";
  }
  cout<<endl;

  for(int i=1; i<=hebcl9->GetNbinsX(); i++)
  {
    cout<<hebcl9->GetBinContent(i)/(1+SurProb)<<", ";
  }
  cout<<endl;

  for(int i=1; i<=hebcl9->GetNbinsX(); i++)
  {
    double staterr = hebcl9->GetBinError(i);
    double syserr = hebcl9->GetBinContent(i)/(1+SurProb)*RelSysErr;
    cout<<TMath::Sqrt(staterr*staterr+syserr*syserr)<<", ";
  }
  cout<<endl<<endl;

}
