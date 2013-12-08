{
  TFile f("../Draw/result.root");
  TH1D *hunb = (TH1D*)f.Get("hunb");

  int NBins = hunb->GetNbinsX();
  int count = 0;
  for(int i=1; i<=NBins; i++)
  {
    cout<<hunb->GetBinLowEdge(i)<<", ";
    count++;
  }
  cout<<hunb->GetBinLowEdge(NBins)+hunb->GetBinWidth(NBins)<<endl;
  cout<<"count = "<<count+1<<endl;
}
