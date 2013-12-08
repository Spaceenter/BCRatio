{
  TFile *f = new TFile("TempResult/RigErrBin.root","READ");

  TF1 *func = (TF1*)f->Get("frerr");

  int NBin;
  double CurrEdgeR, CurrEdgeE;
  double MC12 = 11.1779;
  double StartPoint = 1.2;
  double sigma;

  //Rigidity corrections
  cout<<"Rigidity corrections:"<<endl;
  NBin = 0;
  for(double i=0.2; i<2.8; i=i+0.2)
  {
    CurrEdgeR = TMath::Power(10,i);
    cout<<Form("%4.3f",CurrEdgeR)<<", ";
    NBin++;
  }
  cout<<endl<<"NBin = "<<NBin<<endl<<endl;


  //Rigidity 
  cout<<"Rigidity:"<<endl;
  NBin = 0;
  CurrEdgeR = StartPoint;
  for(;;)
  {
    cout<<Form("%4.3f",CurrEdgeR)<<", ";
    if(CurrEdgeR>TMath::Power(10,3.3)) break;
    NBin++;

    //Update
    double coeff = func->Eval(TMath::Log10(CurrEdgeR));
    if(coeff<0.1) coeff=0.1;
    sigma=3;
    CurrEdgeR += sigma*CurrEdgeR*coeff;
  }
  cout<<endl<<"NBin = "<<NBin<<endl<<endl;

  //Ek/A
  cout<<"Ek/A based on C12:"<<endl;
  NBin = 0;
  CurrEdgeR = StartPoint;
  for(;;)
  {
    CurrEdgeE = MC12/12.0*(TMath::Sqrt(CurrEdgeR*CurrEdgeR*36/MC12/MC12+1)-1);
    cout<<Form("%4.3f",CurrEdgeE)<<", ";
    if(CurrEdgeR>TMath::Power(10,3.3)) break;
    NBin++;

    //Update
    double coeff = func->Eval(TMath::Log10(CurrEdgeR));
    if(coeff<0.1) coeff=0.1;
    if(CurrEdgeR<TMath::Power(10,2.4)) sigma=2;
    else sigma=3;
    CurrEdgeR += sigma*CurrEdgeR*coeff;
  }
  cout<<endl<<"NBin = "<<NBin<<endl<<endl;
}
