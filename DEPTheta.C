double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

void DEPTheta()
{
  gStyle->SetPalette(1);

  TChain *t = new TChain("treeUP","tree");
  t->Add("2417_2448.root");
  //t->Add("../../EnergyCalibration/analysis/V6/2543_noReclustering_Fiducial.root");
  //t->Add("../../EnergyCalibration/analysis/V6/2578_noReclustering_Fiducial.root");

  const int nrHistos = 20;

  TH1F *hSingle[nrHistos];
  TH1F *hMulti[nrHistos];
  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,3500,200,0,15000);

  TF1 *fitSingle[nrHistos];

  double *angle = new double[nrHistos];
  double *res = new double[nrHistos];
  double *resErr = new double[nrHistos];

  for (int i = 0; i < nrHistos; i++) {
     char hSingleName[100];
     char hSingleTitle[100];

     double theta = TMath::Pi()/9.0/nrHistos*i;

     sprintf(hSingleName,"hSingle_%i",i);
     sprintf(hSingleTitle,"%.0f rotation angle",180.0/TMath::Pi()*theta);

     hSingle[i] = new TH1F(hSingleName,hSingleTitle,200,0,6000);

     char SingleCMD[200];
     char MultiCMD[200];
     sprintf(SingleCMD,"(fcsc-2*2075.75)*TMath::Sin(%.2f) + fepclUP*TMath::Cos(%.2f) >> %s",theta,theta,hSingleName);

     t->Draw(SingleCMD);

     char fitSingleName[100];
     sprintf(fitSingleName,"fit_%i",i);

     double cscOff = 0;

     double fitRange1 = (4000-cscOff)*TMath::Sin(theta) + 1500*TMath::Cos(theta);
     double fitRange2 = (5600-cscOff)*TMath::Sin(theta) + 1750*TMath::Cos(theta);
     double Peak = (5000-cscOff)*TMath::Sin(theta) + 1625*TMath::Cos(theta);

     fitSingle[i] = new TF1(fitSingleName,"gaus",fitRange1,fitRange2);

     fitSingle[i]->SetParNames("A1","E1","#sigma");

     fitSingle[i]->SetParameters(10000,Peak,100);

     fitSingle[i]->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
     fitSingle[i]->SetParLimits(2,50,200);

     fitSingle[i]->SetLineWidth(1);
     fitSingle[i]->SetLineColor(kBlue);

     hSingle[i]->Fit(fitSingleName,"r");

     double par[4];
     fitSingle[i]->GetParameters(par);
     double par1 = fitSingle[i]->GetParameter(1);
     double par2 = fitSingle[i]->GetParameter(2);
     double *parErr = fitSingle[i]->GetParErrors();

     angle[i] = theta;
     res[i] = fitSingle[i]->GetParameter(2)/fitSingle[i]->GetParameter(1);
     resErr[i] = TMath::Sqrt(1.0 / par1**2 * parErr[2]**2 + (par[2] / par1**2)**2 * parErr[1]**2);
  }

  t->Draw("(fcsc-2*2075.75):fepclUP>>h2DSingle");

  TGraphErrors *grSingleRes = new TGraphErrors(nrHistos,angle,res,0,resErr);
  grSingleRes->SetMarkerStyle(20);
  grSingleRes->SetMarkerSize(0.8);

  TF1 *fitRes = new TF1("fitRes","pol2",0.05,0.15);
  fitRes->SetLineWidth(1);
  fitRes->SetLineColor(kBlue);

  grSingleRes->Fit("fitRes","r");

  double parRes1 = fitRes->GetParameter(2);
  double parRes2 = fitRes->GetParameter(1);
  double *parResErr = fitRes->GetParErrors();

  cout << "Best resolution: " << -0.5*parRes2/parRes1 << " +- " << TMath::Sqrt(1.0/(2*parRes1**2) * parResErr[1]**2 + (parRes2/(2*parRes1**2))**2 * parResErr[2]**2) << endl;

  TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");

  TCanvas *c3 = new TCanvas();
  c3->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
     c3->cd(i+1);
     hSingle[i]->Draw();
     cout << res[i] << endl;
  }

  TCanvas *c5 = new TCanvas();
  grSingleRes->Draw("AZP");

  return;
}
