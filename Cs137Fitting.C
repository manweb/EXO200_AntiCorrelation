double *PeakPos1;
double *PeakPosErr1;
double *PeakPosMulti;
double *PeakPosErrMulti;

void ProcessRun(int RunID, int ID);

TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();

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

double fitFunctionMulti(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));

  return gauss1 + erf1;
}

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

void Cs137Fitting()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
/*  const int nrRuns = 7;
  int runID[nrRuns] = {2109, 2110, 2111, 2112, 2113, 2117, 2119};
  double runIDCopy[nrRuns] = {2109, 2110, 2111, 2112, 2113, 2117, 2119};
*/
  /*const int nrRuns = 5;
  int runID[nrRuns] = {2410, 2449, 2450, 2469, 2473};
  double runIDCopy[nrRuns] = {2410, 2449, 2450, 2469, 2473};*/

  //int runID[nrRuns] = {2410, 2469};
  //double runIDCopy[nrRuns] = {2410, 2469};

  const int nrRuns = 3;
  int runID[nrRuns] = {2450, 2469, 2473};
  double runIDCopy[nrRuns] = {2450, 2469, 2473};

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPosMulti = new double[nrRuns];
  PeakPosErrMulti = new double[nrRuns];
  
  c1->Print("Cs137_ss.ps[");
  c2->Print("Cs137_ms.ps[");
  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }
  c1->Print("Cs137_ss.ps]");
  c2->Print("Cs137_ms.ps]");

  cout << "Mean (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (multi site) = " << TMath::Mean(nrRuns,PeakPosMulti) << " +- " << TMath::RMS(nrRuns,PeakPosMulti) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPosMulti,0,PeakPosErrMulti);

  gr1->SetTitle("Cs137 peak (single site)");
  gr2->SetTitle("Cs137 peak (multi site)");

  gr1->SetName("Cs137_single");
  gr2->SetName("Cs137_multi");

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetYaxis()->SetRangeUser(1050,1100);
  gr1->GetXaxis()->SetLimits(1560,2480);
  gr1->SetTitle("Th228 Peak");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetYaxis()->SetRangeUser(1050,1100);
  gr2->GetXaxis()->SetLimits(1560,2480);
  gr2->SetTitle("Th228 Peak (multi site)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");

  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");

  /*TFile *fOut = new TFile("../analysis/Cs137_PeakPosition_alpha_220212.root","RECREATE");
  gr1->Write();
  gr2->Write();

  fOut->Close();*/

  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  //sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",RunID);
  sprintf(FileName,"../../EnergyCalibration/analysis/V11/%i_noReclustering_Fiducial.root",RunID);
  //sprintf(FileName,"../analysis/V3/%i_9us.root",RunID);

  cout << "Processing run " << RunID << "..." << endl;

  TFile *f = new TFile(FileName,"READ");
/*  TH1F *hMean = (TH1F*)f->Get("hMean");
  TH1F *hN = (TH1F*)f->Get("hN");
  TH1F *hP = (TH1F*)f->Get("hP");
  TH1F *hMulti = (TH1F*)f->Get("hMulti");
*/
  TTree *t = (TTree*)f->Get("t");

  double ecrec;
  int nsite;

  t->SetBranchAddress("ecrec",&ecrec);
  t->SetBranchAddress("nsite",&nsite);

  char hName_ss[100];
  char hName_ms[100];
  sprintf(hName_ss,"Single Site Spectrum, %i",RunID);
  sprintf(hName_ms,"Multi Site Spectrum, %i",RunID);
  
  TH1F *hMean = new TH1F("hMean",hName_ss,100,0,2500);
  TH1F *hMulti = new TH1F("hMulti",hName_ms,100,0,2500);
  //TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,8000);
  //TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,8000);

  //t->Draw("ecrec_gc>>hMean","nsite == 1");
  //t->Draw("ecrec>>hMean","nsite == 1 && TMath::Abs(zc) < 150");
  //t->Draw("ecrec_gc>>hMulti","nsite > 1");
  //t->Draw("(csc*0.3148-119.5)*TMath::Sin(0.37082) + (ecrec*0.9647+67.4)*TMath::Cos(0.37082)>>hMean","nsite == 1");
  //t->Draw("csc*TMath::Sin(0.1956) + ecrec*TMath::Cos(0.1956)>>hMean","nsite == 1");
  //t->Draw("csc*TMath::Sin(0.2267) + ecrec*TMath::Cos(0.2267)>>hMulti","nsite > 1");
  
  // Reprocesse data (01/30/12)
  //t->Draw("csc*TMath::Sin(0.1916) + epcrec*TMath::Cos(0.1916)>>hMean","nsite == 1");
  //t->Draw("csc*TMath::Sin(0.2164) + epcrec*TMath::Cos(0.2164)>>hMulti","nsite > 1");
  
  // Reprocesse data (02/07/12)
  //t->Draw("csc*TMath::Sin(0.1867) + epcrec*TMath::Cos(0.1867)>>hMean","nsite == 1");
  //t->Draw("csc*TMath::Sin(0.2084) + epcrec*TMath::Cos(0.2084)>>hMulti","nsite > 1");
  
  // Reprocesse data (02/22/12)
  //t->Draw("csc*TMath::Sin(0.1860) + epcrec*TMath::Cos(0.1860)>>hMean","nsite == 1");
  //t->Draw("csc*TMath::Sin(0.2064) + epcrec*TMath::Cos(0.2064)>>hMulti","nsite > 1");
  
  // Reprocesse data (03/21/12)
  t->Draw("csc*TMath::Sin(0.1814) + epcrec*TMath::Cos(0.1814)>>hMean","nsite == 1 && nWires == 3");
  t->Draw("csc*TMath::Sin(0.2036) + epcrec*TMath::Cos(0.2036)>>hMulti","nsite > 1 || nsite == 1 && nWires > 2");

  TF1 *fit = new TF1("fit",fitFunction,850,1300,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameters(800,1100,100,0.8);

  fit->SetParLimits(1,1000,1200);
  fit->SetParLimits(2,20,200);
  fit->SetParLimits(3,0.05,2.0);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunction,850,1350,4);
  //TF1 *fitMulti = new TF1("fitMulti","gaus",530,650);
  fitMulti->SetParNames("A","E","#sigma","A2");

  fitMulti->SetParameters(800,1200,100,0.2);

  fitMulti->SetParLimits(1,1000,1300);
  fitMulti->SetParLimits(2,20,200);
  fitMulti->SetParLimits(3,0.05,2.0);

  hMulti->Fit("fitMulti","r");

  double par[4];
  double *parErr;
  fit->GetParameters(par);
  parErr = fit->GetParErrors();

  PeakPos1[ID] = par[1];
  PeakPosErr1[ID] = parErr[1];

  double parMulti[4];
  double *parErrMulti;
  fitMulti->GetParameters(parMulti);
  parErrMulti = fitMulti->GetParErrors();

  PeakPosMulti[ID] = parMulti[1];
  PeakPosErrMulti[ID] = parErrMulti[1];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << ")" << "  (" << PeakPosMulti[ID] << " +- " << PeakPosErrMulti[ID] << ")" << endl;
  
  c1->cd();
  hMean->Draw("EZP");
  c1->Print("Cs137_ss.ps");
  
  c2->cd();
  hMulti->Draw("EZP");
  c2->Print("Cs137_ms.ps");

  /*char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/%i.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  hMean->Write();
  hMulti->Write();
  fHisto->Close();
  
  c1->cd();
  hMean->Draw("EZP");
  char oName_ss[100];
  sprintf(oName_ss,"/var/www/html/analysis/AntiCorrelation/Cs137/SingleSite/Cs137_%i.png",RunID);
  c1->SaveAs(oName_ss);
  
  c2->cd();
  hMulti->Draw("EZP");
  char oName_ms[100];
  sprintf(oName_ms,"/var/www/html/analysis/AntiCorrelation/Cs137/MultiSite/Cs137_%i.png",RunID);
  c2->SaveAs(oName_ms);*/

  hMean->Delete();
  //hN->Delete();
  //hP->Delete();
  hMulti->Delete();
  fit->Delete();
  fitMulti->Delete();

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",750,1400);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *erf1 = new TF1("erf1",errf,750,1400,3);
  erf1->SetParameters(par[3]*par[0],par[1],par[2]);
  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);

  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  gaus1->Draw("same");
  erf1->Draw("same");

  //TCanvas *c2 = new TCanvas();
  //hN->Draw();

  //TCanvas *c3 = new TCanvas();
  //hP->Draw();

  TCanvas *c4 = new TCanvas();
  hMulti->Draw();*/

  return;
}
