double *PeakPos1;
double *PeakPos2;
double *PeakPosErr1;
double *PeakPosErr2;
double *PeakPos1Multi;
double *PeakPos2Multi;
double *PeakPosErr1Multi;
double *PeakPosErr2Multi;

void ProcessRun(int RunID, int ID);

TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();

double fitFunction(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = par[5];
  //double sigma2 = sigma1*TMath::Sqrt(E2/E1);
  double sigma2 = par[6];
  double A2_erf = A2_gaus*par[3];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  double gauss2 = A2_gaus * TMath::Gaus(x[0],E2,sigma2);
  double erf2 = A2_erf * 0.5 * TMath::Erfc((x[0] - E2) / (TMath::Sqrt(2)*sigma2));

  return gauss1 + erf1 + gauss2 + erf2;
}

double fitFunctionMulti(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = 1.1358*E1;
  double sigma2 = sigma1*TMath::Sqrt(E2/E1);
  double A2_erf = A2_gaus*par[3];

  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  double gauss2 = A2_gaus * TMath::Gaus(x[0],E2,sigma2);
  double erf2 = A2_erf * TMath::Erfc((x[0] - E2) / (TMath::Sqrt(2)*sigma2));

  return gauss1 + erf1 + gauss2 + erf2;
}

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
}

double fitGauss(double *x, double *par)
{
  double A1 = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A2 = par[0]*par[3];
  double E2 = 1.1358*E1;
  double sigma2 = sigma1*TMath::Sqrt(E2/E1);

  return A1*TMath::Gaus(x[0],E1,sigma1) + A2*TMath::Gaus(x[0],E2,sigma2);
}

void Co60Fitting()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
/*  const int nrRuns = 37;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};
  double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875, 1944, 1947, 1953, 1954, 1956, 1957, 1960, 1961, 1962, 1963, 1964, 1965};
*/
/*  const int nrRuns = 25;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875};
  double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690, 1810, 1816, 1822, 1824, 1830, 1833, 1842, 1844, 1848, 1850, 1855, 1857, 1862, 1866, 1869, 1873, 1875};
*/
/*  const int nrRuns = 8;
  int runID[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690};
  double runIDCopy[nrRuns] = {1654, 1658, 1664, 1669, 1676, 1682, 1686, 1690};
*/

/*  const int nrRuns = 17;
  int runID[nrRuns] = {2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2708};
  double runIDCopy[nrRuns] = {2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2708};
*/
/*  const int nrRuns = 9;
  int runID[nrRuns] = {2526, 2538, 2543, 2566, 2578, 2596, 2608, 2620, 2634};
  double runIDCopy[nrRuns] = {2526, 2538, 2543, 2566, 2578, 2596, 2608, 2620, 2634};
*/
  const int nrRuns = 18;
  int runID[nrRuns] = {2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};
  double runIDCopy[nrRuns] = {2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};
  
  /*const int nrRuns = 23;
  int runID[nrRuns] = {2412, 2415, 2416, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};
  double runIDCopy[nrRuns] = {2412, 2415, 2416, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};*/

  PeakPos1 = new double[nrRuns];
  PeakPos2 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];
  PeakPosErr2 = new double[nrRuns];
  PeakPos1Multi = new double[nrRuns];
  PeakPos2Multi = new double[nrRuns];
  PeakPosErr1Multi = new double[nrRuns];
  PeakPosErr2Multi = new double[nrRuns];
  
  c1->Print("Co60_ss.ps[");
  c2->Print("Co60_ms.ps[");
  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }
  c1->Print("Co60_ss.ps]");
  c2->Print("Co60_ms.ps]");

  cout << "Mean1 (single site)  = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean2 (single site)  = " << TMath::Mean(nrRuns,PeakPos2) << " +- " << TMath::RMS(nrRuns,PeakPos2) << endl;
  cout << "Mean1 (multi site) = " << TMath::Mean(nrRuns,PeakPos1Multi) << " +- " << TMath::RMS(nrRuns,PeakPos1Multi) << endl;
  cout << "Mean2 (multi site) = " << TMath::Mean(nrRuns,PeakPos2Multi) << " +- " << TMath::RMS(nrRuns,PeakPos2Multi) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2,0,PeakPosErr2);
  TGraphErrors *gr3 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1Multi,0,PeakPosErr1Multi);
  TGraphErrors *gr4 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2Multi,0,PeakPosErr2Multi);

  gr1->SetTitle("Co60 peak1 (single site)");
  gr2->SetTitle("Co60 peak2 (single site)");
  gr3->SetTitle("Co60 peak1 (multi site)");
  gr4->SetTitle("Co60 peak2 (multi site)");

  gr1->SetName("Co60_p1_single");
  gr2->SetName("Co60_p2_single");
  gr3->SetName("Co60_p1_multi");
  gr4->SetName("Co60_p2_multi");

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);
  gr2->SetMarkerColor(kRed);

  gr1->GetYaxis()->SetRangeUser(1700,2100);
  gr1->GetXaxis()->SetLimits(1570,2720);
  gr1->SetTitle("Co60 Peaks (single site)");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.6);

  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.6);
  gr4->SetMarkerColor(kRed);

  gr3->GetYaxis()->SetRangeUser(1800,2200);
  gr3->GetXaxis()->SetLimits(1570,2720);
  gr3->SetTitle("Co60 Peaks (multi site)");
  gr3->GetXaxis()->SetTitle("run number");
  gr3->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");
  gr2->Draw("Psame");

  TCanvas *c2 = new TCanvas();
  gr3->Draw("AP");
  gr4->Draw("Psame");

  /*TFile *fOut = new TFile("../analysis/Co60_PeakPosition_alpha_220212.root","RECREATE");
  gr1->Write();
  gr2->Write();
  gr3->Write();
  gr4->Write();

  fOut->Close();*/

  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  //sprintf(FileName,"../analysis/V2/%i_noReclustering_Fiducial.root",RunID);
  //sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",RunID);
  sprintf(FileName,"../../EnergyCalibration/analysis/V11/%i_noReclustering_Fiducial.root",RunID);

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
  
  TH1F *hMean = new TH1F("hMean",hName_ss,100,1000,3500);
  TH1F *hMulti = new TH1F("hMulti",hName_ms,100,1000,3500);

  //t->Draw("(csc*0.3148-119.5)*TMath::Sin(0.43218) + (ecrec*0.9647+67.4)*TMath::Cos(0.43218)>>hMean","nsite == 1");
  //t->Draw("ecrec>>hMean","nsite == 1 && TMath::Abs(zc) < 150");
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

  TF1 *fit = new TF1("fit",fitFunction,1650,2300,7);
  fit->SetParNames("A1","E1","#sigma","R1","R2","E2","#sigma2");

  fit->SetParameters(1000,1880,75,0.1,0.8,2140,75);

  fit->SetParLimits(1,1850,1950);
  //fit->SetParLimits(1,1000,1250);
  fit->SetParLimits(2,10,1020);
  fit->SetParLimits(3,0.05,1.0);
  //fit->SetParLimits(3,0.2,5);
  fit->SetParLimits(4,0.6,1.0);
  fit->SetParLimits(5,2100,2300);
  fit->SetParLimits(6,10,130);
  //fit->SetParLimits(5,1200,1400);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunction,1700,2400,7);
  //TF1 *fitMulti = new TF1("fitMulti",fitGauss,1100,1400,4);
  //fitMulti->SetParameters(400,1180,75,0.8);
  fitMulti->SetParNames("A1","E1","#sigma","R1","R2","E2","#sigma2");

  fitMulti->SetParameters(1000,2050,100,0.2,0.8,2350,100);

  //fitMulti->SetParLimits(0,100,2000);
  //fitMulti->SetParLimits(1,1100,1250);
  //fitMulti->SetParLimits(2,10,120);
  //fitMulti->SetParLimits(3,0.5,1);
  //fitMulti->SetParLimits(1,1100,1250);
  fitMulti->SetParLimits(1,1900,2200);
  fitMulti->SetParLimits(2,10,200);
  fitMulti->SetParLimits(3,0.05,1.0);
  fitMulti->SetParLimits(4,0.5,1.5);
  fitMulti->SetParLimits(5,2200,2400);
  fitMulti->SetParLimits(6,10,130);

  hMulti->Fit("fitMulti","r");

  double par[7];
  double *parErr;
  fit->GetParameters(par);
  parErr = fit->GetParErrors();

  PeakPos1[ID] = par[1];
  PeakPos2[ID] = par[5];
  PeakPosErr1[ID] = parErr[1];
  PeakPosErr2[ID] = parErr[5];

  double parMulti[7];
  double *parErrMulti;
  fitMulti->GetParameters(parMulti);
  parErrMulti = fitMulti->GetParErrors();

  PeakPos1Multi[ID] = parMulti[1];
  PeakPos2Multi[ID] = parMulti[5];
  PeakPosErr1Multi[ID] = parErrMulti[1];
  PeakPosErr2Multi[ID] = parErrMulti[5];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << "  " << PeakPos2[ID] << " +- " << PeakPosErr2[ID] << ")  (" << PeakPos1Multi[ID] << " +- " << PeakPosErr1Multi[ID] << "  " << PeakPos2Multi[ID] << " +- " << PeakPosErr2Multi[ID] << ")" << endl;
  
  c1->cd();
  hMean->Draw("EZP");
  c1->Print("Co60_ss.ps");
  
  c2->cd();
  hMulti->Draw("EZP");
  c2->Print("Co60_ms.ps");

  /*char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/%i.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  hMean->Write();
  hMulti->Write();
  fHisto->Close();
  
  c1->cd();
  hMean->Draw("EZP");
  char oName_ss[100];
  sprintf(oName_ss,"/var/www/html/analysis/AntiCorrelation/Co60/SingleSite/Co60_%i.png",RunID);
  c1->SaveAs(oName_ss);
  
  c2->cd();
  hMulti->Draw("EZP");
  char oName_ms[100];
  sprintf(oName_ms,"/var/www/html/analysis/AntiCorrelation/Co60/MultiSite/Co60_%i.png",RunID);
  c2->SaveAs(oName_ms);*/

  hMean->Delete();
  //hN->Delete();
  //hP->Delete();
  fit->Delete();
  fitMulti->Delete();
  hMulti->Delete();

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",1300,2500);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *gaus2 = new TF1("gaus2","[0]*TMath::Gaus(x,[1],[2])",1300,2500);
  gaus2->SetParameters(par[0]*par[4],par[5],par[6]);

  TF1 *erf1 = new TF1("erf1",errf,1300,2500,3);
  erf1->SetParameters(par[0]*par[3],par[1],par[2]);

  TF1 *erf2 = new TF1("erf2",errf,1300,2500,3);
  erf2->SetParameters(par[0]*par[4]*par[3],par[5],par[6]);

  TF1 *gaus1Multi = new TF1("gaus1Multi","gaus",1300,2500);
  gaus1Multi->SetParameters(parMulti[0],parMulti[1],parMulti[2]);

  TF1 *gaus2Multi = new TF1("gaus2Multi","gaus",1300,2500);
  gaus2Multi->SetParameters(parMulti[0]*parMulti[4],parMulti[5],parMulti[6]);

  gaus1->SetLineWidth(1);
  gaus1->SetLineColor(kRed);
  gaus2->SetLineWidth(1);
  gaus2->SetLineColor(kBlue);
  erf1->SetLineWidth(1);
  erf1->SetLineColor(kRed);
  erf1->SetLineStyle(2);
  erf2->SetLineWidth(1);
  erf2->SetLineColor(kBlue);
  erf2->SetLineStyle(2);

  gaus1Multi->SetLineWidth(1);
  gaus1Multi->SetLineColor(kRed);
  gaus2Multi->SetLineWidth(1);
  gaus2Multi->SetLineColor(kBlue);

  TCanvas *c1 = new TCanvas();
  hMean->Draw();
  gaus1->Draw("same");
  gaus2->Draw("same");
  erf1->Draw("same");
  erf2->Draw("same");

  //TCanvas *c2 = new TCanvas();
  //hN->Draw();

  //TCanvas *c3 = new TCanvas();
  //hP->Draw();

  TCanvas *c4 = new TCanvas();
  hMulti->Draw();
  gaus1Multi->Draw("same");
  gaus2Multi->Draw("same");*/

  return;
}
