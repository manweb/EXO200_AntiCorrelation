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

double errf(double *x, double *par)
{
  return par[0] * 0.5 * TMath::Erfc((x[0] - par[1]) / (TMath::Sqrt(2)*par[2]));
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

void Th228Fitting()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
/*  const int nrRuns = 36;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1889, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
*/

// modified list
/*  const int nrRuns = 33;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1702, 1709, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1631, 1638, 1642, 1647, 1702, 1709, 1719, 1731, 1736, 1773, 1786, 1793, 1802, 1915, 1916, 1917, 1918, 1921, 1923, 1925, 1926, 1929, 1930, 1931, 1932, 1937};
*/
/*  const int nrRuns = 15;
  int runID[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};
  double runIDCopy[nrRuns] = {1573, 1587, 1595, 1599, 1605, 1619, 1627, 1638, 1642, 1647, 1697, 1702, 1709, 1714, 1719};
*/

/*  const int nrRuns = 28;
  int runID[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865};
  double runIDCopy[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865};*/
  //int runID[nrRuns] = {2417, 2421, 2422, 2424, 2426};
  //double runIDCopy[nrRuns] = {2417, 2421, 2422, 2424, 2426};

  const int nrRuns = 8;
  int runID[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  double runIDCopy[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};

  /*const int nrRuns = 36;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};*/
  
  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPosMulti = new double[nrRuns];
  PeakPosErrMulti = new double[nrRuns];

  c1->Print("Th228_ss.ps[");
  c2->Print("Th228_ms.ps[");
  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }
  c1->Print("Th228_ss.ps]");
  c2->Print("Th228_ms.ps]");

  cout << "Mean (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (multi site) = " << TMath::Mean(nrRuns,PeakPosMulti) << " +- " << TMath::RMS(nrRuns,PeakPosMulti) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPosMulti,0,PeakPosErrMulti);

  gr1->SetTitle("Th228 peak (single site)");
  gr2->SetTitle("Th228 peak (multi site)");

  gr1->SetName("Th228_single");
  gr2->SetName("Th228_multi");

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetYaxis()->SetRangeUser(4000,4500);
  gr1->GetXaxis()->SetLimits(2400,3150);
  gr1->SetTitle("Th228 Peak");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetYaxis()->SetRangeUser(4000,4500);
  gr2->GetXaxis()->SetLimits(2400,3150);
  gr2->SetTitle("Th228 Peak (multi site)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");

  TCanvas *c1_2 = new TCanvas();
  gr1->Draw("AP");

  TCanvas *c2_2 = new TCanvas();
  gr2->Draw("AP");

  /*TFile *fOut = new TFile("../analysis/Th228_PeakPosition_alpha_220212.root","RECREATE");

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

  TH1F *hMean = new TH1F("hMean",hName_ss,100,3500,5000);
  TH1F *hMulti = new TH1F("hMulti",hName_ms,100,3500,5000);
  //TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,15000);
  //TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,15000);

  //t->Draw("ecrec_gc>>hMean","nsite == 1");
  //t->Draw("ecrec>>hMean","nsite == 1 && TMath::Abs(zc) < 150");
  //t->Draw("ecrec_gc>>hMulti","nsite == 2");
  //t->Draw("(csc*0.3148-119.5)*TMath::Sin(0.5369) + (ecrec*0.9647+67.4)*TMath::Cos(0.5369)>>hMean","nsite == 1");
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

  TF1 *fit = new TF1("fit",fitFunction,4000,4400,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameters(1000,4300,100,0.2);

  fit->SetParLimits(1,4100,4400);
  fit->SetParLimits(2,20,200);
  fit->SetParLimits(3,0.01,0.8);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunction,4130,4550,4);
  //TF1 *fitMulti = new TF1("fitMulti","gaus",2550,2850);
  fitMulti->SetParNames("A","E","#sigma","A2");

  fitMulti->SetParameters(1000,4500,100,0.1);

  fitMulti->SetParLimits(1,4300,4600);
  fitMulti->SetParLimits(2,20,200);
  fitMulti->SetParLimits(3,0.01,0.8);

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
  c1->Print("Th228_ss.ps");
  
  c2->cd();
  hMulti->Draw("EZP");
  c2->Print("Th228_ms.ps");

  /*char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/%i.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  hMean->Write();
  hMulti->Write();
  fHisto->Close();
  
  c1->cd();
  hMean->Draw("EZP");
  char oName_ss[100];
  sprintf(oName_ss,"/var/www/html/analysis/AntiCorrelation/Th228/SingleSite/Th228_%i.png",RunID);
  c1->SaveAs(oName_ss);
  
  c2->cd();
  hMulti->Draw("EZP");
  char oName_ms[100];
  sprintf(oName_ms,"/var/www/html/analysis/AntiCorrelation/Th228/MultiSite/Th228_%i.png",RunID);
  c2->SaveAs(oName_ms);*/

  hMean->Delete();
  //hN->Delete();
  //hP->Delete();
  hMulti->Delete();
  fit->Delete();
  fitMulti->Delete();

/*  TF1 *gaus1 = new TF1("gaus1","[0]*TMath::Gaus(x,[1],[2])",3800,4600);
  gaus1->SetParameters(par[0],par[1],par[2]);

  TF1 *erf1 = new TF1("erf1",errf,3800,4600,3);
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
