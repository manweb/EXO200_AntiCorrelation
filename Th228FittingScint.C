double *PeakPos1;
double *PeakPosErr1;
double *PeakPosMulti;
double *PeakPosErrMulti;

void ProcessRun(int RunID, int ID);

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

void Th228FittingScint()
{
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

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPosMulti = new double[nrRuns];
  PeakPosErrMulti = new double[nrRuns];

  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

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
  gr1->GetXaxis()->SetLimits(2400,2900);
  gr1->SetTitle("Th228 Peak");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetYaxis()->SetRangeUser(4000,4500);
  gr2->GetXaxis()->SetLimits(2400,2900);
  gr2->SetTitle("Th228 Peak (multi site)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");

  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");

  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");

  TFile *fOut = new TFile("../analysis/Th228_scint_PeakPosition_all.root","RECREATE");

  gr1->Write();
  gr2->Write();

  fOut->Close();

  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  sprintf(FileName,"../../EnergyCalibration/analysis/V6/%i_noReclustering_Fiducial.root",RunID);
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

  TH1F *hMean = new TH1F("hMean","Combined energy spectrum (single site)",200,0,11000);
  TH1F *hMulti = new TH1F("hMulti","Combined energy spectrum (multi site)",200,0,11000);
  //TH1F *hMean = new TH1F("hMean","Charge Energy Spectrum (mean purity)",200,0,15000);
  //TH1F *hMulti = new TH1F("hMulti","Multi Site Charge Energy Spectrum (mean purity)",200,0,15000);

  //t->Draw("ecrec_gc>>hMean","nsite == 1");
  //t->Draw("ecrec>>hMean","nsite == 1 && TMath::Abs(zc) < 150");
  //t->Draw("ecrec_gc>>hMulti","nsite == 2");
  t->Draw("csc>>hMean","nsite == 1");
  t->Draw("csc>>hMulti","nsite > 1");

  TF1 *fit = new TF1("fit",fitFunction,3890,4600,4);
  fit->SetParNames("A1","E1","#sigma","A2");

  fit->SetParameters(1000,4250,100,0.8);

  fit->SetParLimits(1,3900,4600);
  fit->SetParLimits(2,50,200);
  fit->SetParLimits(3,0.05,2.0);

  hMean->Fit("fit","r");

  TF1 *fitMulti = new TF1("fitMulti",fitFunction,3890,4500,4);
  //TF1 *fitMulti = new TF1("fitMulti","gaus",2550,2850);
  fitMulti->SetParNames("A","E","#sigma","A2");

  fitMulti->SetParameters(1000,4250,100,0.2);

  fitMulti->SetParLimits(1,3900,4600);
  fitMulti->SetParLimits(2,50,200);
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

  char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/%i_scint.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  hMean->Write();
  hMulti->Write();
  fHisto->Close();

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
