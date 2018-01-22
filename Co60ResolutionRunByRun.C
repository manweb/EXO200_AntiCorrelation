double *Res1;
double *Res1Err;
double *Res2;
double *Res2Err;
double *Peak1Pos;
double *Peak1PosErr;
double *Peak2Pos;
double *Peak2PosErr;

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

void Co60ResolutionRunByRun()
{
  const int nrRuns = 17;
  int runID[nrRuns] = {2526, 2538, 2543, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};
  double runIDCopy[nrRuns] = {2526, 2538, 2543, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2683, 2689, 2708};
  
  Res1 = new double[nrRuns];
  Res1Err = new double[nrRuns];
  Res2 = new double[nrRuns];
  Res2Err = new double[nrRuns];
  Peak1Pos = new double[nrRuns];
  Peak1PosErr = new double[nrRuns];
  Peak2Pos = new double[nrRuns];
  Peak2PosErr = new double[nrRuns];
  
  for (int i = 0; i < nrRuns; i++) {
    ProcessRun(runID[i],i);
  }
  
  cout << "Mean resolution 1 = " << TMath::Mean(nrRuns,Res1) << " +- " << TMath::RMS(nrRuns,Res1) << endl;
  cout << "Mean resolution 2 = " << TMath::Mean(nrRuns,Res2) << " +- " << TMath::RMS(nrRuns,Res2) << endl;
  
  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,Res1,0,Res1Err);
  
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("resolution");
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.8);
  
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,Res2,0,Res2Err);
  
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("resolution");
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.8);
  gr2->SetMarkerColor(kRed);
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AZP");
  gr2->Draw("ZPsame");
  
  TGraphErrors *gr3 = new TGraphErrors(nrRuns,runIDCopy,Peak1Pos,0,Peak1PosErr);
  
  gr3->GetXaxis()->SetTitle("run number");
  gr3->GetYaxis()->SetTitle("peak position");
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.8);
  
  TGraphErrors *gr4 = new TGraphErrors(nrRuns,runIDCopy,Peak2Pos,0,Peak2PosErr);
  
  gr4->GetXaxis()->SetTitle("run number");
  gr4->GetYaxis()->SetTitle("peak position");
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.8);
  gr4->SetMarkerColor(kRed);
  
  TCanvas *c2 = new TCanvas();
  gr3->Draw("AZP");
  gr4->Draw("ZPsame");
  
  return;
}

void ProcessRun(int RunID, int ID)
{
  gStyle->SetPalette(1);
  
  char FileName[100];
  //sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",RunID);
  sprintf(FileName,"../../EnergyCalibration/analysis/alpha/V1/%i_noReclustering_Fiducial.root",RunID);
  
  TFile *f = new TFile(FileName,"READ");
  TTree *t = (TTree*)f->Get("t");
  
  cout << "Processing run " << RunID << endl;
  cout << "    fitting .... ";
  
  // Reprocessed data (01/30/12)
  //double theta = 0.1916;
  //double theta = 0.2164;
  
  // Reprocessed data (02/07/12)
  double theta = 0.1867;
  //double theta = 0.2084;
  
  char hName[100];
  sprintf(hName,"Th228 (single site, #theta = %.4f)",theta);
  TH1F *hMin = new TH1F("hMin",hName,200,0,3500);
  char cmd[100];
  //sprintf(cmd,"((csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f))*0.7192 + 49.73 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.6138 - 3.328 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.5791 + 1.925 >> hMin",theta,theta);
  
  // Reprocessed data (01/30/12)
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6117 + 2.443>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5826 + 8.801 >> hMin",theta,theta);
  
  // Reprocessed data (02/07/12)
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598>> hMin",theta,theta);
  
  // Reprocessed data (02/18/12)
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6133 + 2.057>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5874 + 7.283>> hMin",theta,theta);
  
  t->Draw(cmd,"nsite == 1");
  
  double cscOff = 0;
  
  /*  double fitRange1 = (3500-cscOff)*TMath::Sin(theta) + 900*TMath::Cos(theta);
   double fitRange2 = (6000-cscOff)*TMath::Sin(theta) + 1500*TMath::Cos(theta);
   double Peak1 = (4300-cscOff)*TMath::Sin(theta) + 1120*TMath::Cos(theta);
   double Peak2 = (4800-cscOff)*TMath::Sin(theta) + 1300*TMath::Cos(theta);
   */
  double fitRange1 = 1000;
  double fitRange2 = 1500;
  double Peak1 = 1173;
  double Peak2 = 1332;
  
  TF1 *fitMin = new TF1("fitMin",fitFunction,fitRange1,fitRange2,7);
  fitMin->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  
  fitMin->SetParameters(3000,Peak1,75,0.6,0.8,Peak2,75);
  fitMin->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
  fitMin->SetParLimits(2,10,120);
  fitMin->SetParLimits(3,0.05,1.0);
  fitMin->SetParLimits(4,0.5,1.5);
  fitMin->SetParLimits(5,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
  fitMin->SetParLimits(6,10,130);
  
  hMin->Fit("fitMin","rq");
  
  double E1 = fitMin->GetParameter(1);
  double sigma1 = fitMin->GetParameter(2);
  double E2 = fitMin->GetParameter(5);
  double sigma2 = fitMin->GetParameter(6);
  double E1_err = fitMin->GetParError(1);
  double sigma1_err = fitMin->GetParError(2);
  double E2_err = fitMin->GetParError(5);
  double sigma2_err = fitMin->GetParError(6);
  
  cout << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2) << "    " << sigma2 / E2 << " +- " << TMath::Sqrt(sigma2_err**2 / E2**2 + (sigma2/E2**2)**2 * E2_err**2) << endl;
  
  Res1[ID] = sigma1 / E1;
  Res1Err[ID] = TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2);
  Res2[ID] = sigma2 / E2;
  Res2Err[ID] = TMath::Sqrt(sigma2_err**2 / E2**2 + (sigma2/E2**2)**2 * E2_err**2);
  Peak1Pos[ID] = E1;
  Peak1PosErr[ID] = E1_err;
  Peak2Pos[ID] = E2;
  Peak2PosErr[ID] = E2_err;
  
  return;
}