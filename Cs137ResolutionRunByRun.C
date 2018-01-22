double *Res;
double *ResErr;
double *PeakPos;
double *PeakPosErr;

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

void Cs137ResolutionRunByRun()
{
  const int nrRuns = 3;
  int runID[nrRuns] = {2450, 2469, 2473};
  double runIDCopy[nrRuns] = {2450, 2469, 2473};
  
  Res = new double[nrRuns];
  ResErr = new double[nrRuns];
  PeakPos = new double[nrRuns];
  PeakPosErr = new double[nrRuns];
  
  for (int i = 0; i < nrRuns; i++) {
    ProcessRun(runID[i],i);
  }
  
  cout << "Mean reolution = " << TMath::Mean(nrRuns,Res) << " +- " << TMath::RMS(nrRuns,Res) << endl;
  
  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,Res,0,ResErr);
  
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("resolution");
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.8);
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AZP");
  
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos,0,PeakPosErr);
  
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.8);
  
  TCanvas *c2 = new TCanvas();
  gr2->Draw("AZP");
  
  return;
}

void ProcessRun(int RunID, int ID)
{
  gStyle->SetPalette(1);
  
  char FileName[100];
  sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",RunID);
  
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
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598>> hMin",theta,theta);
  
  t->Draw(cmd,"nsite == 1");
  
  double cscOff = 0;
  
/*  double fitRange1 = (1800-cscOff)*TMath::Sin(theta) + 450*TMath::Cos(theta);
  double fitRange2 = (3300-cscOff)*TMath::Sin(theta) + 750*TMath::Cos(theta);
  double Peak = (2500-cscOff)*TMath::Sin(theta) + 600*TMath::Cos(theta);
*/
  double fitRange1 = 500;
  double fitRange2 = 800;
  double Peak = 662;
  
  TF1 *fitMin = new TF1("fitMin",fitFunction,fitRange1,fitRange2,4);
  
  fitMin->SetParNames("A1","E1","#sigma","A2");
  
  fitMin->SetParameters(10000,Peak,100,0.2);
  
  fitMin->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
  fitMin->SetParLimits(2,20,200);
  fitMin->SetParLimits(3,0.05,2.0);
  
  hMin->Fit("fitMin","rq");
  
  double E1 = fitMin->GetParameter(1);
  double sigma1 = fitMin->GetParameter(2);
  double E1_err = fitMin->GetParError(1);
  double sigma1_err = fitMin->GetParError(2);
  
  cout << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2)  << endl;
  
  Res[ID] = sigma1 / E1;
  ResErr[ID] = TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2);
  PeakPos[ID] = E1;
  PeakPosErr[ID] = E1_err;
  
  return;
}
