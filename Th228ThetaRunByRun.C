double *Theta;
double *ThetaErr;

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

void Th228ThetaRunByRun()
{
  const int nrRuns = 8;
  int runID[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  double runIDCopy[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  
  /*const int nrRuns = 13;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  const int nrRuns = 36;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2817, 2837, 2848, 2865, 2867, 2883, 2897, 2904, 2938, 2943, 2947, 2966, 2981, 2991, 3007, 3019, 3028, 3034, 3053, 3074, 3091, 3099, 3109};*/
  
  Theta = new double[nrRuns];
  ThetaErr = new double[nrRuns];
  
  for (int i = 0; i < nrRuns; i++) {
    ProcessRun(runID[i],i);
  }
  
  cout << "Mean theta = " << TMath::Mean(nrRuns,Theta) << " +- " << TMath::RMS(nrRuns,Theta) << endl;
  
  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,Theta,0,ThetaErr);
  
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("theta");
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.8);
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AZP");
  
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
  cout << "    minimizing ";
  
  const int nrHistos = 20;
  
  TH1F *hSingle[nrHistos];
  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,3500,200,0,15000);
  
  TF1 *fitSingle[nrHistos];
  
  double *angle = new double[nrHistos];
  double *res = new double[nrHistos];
  double *resErr = new double[nrHistos];
  
  double thetaRange1 = 0.13;
  double thetaRange2 = 0.24; //TMath::Pi()/9.0;
  double dTheta = thetaRange2 - thetaRange1;
  double stp = dTheta/nrHistos;
  
  for (int i = 0; i < nrHistos; i++) {
    char hSingleName[100];
    char hSingleTitle[100];
    
    double theta = thetaRange1 + stp*i;
    
    sprintf(hSingleName,"hSingle_%i",i);
    sprintf(hSingleTitle,"%.0f rotation angle",180.0/TMath::Pi()*theta);
    
    hSingle[i] = new TH1F(hSingleName,hSingleTitle,200,0,6000);
    
    char SingleCMD[200];
    sprintf(SingleCMD,"csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f) >> %s",theta,theta,hSingleName);
    
    t->Draw(SingleCMD,"nsite == 1");
    
    char fitSingleName[100];
    sprintf(fitSingleName,"fit_%i",i);
    
    double cscOff = 0;
    
    double fitRange1 = (7500-cscOff)*TMath::Sin(theta) + 2450*TMath::Cos(theta);
    double fitRange2 = (9500-cscOff)*TMath::Sin(theta) + 2850*TMath::Cos(theta);
    double Peak = (8500-cscOff)*TMath::Sin(theta) +2660*TMath::Cos(theta);
    
    fitSingle[i] = new TF1(fitSingleName,fitFunction,fitRange1,fitRange2,4);
    
    fitSingle[i]->SetParNames("A1","E1","#sigma","A2");
    
    fitSingle[i]->SetParameters(10000,Peak,100,0.2);
    
    fitSingle[i]->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
    fitSingle[i]->SetParLimits(2,50,200);
    fitSingle[i]->SetParLimits(3,0.05,2.0);
    
    fitSingle[i]->SetLineWidth(1);
    fitSingle[i]->SetLineColor(kBlue);
    
    hSingle[i]->Fit(fitSingleName,"rq");
    
    double par[4];
    fitSingle[i]->GetParameters(par);
    double par1 = fitSingle[i]->GetParameter(1);
    double par2 = fitSingle[i]->GetParameter(2);
    double *parErr = fitSingle[i]->GetParErrors();
    
    angle[i] = theta;
    res[i] = fitSingle[i]->GetParameter(2)/fitSingle[i]->GetParameter(1);
    resErr[i] = TMath::Sqrt(1.0 / par1**2 * parErr[2]**2 + (par[2] / par1**2)**2 * parErr[1]**2);
    
    cout << ".";
  }
  
  t->Draw("csc:epcrec>>h2DSingle","nsite == 1");
  
  TGraphErrors *grSingleRes = new TGraphErrors(nrHistos,angle,res,0,resErr);
  grSingleRes->SetMarkerStyle(20);
  grSingleRes->SetMarkerSize(0.8);
  
  TF1 *fitRes = new TF1("fitRes","pol2",0.135,0.24);
  fitRes->SetLineWidth(1);
  fitRes->SetLineColor(kBlue);
  
  grSingleRes->Fit("fitRes","rq");
  
  double parRes1 = fitRes->GetParameter(2);
  double parRes2 = fitRes->GetParameter(1);
  double *parResErr = fitRes->GetParErrors();
  
  TF1 *fitRes2 = new TF1("fitRes2","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",0.135,0.24);
  fitRes2->SetParameter(2,-0.5*parRes2/parRes1);
  
  grSingleRes->Fit("fitRes2","rq");
  Theta[ID] = fitRes2->GetParameter(2);
  ThetaErr[ID] = fitRes2->GetParError(2);
  
  cout << " " << Theta[ID] << endl;
  
  TCanvas *c1 = new TCanvas("c1","2D plot");
  h2DSingle->Draw("colz");
  
  TCanvas *c2 = new TCanvas("c2","Spectra at respective angles");
  c2->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
    c2->cd(i+1);
    hSingle[i]->Draw();
  }
  
  TCanvas *c3 = new TCanvas("c3","Resolution vs theta");
  grSingleRes->Draw("AZP");
  
  char oname[100];
  sprintf(oname,"../analysis/output/theta/%i.root",RunID);
  TFile *f = new TFile(oname,"RECREATE");
  
  h2DSingle->Write();
  c2->Write();
  c3->Write();
  
  f->Close();
  
  delete c1;
  delete c2;
  delete c3;
  
  return;
}
