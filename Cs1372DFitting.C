double *PeakPos1;
double *PeakPosErr1;
double *PeakPos2;
double *PeakPosErr2;

Double_t g2(Double_t *x, Double_t *par) {
   Double_t r1 = Double_t((x[0]-par[1])/par[2]);
   Double_t r2 = Double_t((x[1]-par[3])/par[4]);
   return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}
double  fun2(double *x, double *par) {
   double A = par[0];
   double x_0 = par[1];
   double y_0 = par[2];
   double sigma_x = par[3];
   double sigma_y = par[4];
   double a = TMath::Cos(par[5])**2 / (2*sigma_x**2) + TMath::Sin(par[5])**2 / (2*sigma_y**2);
   double b = -1.0 * TMath::Sin(2*par[5]) / (4*sigma_x**2) + TMath::Sin(2*par[5]) / (4*sigma_y**2);
   double c = TMath::Sin(par[5])**2 / (2*sigma_x**2) + TMath::Cos(par[5])**2 / (2*sigma_y**2);

   double val = A*TMath::Exp(-1.0*(a*(x[0] - x_0)**2 + 2*b*(x[0] - x_0)*(x[1] - y_0) + c*(x[1] - y_0)**2));

   //Double_t *p1 = &par[0];
   //Double_t *p2 = &par[5];
   //Double_t result = g2(x,p1);
   //return result;

   return val;
}

void Cs1372DFitting()
{
  gStyle->SetPalette(1);

/*  const int nrRuns = 3;
  int runID[nrRuns] = {2449, 2469, 2473};
  double runIDCopy[nrRuns] = {2449, 2469, 2473};
*/
  const int nrRuns = 5;
  int runID[nrRuns] = {2410, 2449, 2450, 2469, 2473};
  double runIDCopy[nrRuns] = {2410, 2449, 2450, 2469, 2473};

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPos2 = new double[nrRuns];
  PeakPosErr2 = new double[nrRuns];

  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

  cout << "Mean (charge) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (scintillation) = " << TMath::Mean(nrRuns,PeakPos2) << " +- " << TMath::RMS(nrRuns,PeakPos2) << endl;

  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2,0,PeakPosErr2);
  
  gr1->SetTitle("Cs137 peak (charge)");
  gr2->SetTitle("Cs137 peak (scintillation)");
  
  gr1->SetName("Cs137_charge");
  gr2->SetName("Cs137_scintillation");
  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);
  
  gr1->GetYaxis()->SetRangeUser(550,600);
  gr1->GetXaxis()->SetLimits(1560,2480);
  gr1->SetTitle("Cs137 Peak (charge)");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");
  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);
  
  gr2->GetYaxis()->SetRangeUser(2400,2600);
  gr2->GetXaxis()->SetLimits(1560,2480);
  gr2->SetTitle("Cs137 Peak (scintillation)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");
  
  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");
  
  /*TFile *fOut = new TFile("../analysis/Cs137_2D_PeakPosition.root","RECREATE");
  gr1->Write();
  gr2->Write();
  
  fOut->Close();*/
  
  return;
}

void ProcessRun(int RunID, int ID)
{
  char FileName[100];
  sprintf(FileName,"../../EnergyCalibration/analysis/alpha/V2/%i_noReclustering_Fiducial.root",RunID);

  cout << "Processing run " << RunID << "..." << endl;

  TFile *f = new TFile(FileName,"READ");

  TTree *t = (TTree*)f->Get("t");

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,1500,200,0,5000);

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");

  TF2 *f2 = new TF2("f2",fun2,500,700,2000,3000,6);
  f2->SetParameters(40,600,2500,70,200,0.0);
  f2->SetParNames("A","E_I","E_S","#sigma_I","#sigma_S","#Theta");

  f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,540,640);
  f2->SetParLimits(2,2200,2700);
  f2->SetParLimits(3,10,100);
  f2->SetParLimits(4,20,600);

  h2DSingle->Fit("f2","r");

  double par[6];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();

  /*char FormatedOutput[1000];
  sprintf(FormatedOutput,"%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\n",RunID,0,par[1],par[2],par[3],par[4],par[5],parErr[1],parErr[2],parErr[3],parErr[4],parErr[5]);
  
  FILE *ResultFile;
  ResultFile = fopen("Cs137.dat","a");
  fputs(FormatedOutput,ResultFile);
  fclose(ResultFile);*/
  
  PeakPos1[ID] = par[1];
  PeakPosErr1[ID] = parErr[1];
  PeakPos2[ID] = par[2];
  PeakPosErr2[ID] = parErr[2];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << ")" << "  (" << PeakPos2[ID] << " +- " << PeakPosErr2[ID] << ")" << endl;

  /*char fHistoName[100];
  sprintf(fHistoName,"../analysis/Fitting/2D/%i.root",RunID);
  TFile *fHisto = new TFile(fHistoName,"RECREATE");
  h2DSingle->Write();
  fHisto->Close();*/
  
  //TCanvas *c1 = new TCanvas();
  //h2DSingle->Draw("colz");
  //f2->Draw("cont3 same");

  return;
}
