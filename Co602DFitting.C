double *PeakPos1;
double *PeakPosErr1;
double *PeakPos2;
double *PeakPosErr2;
double *PeakPos3;
double *PeakPosErr3;
double *PeakPos4;
double *PeakPosErr4;

Double_t g2(Double_t *x, Double_t *par) {
   double A = par[0];
   double x_0 = par[1];
   double y_0 = par[2];
   double sigma_x = par[3];
   double sigma_y = par[4];
   double a = TMath::Cos(par[5])**2 / (2*sigma_x**2) + TMath::Sin(par[5])**2 / (2*sigma_y**2);
   double b = -1.0 * TMath::Sin(2*par[5]) / (4*sigma_x**2) + TMath::Sin(2*par[5]) / (4*sigma_y**2);
   double c = TMath::Sin(par[5])**2 / (2*sigma_x**2) + TMath::Cos(par[5])**2 / (2*sigma_y**2);

   double val = A*TMath::Exp(-1.0*(a*(x[0] - x_0)**2 + 2*b*(x[0] - x_0)*(x[1] - y_0) + c*(x[1] - y_0)**2));

   return val;
}

Double_t fun2(Double_t *x, Double_t *par) {
   Double_t *p1 = &par[0];
   Double_t *p2 = &par[6];
   Double_t result = g2(x,p1) + g2(x,p2);
   return result;
}

void Co602DFitting()
{
  gStyle->SetPalette(1);

  const int nrRuns = 24;
  int runID[nrRuns] = {2412, 2415, 2416, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2674, 2683, 2689, 2708};
  double runIDCopy[nrRuns] = {2412, 2415, 2416, 2480, 2496, 2526, 2538, 2543, 2555, 2566, 2578, 2596, 2608, 2620, 2634, 2635, 2640, 2646, 2653, 2667, 2674, 2683, 2689, 2708};

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPos2 = new double[nrRuns];
  PeakPosErr2 = new double[nrRuns];

  PeakPos3 = new double[nrRuns];
  PeakPosErr3 = new double[nrRuns];

  PeakPos4 = new double[nrRuns];
  PeakPosErr4 = new double[nrRuns];

  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

  cout << "Mean peak1 (charge) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean peak1 (scintillation) = " << TMath::Mean(nrRuns,PeakPos2) << " +- " << TMath::RMS(nrRuns,PeakPos2) << endl;
  cout << "Mean peak2 (charge) = " << TMath::Mean(nrRuns,PeakPos3) << " +- " << TMath::RMS(nrRuns,PeakPos3) << endl;
  cout << "Mean peak2 (scintillation) = " << TMath::Mean(nrRuns,PeakPos4) << " +- " << TMath::RMS(nrRuns,PeakPos4) << endl;
  
  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2,0,PeakPosErr2);
  TGraphErrors *gr3 = new TGraphErrors(nrRuns,runIDCopy,PeakPos3,0,PeakPosErr3);
  TGraphErrors *gr4 = new TGraphErrors(nrRuns,runIDCopy,PeakPos4,0,PeakPosErr4);
  
  gr1->SetTitle("Co60 peak1 (charge)");
  gr2->SetTitle("Co60 peak1 (scintillation)");
  gr3->SetTitle("Co60 peak2 (charge)");
  gr4->SetTitle("Co60 peak2 (scintillation)");
  
  gr1->SetName("Co60_p1_charge");
  gr2->SetName("Co60_p1_scint");
  gr3->SetName("Co60_p2_charge");
  gr4->SetName("Co60_p2_scint");
  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);
  
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.6);
  gr3->SetMarkerColor(kRed);
  
  gr1->GetYaxis()->SetRangeUser(1700,2100);
  gr1->GetXaxis()->SetLimits(1570,2720);
  gr1->SetTitle("Co60 Peaks (charge)");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");
  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);
  
  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.6);
  gr4->SetMarkerColor(kRed);
  
  gr3->GetYaxis()->SetRangeUser(3000,5000);
  gr3->GetXaxis()->SetLimits(1570,2720);
  gr3->SetTitle("Co60 Peaks (scintillation)");
  gr3->GetXaxis()->SetTitle("run number");
  gr3->GetYaxis()->SetTitle("peak position");
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");
  gr3->Draw("Psame");
  
  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");
  gr4->Draw("Psame");
  
  /*TFile *fOut = new TFile("../analysis/Co60_2D_PeakPosition.root","RECREATE");
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
  sprintf(FileName,"../../EnergyCalibration/analysis/alpha/V2/%i_noReclustering_Fiducial.root",RunID);

  cout << "Processing run " << RunID << "..." << endl;

  TFile *f = new TFile(FileName,"READ");

  TTree *t = (TTree*)f->Get("t");

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,2000,200,0,10000);

  t->Draw("csc:epcrec>>h2DSingle","nsite == 1");

  TF2 *f2 = new TF2("f2",fun2,1000,1500,3500,6000,12);
  double params[12] = {40,1100,4000,70,200,0.0,40,1300,4600,70,200,0.0};
  f2->SetParameters(params);
  //f2->SetParNames("A1","E_I1","E_S1","#sigma_I1","#sigma_S1","#Theta1","A2","E_I2","E_S2","#sigma_I2","#sigma_S2","#Theta2");

  f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,1000,1200);
  f2->SetParLimits(2,3600,4400);
  f2->SetParLimits(3,10,100);
  f2->SetParLimits(4,20,600);
  f2->SetParLimits(6,0,1000);
  f2->SetParLimits(7,1250,1400);
  f2->SetParLimits(8,4400,5000);
  f2->SetParLimits(9,10,100);
  f2->SetParLimits(10,20,600);

  h2DSingle->Fit("f2","r");

  double par[12];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();

  /*char FormatedOutput[2000];
  sprintf(FormatedOutput,"%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\n",RunID,0,par[1],par[2],par[3],par[4],par[5],par[7],par[8],par[9],par[10],par[11],parErr[1],parErr[2],parErr[3],parErr[4],parErr[5],parErr[7],parErr[8],parErr[9],parErr[10],parErr[11]);
  
  FILE *ResultFile;
  ResultFile = fopen("Co60.dat","a");
  fputs(FormatedOutput,ResultFile);
  fclose(ResultFile);*/
  
  PeakPos1[ID] = par[1];
  PeakPosErr1[ID] = parErr[1];
  PeakPos2[ID] = par[2];
  PeakPosErr2[ID] = parErr[2];
  PeakPos3[ID] = par[7];
  PeakPosErr3[ID] = parErr[7];
  PeakPos4[ID] = par[8];
  PeakPosErr4[ID] = parErr[8];

  cout << "End of run " << RunID << "  (" << PeakPos1[ID] << " +- " << PeakPosErr1[ID] << ")" << "  (" << PeakPos2[ID] << " +- " << PeakPosErr2[ID] << ")" << endl;
  cout << "End of run " << RunID << "  (" << PeakPos3[ID] << " +- " << PeakPosErr3[ID] << ")" << "  (" << PeakPos4[ID] << " +- " << PeakPosErr4[ID] << ")" << endl;
  
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
