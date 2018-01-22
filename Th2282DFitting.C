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

void Th2282DFitting()
{
  gStyle->SetPalette(1);

  const int nrRuns = 8;
  int runID[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};
  double runIDCopy[nrRuns] = {2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448};

/*  const int nrRuns = 67;
  int runID[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865, 2866, 2867, 2883, 2884, 2885, 2886, 2895, 2897, 2898, 2904, 2923, 2927, 2933, 2938, 2943, 2947, 2955, 2966, 2971, 2981, 2986, 2991, 2995, 3007, 3018, 3019, 3024, 3028, 3032, 3034, 3048, 3053, 3057, 3074};
  double runIDCopy[nrRuns] = {2417, 2418, 2421, 2422, 2423, 2424, 2426, 2431, 2432, 2433, 2434, 2447, 2448, 2714, 2719, 2725, 2732, 2737, 2743, 2748, 2761, 2766, 2771, 2785, 2804, 2811, 2817, 2828, 2837, 2843, 2848, 2858, 2865, 2866, 2867, 2883, 2884, 2885, 2886, 2895, 2897, 2898, 2904, 2923, 2927, 2933, 2938, 2943, 2947, 2955, 2966, 2971, 2981, 2986, 2991, 2995, 3007, 3018, 3019, 3024, 3028, 3032, 3034, 3048, 3053, 3057, 3074};*/

  PeakPos1 = new double[nrRuns];
  PeakPosErr1 = new double[nrRuns];

  PeakPos2 = new double[nrRuns];
  PeakPosErr2 = new double[nrRuns];

  for (int i = 0; i < nrRuns; i++) {
     ProcessRun(runID[i],i);
  }

  cout << "Mean (single site) = " << TMath::Mean(nrRuns,PeakPos1) << " +- " << TMath::RMS(nrRuns,PeakPos1) << endl;
  cout << "Mean (multi site) = " << TMath::Mean(nrRuns,PeakPos2) << " +- " << TMath::RMS(nrRuns,PeakPos2) << endl;
  
  TGraphErrors *gr1 = new TGraphErrors(nrRuns,runIDCopy,PeakPos1,0,PeakPosErr1);
  TGraphErrors *gr2 = new TGraphErrors(nrRuns,runIDCopy,PeakPos2,0,PeakPosErr2);
  
  gr1->SetTitle("Th228 peak (charge)");
  gr2->SetTitle("Th228 peak (scintillation)");
  
  gr1->SetName("Th228_charge");
  gr2->SetName("Th228_scint");
  
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);
  
  gr1->GetYaxis()->SetRangeUser(2200,3000);
  gr1->GetXaxis()->SetLimits(2400,2900);
  gr1->SetTitle("Th228 Peak (charge)");
  gr1->GetXaxis()->SetTitle("run number");
  gr1->GetYaxis()->SetTitle("peak position");
  
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);
  
  gr2->GetYaxis()->SetRangeUser(7500,9000);
  gr2->GetXaxis()->SetLimits(2400,2900);
  gr2->SetTitle("Th228 Peak (scintillation)");
  gr2->GetXaxis()->SetTitle("run number");
  gr2->GetYaxis()->SetTitle("peak position");
  
  TCanvas *c1 = new TCanvas();
  gr1->Draw("AP");
  
  TCanvas *c2 = new TCanvas();
  gr2->Draw("AP");
  
  /*TFile *fOut = new TFile("../analysis/Th228_2D_PeakPosition.root","RECREATE");
  
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

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,3500,200,0,15000);

  t->Draw("csc:epcrec>>h2DSingle","nsite == 1");

  if (h2DSingle->GetEntries() < 18000) {return;}
  
  TF2 *f2 = new TF2("f2",fun2,2450,2900,7500,9800,6);
  f2->SetParameters(40,2640,8500,70,550,-0.13);
  f2->SetParNames("A","E_I","E_S","#sigma_I","#sigma_S","#Theta");

  f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,2560,2800);
  f2->SetParLimits(2,8000,9000);
  f2->SetParLimits(3,10,200);
  f2->SetParLimits(4,20,650);
  f2->SetParLimits(5,-0.2,-0.05);

  h2DSingle->Fit("f2","r");

  double par[6];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();

  /*char FormatedOutput[1000];
  sprintf(FormatedOutput,"%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\n",RunID,0,par[1],par[2],par[3],par[4],par[5],parErr[1],parErr[2],parErr[3],parErr[4],parErr[5]);
  
  FILE *ResultFile;
  ResultFile = fopen("Th228.dat","a");
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
