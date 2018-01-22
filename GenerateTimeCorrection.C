double theta_ss = 0.1814;
double theta_ms = 0.2036;
double p0_ss = 1.305;
double p1_ss = 0.6174;
double p0_ms = 2.259;
double p1_ms = 0.5943;

void SetGraphStyle(TGraphErrors *gr, char *title, char *xAxis, char *yAxis, int MarkerSize, int MarkerStyle, int MarkerColor);

TCanvas *c1 = new TCanvas();
TCanvas *c2 = new TCanvas();

double fitFunctionTh228(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  
  return gauss1 + erf1;
}

double fitFunctionCo60(double *x, double *par)
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

void GenerateTimeCorrection()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  const int nrRunsCo60_S2 = 6;
  const int nrRunsCo60_S5 = 6;
  const int nrRunsCo60_S8 = 5;
  const int nrRunsCo60_S11 = 0;
  const int nrRunsCo60_S17 = 0;
  const int nrRunsTh228_S2 = 16;
  const int nrRunsTh228_S5 = 79;
  const int nrRunsTh228_S8 = 17;
  const int nrRunsTh228_S11 = 2;
  const int nrRunsTh228_S17 = 2;
  const int nrRunsTh228_S5_1 = 4;
  const int nrRunsTh228_S5_2 = 73;
  
  int runIDCo60_S2[nrRunsCo60_S2] = {2538,2566,2608,2640,2667,2689};
  int runIDCo60_S5[nrRunsCo60_S5] = {2543,2578,2620,2635,2646,2708};
  int runIDCo60_S8[nrRunsCo60_S8] = {2526,2596,2634,2653,2683};
  int runIDCo60_S11[nrRunsCo60_S11] = {};
  int runIDCo60_S17[nrRunsCo60_S17] = {};
  int runIDTh228_S2[nrRunsTh228_S2] = {2422,2432,2433,2714,2732,2748,2771,2811,2843,2923,2943,2971,2995,3019,3048,3091};
  int runIDTh228_S5[nrRunsTh228_S5] = {2418,2423,2424,2426,2434,2447,2719,2737,2761,2785,2817,2837,2848,2883,2904,2938,2947,2966,2981,2991,3007,3024,3028,3034,3053,3074,3099,3109,3120,3137,3142,3147,3152,3153,3158,3169,3173,3178,3184,3195,3199,3204,3214,3218,3223,3231,3236,3246,3254,3258,3262,3266,3270,3281,3283,3295,3298,3342,3345,3347,3351,3352,3355,3358,3361,3364,3367,3376,3379,3382,3383,3386,3389,3392,3395,3404,3407,3410,3413};
  int runIDTh228_S8[nrRunsTh228_S8] = {2421,2431,2725,2743,2766,2804,2828,2858,2897,2898,2933,2955,2986,3018,3057,3222,3308};
  int runIDTh228_S11[nrRunsTh228_S11] = {2448,2865};
  int runIDTh228_S17[nrRunsTh228_S17] = {2866,2867};
  int runIDTh228_S5_1[nrRunsTh228_S5_1] = {2424,2426,2434,2447};
  int runIDTh228_S5_2[nrRunsTh228_S5_2] = {2719,2737,2761,2785,2817,2837,2848,2883,2904,2938,2947,2966,2981,2991,3007,3024,3028,3034,3053,3074,3099,3109,3120,3137,3142,3147,3152,3153,3158,3169,3173,3178,3184,3195,3199,3204,3214,3218,3223,3231,3236,3246,3254,3258,3262,3266,3270,3281,3283,3295,3298,3342,3345,3347,3351,3352,3355,3358,3361,3364,3367,3376,3379,3382,3383,3386,3389,3392,3395,3404,3407,3410,3413};
  
  // 3310, 3311, 3313
  
  TGraphErrors *grCo60_S2 = new TGraphErrors(nrRunsCo60_S2);
  TGraphErrors *grCo60_S5 = new TGraphErrors(nrRunsCo60_S5);
  TGraphErrors *grCo60_S8 = new TGraphErrors(nrRunsCo60_S8);
  TGraphErrors *grCo60_S11 = new TGraphErrors(nrRunsCo60_S11);
  TGraphErrors *grCo60_S17 = new TGraphErrors(nrRunsCo60_S17);
  TGraphErrors *grTh228_S2 = new TGraphErrors(nrRunsTh228_S2);
  TGraphErrors *grTh228_S5 = new TGraphErrors(nrRunsTh228_S5);
  TGraphErrors *grTh228_S8 = new TGraphErrors(nrRunsTh228_S8);
  TGraphErrors *grTh228_S11 = new TGraphErrors(nrRunsTh228_S11);
  TGraphErrors *grTh228_S17 = new TGraphErrors(nrRunsTh228_S17);
  TGraphErrors *grCorrection = new TGraphErrors(nrRunsTh228_S5_1 + nrRunsTh228_S5_2 + nrRunsCo60_S5);
  
  SetGraphStyle(grCo60_S2,"Peak positions Co60 S2","unix time","peak position",0.8,24,kBlack);
  SetGraphStyle(grCo60_S5,"Peak positions Co60 S5","unix time","peak position",0.8,25,kBlack);
  SetGraphStyle(grCo60_S8,"Peak positions Co60 S8","unix time","peak position",0.8,26,kBlack);
  SetGraphStyle(grCo60_S11,"Peak positions Co60 S11","unix time","peak position",0.8,27,kBlack);
  SetGraphStyle(grCo60_S17,"Peak positions Co60 S17","unix time","peak position",0.8,28,kBlack);
  SetGraphStyle(grTh228_S2,"Peak positions Th228 S2","unix time","peak position",0.8,20,kBlack);
  SetGraphStyle(grTh228_S5,"Peak positions Th228 S5","unix time","peak position",0.8,21,kBlack);
  SetGraphStyle(grTh228_S8,"Peak positions Th228 S8","unix time","peak position",0.8,22,kBlack);
  SetGraphStyle(grTh228_S11,"Peak positions Th228 S11","unix time","peak position",0.8,23,kBlack);
  SetGraphStyle(grTh228_S17,"Peak positions Th228 S17","unix time","peak position",0.8,34,kBlack);
  SetGraphStyle(grCorrection,"Peak positions S5","unix time","peak position",0.8,20,kBlack);
  
  grCorrection->SetName("ThCoS5PeakPosTimeCor");
  
  ProcessCo60Run(grCo60_S2,nrRunsCo60_S2,runIDCo60_S2);
  ProcessCo60Run(grCo60_S5,nrRunsCo60_S5,runIDCo60_S5);
  ProcessCo60Run(grCo60_S8,nrRunsCo60_S8,runIDCo60_S8);
  ProcessCo60Run(grCo60_S11,nrRunsCo60_S11,runIDCo60_S11);
  ProcessCo60Run(grCo60_S17,nrRunsCo60_S17,runIDCo60_S17);
  
  ProcessTh228Run(grTh228_S2,nrRunsTh228_S2,runIDTh228_S2);
  ProcessTh228Run(grTh228_S5,nrRunsTh228_S5,runIDTh228_S5);
  ProcessTh228Run(grTh228_S8,nrRunsTh228_S8,runIDTh228_S8);
  ProcessTh228Run(grTh228_S11,nrRunsTh228_S11,runIDTh228_S11);
  ProcessTh228Run(grTh228_S17,nrRunsTh228_S17,runIDTh228_S17);
  
  /*c1->Print("TimeCorrectionFits.ps[");
  ProcessTh228Run(grCorrection,nrRunsTh228_S5_1,runIDTh228_S5_1);
  ProcessCo60Run(grCorrection,nrRunsCo60_S5,runIDCo60_S5,nrRunsTh228_S5_1);
  ProcessTh228Run(grCorrection,nrRunsTh228_S5_2,runIDTh228_S5_2,nrRunsTh228_S5_1 + nrRunsCo60_S5);
  c1->Print("TimeCorrectionFits.ps]");*/
  
  //cout << grCorrection->GetRMS(2) << endl;
  
  TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
  l->AddEntry(grCo60_S2,"Co60 S2","lp");
  l->AddEntry(grCo60_S5,"Co60 S5","lp");
  l->AddEntry(grCo60_S8,"Co60 S8","lp");
  l->AddEntry(grTh228_S2,"Th228 S2","lp");
  l->AddEntry(grTh228_S5,"Th228 S5","lp");
  l->AddEntry(grTh228_S8,"Th228 S8","lp");
  l->AddEntry(grTh228_S11,"Th228 S11","lp");
  l->AddEntry(grTh228_S17,"Th228 S17","lp");
  
  grCo60_S2->GetXaxis()->SetLimits(1315000000,1333000000);
  
  TCanvas *c1 = new TCanvas("c1","Th228 S5");
  grCo60_S2->Draw("AZP");
  grCo60_S5->Draw("ZPsame");
  grCo60_S8->Draw("ZPsame");
  grTh228_S2->Draw("ZPsame");
  grTh228_S5->Draw("ZPsame");
  grTh228_S8->Draw("ZPsame");
  grTh228_S11->Draw("PZsame");
  grTh228_S17->Draw("ZPsame");
  l->Draw("same");
  
  /*TCanvas *c2 = new TCanvas("c2","Correction");
  grCorrection->Draw("AZPL");
  
  TFile *f = new TFile("ThCoS5PeakPosTimeCor.root","RECREATE");
  grCorrection->Write();
  f->Close();*/
  
  return;
}

void ProcessCo60Run(TGraphErrors *gr, int nrRuns, int *runID, int offset = 0)
{
  TFile *f = 0;
  TTree *t = 0;
  
  TH1F *h_ss = 0;
  
  for (int i = 0; i < nrRuns; i++) {
    
    char FileName[100];
    //sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",runID[i]);
    sprintf(FileName,"../../EnergyCalibration/analysis/V11/%i_noReclustering_Fiducial.root",runID[i]);
  
    f = new TFile(FileName,"READ");
    t = (TTree*)f->Get("t");
  
    cout << "Processing run " << runID[i] << endl;
    cout << "    fitting .... ";
  
    // Reprocessed data (01/30/12)
    /*double theta_ss = 0.1916;
    double theta_ms = 0.2164;
    double p0_ss = 2.443;
    double p1_ss = 0.6117;
    double p0_ms = 8.801;
    double p1_ms = 0.5826;*/
  
    // Reprocessed data (02/07/12)
    /*double theta_ss = 0.1867;
    double theta_ms = 0.2084;
    double p0_ss = 1.196;
    double p1_ss = 0.6152;
    double p0_ms = -0.1598;
    double p1_ms = 0.5936;*/
    
    // Reprocessed data (02/18/12)
    /*double theta_ss = 0.1867;
    double theta_ms = 0.2084;
    double p0_ss = 2.057;
    double p1_ss = 0.6133;
    double p0_ms = 7.283;
    double p1_ms = 0.5874;*/
    
    // Reprocessed data (02/22/12)
    /*double theta_ss = 0.1860;
    double theta_ms = 0.2064;
    double p0_ss = 3.81;
    double p1_ss = 0.6114;
    double p0_ms = 7.6951;
    double p1_ms = 5.9;*/
  
    TH1F *h_ss = new TH1F("h_ss",Form("Co60 (single site, #theta = %.4f)",theta_ss),100,500,2000);
    //TH1F *h_ms = new TH1F("h_ms",Form("Co60 (multi site, #theta = %.4f)",theta_ms),100,500,2000);
    
    char *cmd_ss = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.6f + %.6f >> h_ss",theta_ss,theta_ss,p1_ss,p0_ss);
    //char *cmd_ms = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.6f + %.6f >> h_ms",theta_ms,theta_ms,p1_ms,p0_ms);
  
    t->Draw(cmd_ss,"nsite == 1 && nWires <= 2","goff");
    //t->Draw(cmd_ms,"nsite > 1 || nsite == 1 && nWires > 2","goff");
    
    double fitRange1 = 1050;
    double fitRange2 = 1430;
    double Peak1 = 1173;
    double Peak2 = 1332;
  
    TF1 *fitMin = new TF1("fitMin",fitFunctionCo60,fitRange1,fitRange2,7);
    fitMin->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  
    fitMin->SetParameters(3000,Peak1,75,0.6,0.8,Peak2,75);
    fitMin->SetParLimits(1,Peak1-Peak1/20.0,Peak1+Peak1/20.0);
    fitMin->SetParLimits(2,10,120);
    fitMin->SetParLimits(3,0.05,1.0);
    fitMin->SetParLimits(4,0.5,1.5);
    fitMin->SetParLimits(5,Peak2-Peak2/20.0,Peak2+Peak2/20.0);
    fitMin->SetParLimits(6,10,130);
  
    h_ss->Fit("fitMin","rq");
  
    double E1 = fitMin->GetParameter(1);
    double sigma1 = fitMin->GetParameter(2);
    double E2 = fitMin->GetParameter(5);
    double sigma2 = fitMin->GetParameter(6);
    double E1_err = fitMin->GetParError(1);
    double sigma1_err = fitMin->GetParError(2);
    double E2_err = fitMin->GetParError(5);
    double sigma2_err = fitMin->GetParError(6);
    
    double time_err = 0;
    double time = GetTime(&time_err,runID[i]);
    cout << E2 << " +- " << E2_err << "    " << time << " +- " << time_err << endl;
    
    gr->SetPoint(i+offset,time,E2/1332.492*2614.511);  // <-- Get precise values
    gr->SetPointError(i+offset,time_err,E2_err/1332.492*2614.511);  // <-- Get precise values
    
    //c1->cd();
    //h_ss->Draw("EZP");
    //c1->Print("TimeCorrectionFits.ps");
    
    /*c1->cd();
    h_ss->Draw("EZP");
    char oName_ss[100];
    sprintf(oName_ss,"/var/www/html/analysis/CorrectionFunction/SingleSite/%i.png",runID[i]);
    c1->SaveAs(oName_ss);
    
    c2->cd();
    h_ms->Draw("EZP");
    char oName_ms[100];
    sprintf(oName_ms,"/var/www/html/analysis/CorrectionFunction/MultiSite/%i.png",runID[i]);
    c2->SaveAs(oName_ms);*/
    
    //delete f;
    //delete t;
    //delete h_ss;
    //delete fitMin;
    
    f->Close();
  }
  
  //delete f;
  //delete t;
  //delete h_ss;
  
  return;
}

void ProcessTh228Run(TGraphErrors *gr, int nrRuns, int *runID, int offset = 0)
{
  TFile *f = 0;
  TTree *t = 0;
  
  TH1F *h_ss = 0;
  
  for (int i = 0; i < nrRuns; i++) {
    
    char FileName[100];
    //sprintf(FileName,"../../EnergyCalibration/analysis/V9/%i_noReclustering_Fiducial.root",runID[i]);
    sprintf(FileName,"../../EnergyCalibration/analysis/V11/%i_noReclustering_Fiducial.root",runID[i]);
  
    f = new TFile(FileName,"READ");
    t = (TTree*)f->Get("t");
  
    cout << "Processing run " << runID[i] << endl;
    cout << "    fitting .... ";
    
    // Reprocessed data (01/30/12)
    /*double theta_ss = 0.1916;
    double theta_ms = 0.2164;
    double p0_ss = 2.443;
    double p1_ss = 0.6117;
    double p0_ms = 8.801;
    double p1_ms = 0.5826;*/
    
    // Reprocessed data (02/07/12)
    /*double theta_ss = 0.1867;
    double theta_ms = 0.2084;
    double p0_ss = 1.196;
    double p1_ss = 0.6152;
    double p0_ms = -0.1598;
    double p1_ms = 0.5936;*/
    
    // Reprocessed data (02/18/12)
    /*double theta_ss = 0.1867;
    double theta_ms = 0.2084;
    double p0_ss = 2.057;
    double p1_ss = 0.6133;
    double p0_ms = 7.283;
    double p1_ms = 0.5874;*/
    
    // Reprocessed data (02/22/12)
    /*double theta_ss = 0.1860;
    double theta_ms = 0.2064;
    double p0_ss = 3.81;
    double p1_ss = 0.6114;
    double p0_ms = 7.6951;
    double p1_ms = 5.9;*/
    
    h_ss = new TH1F("h_ss",Form("Th228 (single site, #theta = %.4f)",theta_ss),100,2000,3500);
    //TH1F *h_ms = new TH1F("h_ms",Form("Th228 (multi site, #theta = %.4f)",theta_ms),100,2000,3500);
    
    char *cmd_ss = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.6f + %.6f >> h_ss",theta_ss,theta_ss,p1_ss,p0_ss);
    //char *cmd_ms = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.6f + %.6f >> h_ms",theta_ms,theta_ms,p1_ms,p0_ms);
    
    t->Draw(cmd_ss,"nsite == 1 && nWires <= 2","goff");
    //t->Draw(cmd_ms,"nsite > 1 || nsite == 1 && nWires > 2","goff");
    
    double fitRange1 = 2480;
    double fitRange2 = 2720;
    double Peak = 2615;
    
    TF1 *fitMin = new TF1("fitMin",fitFunctionTh228,fitRange1,fitRange2,4);
  
    fitMin->SetParNames("A1","E1","#sigma","A2");
  
    fitMin->SetParameters(10000,Peak,50,0.2);
  
    fitMin->SetParLimits(1,Peak-Peak/30.0,Peak+Peak/30.0);
    fitMin->SetParLimits(2,30,200);
    fitMin->SetParLimits(3,0.01,2.0);
  
    h_ss->Fit("fitMin","rq");
  
    double E1 = fitMin->GetParameter(1);
    double sigma1 = fitMin->GetParameter(2);
    double E1_err = fitMin->GetParError(1);
    double sigma1_err = fitMin->GetParError(2);
    
    double time_err = 0;
    double time = GetTime(&time_err,runID[i]);
    cout << E1 << " +- " << E1_err << "    " << time << " +- " << time_err << endl;
    
    gr->SetPoint(i+offset,time,E1);  // <-- Get precise values
    gr->SetPointError(i+offset,time_err,E1_err);  // <-- Get precise values
    
    //c1->cd();
    //h_ss->Draw("EZP");
    //c1->Print("TimeCorrectionFits.ps");
    
    /*c1->cd();
    h_ss->Draw("EZP");
    char oName_ss[100];
    sprintf(oName_ss,"/var/www/html/analysis/CorrectionFunction/SingleSite/%i.png",runID[i]);
    c1->SaveAs(oName_ss);
    
    c2->cd();
    h_ms->Draw("EZP");
    char oName_ms[100];
    sprintf(oName_ms,"/var/www/html/analysis/CorrectionFunction/MultiSite/%i.png",runID[i]);
    c2->SaveAs(oName_ms);*/
    
    //delete f;
    //delete t;
    //delete h_ss;
    //delete fitMin;
    
    f->Close();
  }
  
  //delete f;
  //delete t;
  //delete h_ss;
  
  return;
}

double GetTime(double *dt, int runID)
{
  char fname[100];
  sprintf(fname,"/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/%i/recon0000%i-*.root",runID,runID);
  
  TChain *t = new TChain("tree");
  
  t->Add(fname);
  
  int nentries = t->GetEntries();
  EXOEventData *ED = 0;
  
  t->SetBranchAddress("EventBranch",&ED);
  
  t->GetEntry(0);
  double tStart = ED->fEventHeader.fTriggerSeconds;
  t->GetEntry(nentries - 1);
  double tEnd = ED->fEventHeader.fTriggerSeconds;
  
  *dt = (tEnd - tStart) / 2.0;
  
  delete t;
  
  return tStart + (tEnd - tStart) / 2.0;
}

void SetGraphStyle(TGraphErrors *gr, char *title, char *xAxis, char *yAxis, double MarkerSize, int MarkerStyle, int MarkerColor)
{
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle(xAxis);
  gr->GetYaxis()->SetTitle(yAxis);
  gr->SetMarkerSize(MarkerSize);
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  
  return;
}
