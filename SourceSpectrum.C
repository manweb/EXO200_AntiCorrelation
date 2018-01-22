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

double fitFunctionCs137(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  
  return gauss1 + erf1;
}

void SourceSpectrum()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  TChain *tTh228 = new TChain("t","tree");
  TChain *tCo60 = new TChain("t","tree");
  TChain *tCs137 = new TChain("t","tree");
  
  // Reprocessed data (02/07/12)
  tTh228->Add("../../EnergyCalibration/analysis/V9/2424_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2426_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2431_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2432_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2433_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2434_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2447_noReclustering_Fiducial.root");
  tTh228->Add("../../EnergyCalibration/analysis/V9/2448_noReclustering_Fiducial.root");
  
  tCo60->Add("../../EnergyCalibration/analysis/V9/2526_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2538_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2543_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2566_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2578_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2596_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2608_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2620_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2634_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2635_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2640_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2646_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2653_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2667_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2683_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2689_noReclustering_Fiducial.root");
  tCo60->Add("../../EnergyCalibration/analysis/V9/2708_noReclustering_Fiducial.root");
  
  tCs137->Add("../../EnergyCalibration/analysis/V9/2450_noReclustering_Fiducial.root");
  tCs137->Add("../../EnergyCalibration/analysis/V9/2469_noReclustering_Fiducial.root");
  tCs137->Add("../../EnergyCalibration/analysis/V9/2473_noReclustering_Fiducial.root");
  
  double theta_ss = 0.1867;
  double theta_ms = 0.2084;
  
  TH1F *hTh228_ss = new TH1F("hTh228_ss","Th228 (single site, #theta = 0.1867)",200,0,3500);
  TH1F *hTh228_ms = new TH1F("hTh228_ms","Th228 (multi site, #theta = 0.2084)",200,0,3500);
  TH1F *hCo60_ss = new TH1F("hCo60_ss","Co60 (single site, #theta = 0.1867)",200,0,2000);
  TH1F *hCo60_ms = new TH1F("hCo60_ms","Co60  (multi site, #theta = 0.2084)",200,0,2000);
  TH1F *Cs137_ss = new TH1F("hCs137_ss","Cs137 (single site, #theta = 0.1867)",200,0,1500);
  TH1F *Cs137_ms = new TH1F("hCs137_ms","Cs137 (multi site, #theta = 0.2084)",200,0,1500);
  
  hTh228_ss->GetXaxis()->SetTitle("energy [keV]");
  hTh228_ms->GetXaxis()->SetTitle("energy [keV]");
  hCo60_ss->GetXaxis()->SetTitle("energy [keV]");
  hCo60_ms->GetXaxis()->SetTitle("energy [keV]");
  hCs137_ss->GetXaxis()->SetTitle("energy [keV]");
  hCs137_ms->GetXaxis()->SetTitle("energy [keV]");
  
  char cmdTh228_ss[100];
  char cmdTh228_ms[100];
  char cmdCo60_ss[100];
  char cmdCo60_ms[100];
  char cmdCs137_ss[100];
  char cmdCs137_ms[100];
  
  sprintf(cmdTh228_ss,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196 >> hTh228_ss",theta_ss,theta_ss);
  sprintf(cmdTh228_ms,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598 >> hTh228_ms",theta_ms,theta_ms);
  sprintf(cmdCo60_ss,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196 >> hCo60_ss",theta_ss,theta_ss);
  sprintf(cmdCo60_ms,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598 >> hCo60_ms",theta_ms,theta_ms);
  sprintf(cmdCs137_ss,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196 >> hCs137_ss",theta_ss,theta_ss);
  sprintf(cmdCs137_ms,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598 >> hCs137_ms",theta_ms,theta_ms);
  
  tTh228->Draw(cmdTh228_ss, "nsite == 1");
  tTh228->Draw(cmdTh228_ms, "nsite > 1");
  tCo60->Draw(cmdCo60_ss, "nsite == 1");
  tCo60->Draw(cmdCo60_ms, "nsite > 1");
  tCs137->Draw(cmdCs137_ss, "nsite == 1");
  tCs137->Draw(cmdCs137_ms, "nsite > 1");
  
  TF1 *fitTh228_ss = new TF1("fitTh228_ss",fitFunctionTh228,2400,2800,4);
  fitTh228_ss->SetParNames("A1","E1","#sigma","A2");
  fitTh228_ss->SetParameters(10000,2614,100,0.2);
  fitTh228_ss->SetParLimits(1,2500,2700);
  fitTh228_ss->SetParLimits(2,20,200);
  fitTh228_ss->SetParLimits(3,0.05,2.0);
  
  TF1 *fitTh228_ms = new TF1("fitTh228_ms",fitFunctionTh228,2300,2750,4);
  fitTh228_ms->SetParNames("A1","E1","#sigma","A2");
  fitTh228_ms->SetParameters(10000,2614,100,0.2);
  fitTh228_ms->SetParLimits(1,2500,2700);
  fitTh228_ms->SetParLimits(2,20,200);
  fitTh228_ms->SetParLimits(3,0.05,2.0);
  
  TF1 *fitCo60_ss = new TF1("fitCo60_ss",fitFunctionCo60,1000,1500,7);
  fitCo60_ss->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  fitCo60_ss->SetParameters(3000,1173,75,0.6,0.8,1332,75);
  fitCo60_ss->SetParLimits(1,1050,1250);
  fitCo60_ss->SetParLimits(2,10,120);
  fitCo60_ss->SetParLimits(3,0.05,1.0);
  fitCo60_ss->SetParLimits(4,0.5,1.5);
  fitCo60_ss->SetParLimits(5,1250,1380);
  fitCo60_ss->SetParLimits(6,10,130);
  
  TF1 *fitCo60_ms = new TF1("fitCo60_ms",fitFunctionCo60,1000,1450,7);
  fitCo60_ms->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");
  fitCo60_ms->SetParameters(3000,1173,75,0.6,0.8,1332,75);
  fitCo60_ms->SetParLimits(1,1050,1250);
  fitCo60_ms->SetParLimits(2,10,120);
  fitCo60_ms->SetParLimits(3,0.05,1.0);
  fitCo60_ms->SetParLimits(4,0.5,1.5);
  fitCo60_ms->SetParLimits(5,1250,1380);
  fitCo60_ms->SetParLimits(6,10,130);
  
  TF1 *fitCs137_ss = new TF1("fitCs137_ss",fitFunctionCs137,500,800,4);
  fitCs137_ss->SetParNames("A1","E1","#sigma","A2");
  fitCs137_ss->SetParameters(10000,662,100,0.2);
  fitCs137_ss->SetParLimits(1,600,700);
  fitCs137_ss->SetParLimits(2,10,200);
  fitCs137_ss->SetParLimits(3,0.05,2.0);
  
  TF1 *fitCs137_ms = new TF1("fitCs137_ms",fitFunctionCs137,500,800,4);
  fitCs137_ms->SetParNames("A1","E1","#sigma","A2");
  fitCs137_ms->SetParameters(10000,662,100,0.2);
  fitCs137_ms->SetParLimits(1,600,700);
  fitCs137_ms->SetParLimits(2,10,200);
  fitCs137_ms->SetParLimits(3,0.05,2.0);
  
  hTh228_ss->Fit("fitTh228_ss","r");
  hTh228_ms->Fit("fitTh228_ms","r");
  hCo60_ss->Fit("fitCo60_ss","r");
  hCo60_ms->Fit("fitCo60_ms","r");
  hCs137_ss->Fit("fitCs137_ss","r");
  hCs137_ms->Fit("fitCs137_ms","r");
  
  TCanvas *c1 = new TCanvas("c1","Th228",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  hTh228_ss->Draw();
  c1->cd(2);
  hTh228_ms->Draw();
  
  TCanvas *c2 = new TCanvas("c2","Co60",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  hCo60_ss->Draw();
  c2->cd(2);
  hCo60_ms->Draw();
  
  TCanvas *c3 = new TCanvas("c3","Cs137",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  hCs137_ss->Draw();
  c3->cd(2);
  hCs137_ms->Draw();
  
  return;
}