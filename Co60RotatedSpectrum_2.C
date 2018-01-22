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
  
  return gauss1 + erf1 + gauss2 + erf2 + par[7]*TMath::Exp(-x[0]/par[8]);
}

void Co60RotatedSpectrum()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/V10/2526_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2538_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2543_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2566_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2578_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2596_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2608_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2620_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2634_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2635_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2640_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2646_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2653_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2667_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2683_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2689_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V10/2708_noReclustering_Fiducial.root");
  
  double theta_ss = 0.1814;
  double theta_ms = 0.2036;
  
  double p0_ss = -1.221;
  double p1_ss = 0.6188;
  
  double p0_ms = -0.4439;
  double p1_ms = 0.5958;
  
  double p0_1u_ss = 0.8928;
  double p1_1u_ss = 0.6159;
  
  //double p0_2u_ss = p0_1u_ss;
  //double p1_2u_ss = 1.00743*p1_1u_ss;
  
  double p0_2u_ss = -1.136;
  double p1_2u_ss = 0.6204;
  
  TH1F *h_ss = new TH1F("h_ss",Form("Co60 (single site, #theta = %.4f)",theta_ss),200,500,2000);
  char *cmd_ss = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss",theta_ss,theta_ss,p1_ss,p0_ss);
  
  TH1F *h_ms = new TH1F("h_ms",Form("Co60 (multi site, #theta = %.4f)",theta_ms),200,500,2000);
  char *cmd_ms = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ms",theta_ms,theta_ms,p1_ms,p0_ms);
  
  TH1F *h_ss_1 = new TH1F("h_ss_1",Form("Co60 (single site, #theta = %.4f)",theta_ss),200,500,2000);
  char *cmd_ss_1 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_1",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ss_2 = new TH1F("h_ss_2",Form("Co60 (single site, #theta = %.4f)",theta_ss),200,500,2000);
  char *cmd_ss_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_2",theta_ss,theta_ss,p1_2u_ss,p0_2u_ss);
  
  TH1F *h_ss_3 = new TH1F("h_ss_3",Form("Co60 (single site, #theta = %.4f)",theta_ss),200,500,2000);
  char *cmd_ss_3 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_3",theta_ss,theta_ss,p1_ss,p0_ss);
  
  h_ss_1->SetLineColor(kRed);
  h_ss_2->SetLineColor(kGreen);
  
  t->Draw(cmd_ss,"nsite == 1 && nWires <= 2");
  t->Draw(cmd_ms,"nsite > 1 || nsite == 1 && nWires > 2");
  t->Draw(cmd_ss_1,"nsite == 1 && nWires == 1");
  t->Draw(cmd_ss_2,"nsite == 1 && nWires == 2");
  
  TH1F *hSum = h_ss_1->Clone("hSum");
  hSum->Add(h_ss_2);
  
  hSum->SetLineColor(kBlack);
  
  double fitRange1 = 500;
  double fitRange2 = 1500;
  double Peak1 = 1173;
  double Peak2 = 1333;
  
  TF1 *fit_ss = new TF1("fit_ss",fitFunction,fitRange1,fitRange2,9);
  
  fit_ss->SetParNames("A1","E1","#sigma","R1","R2","E2","#sigma2","p1","p0");
  
  fit_ss->SetParameters(1000,Peak1,75,0.6,0.8,Peak2,75,8500,200);
  
  fit_ss->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
  fit_ss->SetParLimits(2,20,200);
  fit_ss->SetParLimits(3,0.05,1.0);
  fit_ss->SetParLimits(4,0.05,2.0);
  fit_ss->SetParLimits(5,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
  fit_ss->SetParLimits(6,20,200);
  
  fit_ss->SetLineColor(kBlue);
  fit_ss->SetLineWidth(2);
  
  TF1 *fit_ms = new TF1("fit_ms",fitFunction,fitRange1,fitRange2,9);
  
  fit_ms->SetParNames("A1","E1","#sigma","R1","R2","E2","#sigma2","p1","p0");
  
  fit_ms->SetParameters(1000,Peak1,75,0.6,0.8,Peak2,75,8500,200);
  
  fit_ms->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
  fit_ms->SetParLimits(2,20,200);
  fit_ms->SetParLimits(3,0.05,1.0);
  fit_ms->SetParLimits(4,0.05,2.0);
  fit_ms->SetParLimits(5,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
  fit_ms->SetParLimits(6,20,200);
  
  fit_ms->SetLineColor(kBlue);
  fit_ms->SetLineWidth(2);
  
  h_ss->Fit("fit_ss","r");
  h_ms->Fit("fit_ms","r");
  
  double E1 = fit_ss->GetParameter(1);
  double sigma1 = fit_ss->GetParameter(2);
  double E1_err = fit_ss->GetParError(1);
  double sigma1_err = fit_ss->GetParError(2);
  
  double E2 = fit_ss->GetParameter(5);
  double sigma2 = fit_ss->GetParameter(6);
  double E2_err = fit_ss->GetParError(5);
  double sigma2_err = fit_ss->GetParError(6);
  
  cout << "Resolution from best theta: " << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2)  << endl;
  
  cout << "Resolution from best theta: " << sigma2 / E2 << " +- " << TMath::Sqrt(sigma2_err**2 / E2**2 + (sigma2/E2**2)**2 * E2_err**2)  << endl;
  
  TF1 *fit_lin_ss = new TF1("fit_lin_ss","[1]*TMath::Exp(-x/[0])",500,2000);
  fit_lin_ss->SetParameters(fit_ss->GetParameter(8),fit_ss->GetParameter(7));
  fit_lin_ss->SetLineColor(kRed);
  fit_lin_ss->SetLineWidth(2);
  
  TF1 *fit_lin_ms = new TF1("fit_lin_ms","[1]*TMath::Exp(-x/[0])",500,2000);
  fit_lin_ms->SetParameters(fit_ms->GetParameter(8),fit_ms->GetParameter(7));
  fit_lin_ms->SetLineColor(kRed);
  fit_lin_ms->SetLineWidth(2);
  
  TCanvas *c1 = new TCanvas();
  h_ss->Draw("EZP");
  //hSum->Draw("EZP");
  h_ss_1->Draw("EZPsame");
  h_ss_2->Draw("EZPsame");
  fit_lin_ss->Draw("same");
  
  TCanvas *c2 = new TCanvas();
  h_ms->Draw("EZP");
  fit_lin_ms->Draw("same");
}