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

void Cs137RotatedSpectrum()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/V11/2450_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2469_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2473_noReclustering_Fiducial.root");
  
  /*double theta_ss = 0.1814;
  double theta_ms = 0.2036;
  
  double p0_ss = -1.221;
  double p1_ss = 0.6188;
  
  double p0_ms = -0.4439;
  double p1_ms = 0.5958;
  
  double p0_1u_ss = 0.8928;
  double p1_1u_ss = 0.6159;
  
  double p0_2u_ss = p0_1u_ss;
  double p1_2u_ss = 1.00743*p1_1u_ss;*/
  
  //double p0_2u_ss = -1.136;
  //double p1_2u_ss = 0.6204;
  
  // #### Calibration 03/22/12 ###############################################################
  
  double theta_ss = 0.1814;
  double theta_ms = 0.2036;
  
  double p0_ss = 1.305;
  double p1_ss = 0.6174;
  
  double p0_ms = 2.259;
  double p1_ms = 0.5943;
  
  double p0_1u_ss = 3.352;
  double p1_1u_ss = 0.6147;
  
  //double p0_2u_ss = p0_1u_ss;
  //double p1_2u_ss = 1.00683*p1_1u_ss;
  
  double p0_2u_ss = 2.191;
  double p1_2u_ss = 0.6189;
  
  // #########################################################################################
  
  TH1F *h_ss = new TH1F("h_ss",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,0,1500);
  char *cmd_ss = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss",theta_ss,theta_ss,p1_ss,p0_ss);
  
  TH1F *h_ms = new TH1F("h_ms",Form("Th228 (multi site, #theta = %.4f)",theta_ms),200,0,1500);
  char *cmd_ms = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ms",theta_ms,theta_ms,p1_ms,p0_ms);
  
  TH1F *h_ss_1 = new TH1F("h_ss_1",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,0,1500);
  char *cmd_ss_1 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_1",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ss_2 = new TH1F("h_ss_2",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,0,1500);
  char *cmd_ss_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_2",theta_ss,theta_ss,p1_2u_ss,p0_2u_ss);
  
  TH1F *h_ss_3 = new TH1F("h_ss_3",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,0,1500);
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
  
  double fitRange1_ss = 530;
  double fitRange2_ss = 750;
  double fitRange1_ms = 530;
  double fitRange2_ms = 750;
  double Peak = 662;
  
  TF1 *fit_ss = new TF1("fit_ss",fitFunction,fitRange1_ss,fitRange2_ss,4);
  
  fit_ss->SetParNames("A1","E1","#sigma","A2");
  
  fit_ss->SetParameters(10000,Peak,100,0.2);
  
  fit_ss->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
  fit_ss->SetParLimits(2,20,200);
  fit_ss->SetParLimits(3,0.05,2.0);
  
  fit_ss->SetLineColor(kBlue);
  fit_ss->SetLineWidth(2);
  
  TF1 *fit_ms = new TF1("fit_ms",fitFunction,fitRange1_ms,fitRange2_ms,4);
  
  fit_ms->SetParNames("A1","E1","#sigma","A2");
  
  fit_ms->SetParameters(10000,Peak,100,0.2);
  
  fit_ms->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
  fit_ms->SetParLimits(2,20,200);
  fit_ms->SetParLimits(3,0.05,2.0);
  
  fit_ms->SetLineColor(kBlue);
  fit_ms->SetLineWidth(2);
  
  h_ss->Fit("fit_ss","r");
  h_ms->Fit("fit_ms","r");
  
  double E_ss = fit_ss->GetParameter(1);
  double sigma_ss = fit_ss->GetParameter(2);
  double E_err_ss = fit_ss->GetParError(1);
  double sigma_err_ss = fit_ss->GetParError(2);
  
  double E_ms = fit_ms->GetParameter(1);
  double sigma_ms = fit_ms->GetParameter(2);
  double E_err_ms = fit_ms->GetParError(1);
  double sigma_err_ms = fit_ms->GetParError(2);
  
  cout << "Resolution from best theta: " << sigma_ss / E_ss << " +- " << TMath::Sqrt(sigma_err_ss**2 / E_ss**2 + (sigma_ss/E_ss**2)**2 * E_err_ss**2)  << endl;
  
  cout << "Resolution from best theta: " << sigma_ms / E_ms << " +- " << TMath::Sqrt(sigma_err_ms**2 / E_ms**2 + (sigma_ms/E_ms**2)**2 * E_err_ms**2)  << endl;
  
  TCanvas *c1 = new TCanvas();
  h_ss->Draw("EZP");
  //hSum->Draw("EZP");
  h_ss_1->Draw("EZPsame");
  h_ss_2->Draw("EZPsame");
  
  TCanvas *c2 = new TCanvas();
  h_ms->Draw("EZP");
}