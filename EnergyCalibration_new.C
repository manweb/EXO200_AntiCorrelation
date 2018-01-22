void EnergyCalibration_new()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  const int nDataPoints = 4;
  
  // true source energies
  double Etrue[nDataPoints] = {661.657, 1173.228, 1332.492, 2614.511};
  
  // measured values for single and multi site
  /*double Erec_single[nDataPoints] = {1069.43, 1899.55, 2159.59, 4225.79};
  double ErecErr_single[nDataPoints] = {6.04, 7.91, 8.85, 8.24};
  
  double Erec_multi[nDataPoints] = {1107.05, 1974.8, 2243.29, 4385.87};
  double ErecErr_multi[nDataPoints] = {6.95, 8.94, 10.16, 8.69};
  
  // measured values for single and nWires == 1
  double Erec_single[nDataPoints] = {1070.66, 1904.86, 2166.25, 4242.23};
  double ErecErr_single[nDataPoints] = {6.66, 7.57, 8.71, 8.37};
  
  // measured values for single and nWires == 2
  double Erec_multi[nDataPoints] = {1067.53, 1893.59, 2152.39, 4214.68};
  double ErecErr_multi[nDataPoints] = {4.95, 8.51, 9.03, 9.00};
  
  // measured values for single and nWires == 3
  double Erec_multi[nDataPoints] = {1062.9, 1884.07, 2141.6, 4195.86};
  double ErecErr_multi[nDataPoints] = {9.96, 7.86, 8.72, 10.24};*/
  
// #### Calibration 03/22/12 ############################################################
  
  // measured values for single and multi site
  double Erec_single[nDataPoints] = {1065.65, 1899.37, 2159.38, 4227.14};
  double ErecErr_single[nDataPoints] = {3.97, 3.85, 3.47, 6.70};
  
  double Erec_multi[nDataPoints] = {1103.18, 1974.57, 2243.1, 4387.59};
  double ErecErr_multi[nDataPoints] = {4.52, 4.51, 4.51, 7.59};
  
  // measured values for single and nWires == 1
  //double Erec_single[nDataPoints] = {1066.59, 1904.68, 2166.09, 4243.38};
  //double ErecErr_single[nDataPoints] = {4.58, 4.20, 4.18, 6.77};
  
  // measured values for single and nWires == 2
  //double Erec_multi[nDataPoints] = {1063.63, 1893.43, 2152.18, 4215.97};
  //double ErecErr_multi[nDataPoints] = {2.84, 3.86, 3.12, 6.61};
  
  // measured values for single and nWires == 3
  //double Erec_multi[nDataPoints] = {1057.0, 1883.98, 2141.59, 4197.77};
  //double ErecErr_multi[nDataPoints] = {3.62, 6.39, 5.92, 6.70};
  
// ######################################################################################
  
  // create various graphs
  TGraphErrors *gr1 = new TGraphErrors(nDataPoints,Erec_single,Etrue,ErecErr_single,0); // E_true vs E_rec (single site)
  TGraphErrors *gr2 = new TGraphErrors(nDataPoints,Etrue,Erec_single,0,ErecErr_single); // E_rec vs E_true (singel site)
  TGraphErrors *gr3 = new TGraphErrors(nDataPoints,Erec_multi,Etrue,ErecErr_multi,0); // E_true vs E_rec (multi site)
  TGraphErrors *gr4 = new TGraphErrors(nDataPoints,Etrue,Erec_multi,0,ErecErr_multi); // E_rec vs E_true (multi site)
  
  SetGraphProperties(gr1,"E_{true} vs. E_{rec} (single site)","rec. energy","true energy (keV)");
  SetGraphProperties(gr2,"E_{rec} vs. E_{true} (single site)","true energy (keV)","rec. energy");
  SetGraphProperties(gr3,"E_{true} vs. E_{rec} (multi site)","rec. energy","true energy (keV)");
  SetGraphProperties(gr4,"E_{rec} vs. E_{true} (multi site)","true energy (keV)","rec. energy");
  
  double par1_ss;
  double par2_ss;
  double par1_ms;
  double par2_ms;
  double par1_err_ss;
  double par2_err_ss;
  double par1_err_ms;
  double par2_err_ms;
  double cor_ss;
  double cor_ms;
  
  TGraphErrors *residuals_ss = LinearFits(gr1,gr2,&par1_ss,&par2_ss,&par1_err_ss,&par2_err_ss,&cor_ss);
  TGraphErrors *residuals_ms = LinearFits(gr3,gr4,&par1_ms,&par2_ms,&par1_err_ms,&par2_err_ms,&cor_ms);
  
  residuals_ss->GetYaxis()->SetRangeUser(-10.0,10.0);
  residuals_ms->GetYaxis()->SetRangeUser(-10.0,10.0);
  
  SetGraphProperties(residuals_ss,"Residuals (single site)","energy (keV)","residual (%)");
  SetGraphProperties(residuals_ms,"Residuals (multi site)","energy (keV)","residual (%)");
  
  TF1 *sigma_p_ss = new TF1("sigma_p_ss",Sigma,0,4500,6);
  TF1 *sigma_n_ss = new TF1("sigma_n_ss",Sigma,0,4500,6);
  TF1 *sigma_p_ms = new TF1("sigma_p_ms",Sigma,0,4500,6);
  TF1 *sigma_n_ms = new TF1("sigma_n_ms",Sigma,0,4500,6);
  
  sigma_p_ss->SetLineWidth(1);
  sigma_n_ss->SetLineWidth(1);
  sigma_p_ms->SetLineWidth(1);
  sigma_n_ms->SetLineWidth(1);
  
  sigma_p_ss->SetParameters(0,0,par1_err_ss,par2_err_ss,cor_ss,1);
  sigma_n_ss->SetParameters(0,0,par1_err_ss,par2_err_ss,cor_ss,-1);
  sigma_p_ms->SetParameters(0,0,par1_err_ms,par2_err_ms,cor_ms,1);
  sigma_n_ms->SetParameters(0,0,par1_err_ms,par2_err_ms,cor_ms,-1);
  
  TCanvas *c1 = new TCanvas("c1","Etrue vs Erec (single site)");
  gr1->Draw("AZP");
  
  TCanvas *c2 = new TCanvas("c2","Erec vs Etrue (single site)");
  gr2->Draw("AZP");
  
  TCanvas *c3 = new TCanvas("c3","Etrue vs Erec (multi site)");
  gr3->Draw("AZP");
  
  TCanvas *c4 = new TCanvas("c4","Erec vs Etrue (multi site)");
  gr4->Draw("AZP");
  
  TCanvas *c5 = new TCanvas("c5","Residuals (single site)");
  residuals_ss->Draw("AZP");
  sigma_p_ss->Draw("same");
  sigma_n_ss->Draw("same");
  
  residuals_ss->GetXaxis()->SetLimits(0,3000);
  c5->Update();
  
  TCanvas *c6 = new TCanvas("c6","Residuals (multi site)");
  residuals_ms->Draw("AZP");
  sigma_p_ms->Draw("same");
  sigma_n_ms->Draw("same");
  
  residuals_ms->GetXaxis()->SetLimits(0,3000);
  c6->Update();
  
  /*gr3->SetMarkerStyle(21);
  
  TCanvas *c7 = new TCanvas("c7","Etrue vs Erec, 1 and 2 u wires");
  gr1->Draw("AZP");
  gr3->Draw("ZPsame");*/

}

TGraphErrors *LinearFits(TGraphErrors *gr1, TGraphErrors *gr2, double *a, double *b, double *a_err, double *b_err, double *cor)
{
  TF1 *fit1 = new TF1("fit1","[1]*x+[0]");
  fit1->SetLineWidth(1);
  
  TF1 *fit2 = new TF1("fit2","[1]*x+[0]");
  fit2->SetLineWidth(1);
  
  // first fit the graph E_rec vs E_true
  gr2->Fit("fit1");
  
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  double var1 = fitter->GetCovarianceMatrixElement(0,0);
  double var2 = fitter->GetCovarianceMatrixElement(1,1);
  double correlation = fitter->GetCovarianceMatrixElement(0,1);
  
  double par1 = fit1->GetParameter(0);
  double par2 = fit1->GetParameter(1);
  
  // create initial parameters for the fit of E_true vs E_rec
  double par1_init = -1.0*par1/par2;
  double par2_init = 1.0/par2;
  
  // fit the graph E_true vs E_rec
  //fit2->SetParameters(par1_init,par2_init);
  fit2->SetParameters(5.0,0.6);
  gr1->Fit("fit2");
  
  gMinuit->mnmatu(1);
  
  // calculate residuals
  int nDataPoints = gr1->GetN();
  double *Erec = gr1->GetX();
  double *Etrue = gr1->GetY();
  double *Erec_err = gr1->GetEX();
  
  double *residual = new double[nDataPoints];
  double *residual_err = new double[nDataPoints];
  
  for (int i = 0; i < nDataPoints; i++) {
    residual[i] = (Erec[i] - fit1->Eval(Etrue[i])) / fit1->Eval(Etrue[i]) * 100.0;
    residual_err[i] = Erec_err[i] / fit1->Eval(Etrue[i]) * 100.0;
  }
  
  TGraphErrors *grResiduals = new TGraphErrors(nDataPoints,Etrue,residual,0,residual_err);
  
  // return fit result
  *a = fit1->GetParameter(1);
  *b = fit1->GetParameter(0);
  *a_err = var2;
  *b_err = var1;
  *cor = correlation;
  
  return grResiduals;
}

double Sigma(double *x, double *par)
{
  return par[0]*x[0] + par[1] + par[5]*TMath::Sqrt(x[0]**2 * par[2] + par[3] + x[0]*par[4]) / x[0] * 100.0;
}

void SetGraphProperties(TGraphErrors *gr, char *title, char *xAxis, char *yAxis, int MarkerStyle = 20, double MarkerSize = 0.8)
{
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle(xAxis);
  gr->GetYaxis()->SetTitle(yAxis);
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerSize(MarkerSize);
  
  return;
}