
double fitFunction(double *x, double *par)
{
  double w = 18.7*0.9278 / 1000.0;

  //return TMath::Sqrt((w*par[0])*(w*par[0]) + (w*x[0]*par[1])*(w*x[0]*par[1]) + (par[2]*x[0])*(par[2]*x[0])) / x[0];
  //return TMath::Sqrt((w*par[0])*(w*par[0]) + w*par[1]*x[0] + (par[2]*x[0])*(par[2]*x[0])) / x[0];
  //return TMath::Sqrt(par[0] + par[1]*x[0] + par[2]*x[0]*x[0]) / x[0];
  //return par[2]/x[0] + par[1]/TMath::Sqrt(x[0]) + par[0];
  
  //return TMath::Sqrt(par[0]**2/x[0] + par[1]**2/x[0]**2);
  return TMath::Sqrt(par[0]**2/x[0] + par[1]**2/x[0]**2 + par[2]**2);
}

void EnergyResolution()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  const int nData = 4;
  double E[nData] = {661.657, 1173.228, 1332.492, 2614.511};
  //double R[nData] = {0.061, 0.0358, 0.0325, 0.0206}; // new shaping times
  //double Rerr[nData] = {0.0014, 0.00105, 0.0006, 0.0004}; // new shaping times

  //double R[nData] = {0.0607939, 0.0362046, 0.0320747, 0.0203501};
  //double Rerr[nData] = {0.00054, 0.00031, 0.00019, 0.00007};

  // calibrated (run by run to take systematics into account), theta = 0.1956
  /*double R[nData] = {0.064546233, 0.036823033, 0.032884422, 0.020087913};
  double Rerr[nData] = {0.001917017, 0.001297902, 0.001110586, 0.000612329};
  
  // calibrated, multi site (run by run to take systematics into account), theta = 0.2267
  double R[nData] = {0.078305133, 0.043414144, 0.0381218, 0.02373635};
  double Rerr[nData] = {0.002591615, 0.001513609, 0.001416906, 0.000864547};*/
  
  // Reprocessed data (01/30/12)
  // calibrated (run by run to take systematics into account), theta = 0.1916
  /*double R[nData] = {0.0635084, 0.036567456, 0.032714878, 0.019993588};
  double Rerr[nData] = {0.001919157, 0.001621318, 0.000903497, 0.000581574};
  
  // calibrated, multi site (run by run to take systematics into account), theta = 0.2164
  double R[nData] = {0.076727567, 0.0424315, 0.036940556, 0.022916663};
  double Rerr[nData] = {0.002117426, 0.001223671, 0.001447468, 0.000751454};*/
  
  // Reprocessed data (02/07/12)
  // calibrated (run by run to take systematics into account), theta = 0.1867
  //double R[nData] = {0.059052033, 0.035700989, 0.031326622, 0.018434475};
  //double Rerr[nData] = {0.0030079, 0.000710393, 0.001705939, 0.000644516};
  
  // calibrated, multi site (run by run to take systematics into account), theta = 0.2084
  //double R[nData] = {0.067239567, 0.040610133, 0.033665444, 0.020992375};
  //double Rerr[nData] = {0.002555952, 0.001297748, 0.000805116, 0.000840202};
  
  // uncalibrated (run by run to take systematics into account), theta = 0.1956
  //double R[nData] = {0.061021533, 0.036059122, 0.031771844, 0.01987885};
  //double Rerr[nData] = {0.001971915, 0.001189538, 0.001133059, 0.000583484};
  
  // calibrated (ionization, scintillation, anticorrelation), run by run, individual angles
  //double R[nData] = {0.055478933, 0.034435511, 0.030419689, 0.019862813};
  //double Rerr[nData] = {0.001743358, 0.00121196, 0.001099161, 0.000589919};
  
  // reprocessed, correction 2 u wire (03/21/12)
  /*double R2[nData] = {0.0573855, 0.033145, 0.0293416, 0.0161868};
  double R2err[nData] = {0.000498274, 0.00020741, 0.000110952, 7.63949e-05};
  
  // reprocessed (03/21/12)
  double R[nData] = {0.0574296, 0.0333193, 0.0295098, 0.0165676};
  double Rerr[nData] = {0.000497416, 0.000203371, 0.000112946, 7.96293e-05};
  
  double R[nData] = {0.0638152, 0.0370655, 0.0327581, 0.0183353};
  double Rerr[nData] = {0.00070252, 0.000205626, 0.000103656, 4.69745e-05};*/
  
// #### Calibration 03/22/12 #################################################
  
  // reprocessed, correction 2 u wire (03/21/12)
  double R2[nData] = {0.0568267, 0.03256, 0.0287372, 0.0162697};
  double R2err[nData] = {0.000483073, 0.000195144, 0.000107071, 6.62179e-05};
  
  // reprocessed (03/21/12)
  double R[nData] = {0.0570073, 0.0328359, 0.0289199, 0.0166462};
  double Rerr[nData] = {0.000484691, 0.000199294, 0.000108166, 6.90148e-05};
  
  //double R2[nData] = {0.0637251, 0.0366017, 0.0321563, 0.0184252};
  //double R2err[nData] = {0.000695513, 0.000200046, 9.91778e-05, 4.12159e-05};
  
// ###########################################################################

  TGraphErrors *gr = new TGraphErrors(nData,E,R,0,Rerr);
  //TGraph *gr = new TGraph(nData,E,R);
  TGraphErrors *gr2 = new TGraphErrors(nData,E,R2,0,R2err);

  gr->GetXaxis()->SetTitle("energy [keV]");
  gr->GetYaxis()->SetTitle("#sigma/E");

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.6);
  
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(0.6);

  TF1 *fit = new TF1("fit",fitFunction,100,3500,3);
  //fit->SetParNames("#sigma_{e}","F","k");
  //fit->SetParameters(1595, -0.1, 0.001);

  //fit->FixParameter(0,744);
  //fit->SetParLimits(2,0,1000);
  
  TF1 *fit2 = new TF1("fit2",fitFunction,100,3500,3);

  fit->SetLineWidth(1);
  
  fit2->SetLineWidth(1);
  fit2->SetLineStyle(2);

  gr->Fit("fit","r");
  gr2->Fit("fit2","r");

  fitter = TVirtualFitter::GetFitter();
  TMatrixD* COVM = new TMatrixD(2,2,fitter->GetCovarianceMatrix());
  double cov = COVM[0][1];
  cout << cov << endl;

  gMinuit->mnmatu(1);

  TCanvas *c1 = new TCanvas();
  gr->Draw("AZP");
  gr2->Draw("ZPsame");

return;
}
