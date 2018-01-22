double LineFitError(double *x, double *par)
{
  double sy1 = par[6];
  double sy2 = par[7];
  double sy3 = par[8];
  double sy4 = par[9];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];

  double bracketX = 1.0/4.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4)); // [x]
  double bracketX2 = 1.0/4.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4)); // [x2]
  double bracket1 = 1.0/4.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/4.0 / bracket1;
  double sb2 = 1.0/4.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[0] + par[1]*x[0] + par[10]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
}

double LineFitErrorMulti(double *x, double *par)
{
  double sy1 = par[5];
  double sy2 = par[6];
  double sy3 = par[7];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];

  double bracketX = 1.0/3.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3)); // [x]
  double bracketX2 = 1.0/3.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3)); // [x2]
  double bracket1 = 1.0/3.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/3.0 / bracket1;
  double sb2 = 1.0/3.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[0] + par[1]*x[0] + par[8]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2);
}

double LineFitResidual(double *x, double *par)
{
  double sy1 = par[6];
  double sy2 = par[7];
  double sy3 = par[8];
  double sy4 = par[9];

  // create covariance of the line fit
  double x1 = par[2];
  double x2 = par[3];
  double x3 = par[4];
  double x4 = par[5];

  double bracketX = 1.0/4.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3) + x4/(sy4*sy4)); // [x]
  double bracketX2 = 1.0/4.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3) + x4*x4/(sy4*sy4)); // [x2]
  double bracket1 = 1.0/4.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3) + 1.0/(sy4*sy4)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/4.0 / bracket1;
  double sb2 = 1.0/4.0 * bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return par[10]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2) / x[0] * 100.0;
}

double LineFitResidualMulti(double *x, double *par)
{
  double sy1 = par[7];
  double sy2 = par[8];
  double sy3 = par[9];

  // create covariance of the line fit
  double x1 = par[4];
  double x2 = par[5];
  double x3 = par[6];

  double bracketX = 1.0/3.0 * (x1/(sy1*sy1) + x2/(sy2*sy2) + x3/(sy3*sy3)); // [x]
  double bracketX2 = 1.0/3.0 * (x1*x1/(sy1*sy1) + x2*x2/(sy2*sy2) + x3*x3/(sy3*sy3)); // [x2]
  double bracket1 = 1.0/3.0 * (1.0/(sy1*sy1) + 1.0/(sy2*sy2) + 1.0/(sy3*sy3)); // [1]

  double cov = -1.0 * (bracketX / bracket1);
  double sa2 = 1.0/3.0 / bracket1;
  double sb2 = 1.0/3.0* bracket1 / (bracketX2*bracket1 - bracketX*bracketX);

  return (par[2] + par[3]*x[0] + par[10]*TMath::Sqrt((x[0]+cov)*(x[0]+cov)*sb2 + sa2) - par[0] - par[1]*x[0]) / x[0] * 100.0;
}

void EnergyCalibration_4()
{
/*  double Etrue_single[4] = {511.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {468.03, 1203.23, 1377.2, 2751.33};
  double ErecErr_single[4] = {5.58, 8.74119, 7.51124, 15.9804};

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1702.1};
  double ErecErr_pp[1] = {12.67};
*/
/*  double Etrue_single[4] = {511.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {485.54, 1201.17, 1374.91, 2748.86}; // no reclustering (fiducial volume)
  //double Erec_single[4] = {485.54, 1195.45, 1371.13, 2760.85}; // reclustering (fiducial volume)
  //double Erec_single[4] = {463.86, 1201.17, 1374.91, 2748.86}; // old value for 511 peak
  double ErecErr_single[4] = {5.04, 11.01, 9.3, 15.04}; // no reclustering (fiducial volume)
  //double ErecErr_single[4] = {5.04, 11.32, 10.8, 20.43}; // reclustering (fiducial volume)
  //double ErecErr_single[4] = {5.46, 11.01, 9.3, 15.04}; // old value for 511 peak

  double Etrue_multi[3] = {1173.2, 1332.5, 2614.0};
  double Erec_multi[3] = {1155.17, 1312.04, 2707.93}; // gauss + erf fit
  //double Erec_multi[3] = {1139.13, 1323.46, 2687.8}; // gauss fit on sum spectrum
  //double Erec_multi[3] = {1152.15, 1308.61, 2678.98}; // reclustering
  double ErecErr_multi[3] = {7.43412, 8.44368, 12.0313}; // gauss + erf fit
  //double ErecErr_multi[3] = {7.68, 9.42, 5.4}; // gauss fit on sum spectrum
  //double ErecErr_multi[3] = {6.34031, 7.20133, 9.97868}; // reclustering

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1708.26};
  double ErecErr_pp[1] = {10.15};
*/
// **** New shaping times *************************************************************************

/*  double Etrue_single[4] = {511.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {459.88, 1137.75, 1305.1, 2640.21}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {3.13, 1.48, 2.62, 8.73}; // no reclustering (fiducial volume)

  double Etrue_multi[3] = {1173.2, 1332.5, 2614.0};
  double Erec_multi[3] = {1097.78, 1246.86, 2607.92}; // gauss + erf fit
  double ErecErr_multi[3] = {6.16, 7.0, 12.79}; // gauss + erf fit

  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {1622.0};
  double ErecErr_pp[1] = {8.53};
*/
// ************************************************************************************************

// **** MC new shaping times **********************************************************************

/*  double Etrue_single[4] = {511.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {499.04, 1136.07, 1288.54, 2514.79}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {0.539, 0.134, 0.145, 0.502}; // no reclustering (fiducial volume)

  double Etrue_multi[3] = {1173.2, 1332.5, 2614.0};
  double Erec_multi[3] = {1135.83, 1290.08, 2548.01}; // gauss + erf fit
  double ErecErr_multi[3] = {0.081, 0.091, 0.578}; // gauss + erf fit
*/
// DEP from anticorrelation
  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {2551.84};
  double ErecErr_pp[1] = {10.0};

// DEP from scintillation
/*  double Etrue_pp[1] = {1592.0};
  double Erec_pp[1] = {4873.9};
  double ErecErr_pp[1] = {8.68};

  double Etrue_pp2[1] = {1592.0};
  double Erec_pp2[1] = {5041.21};
  double ErecErr_pp2[1] = {9.11};*/
// ************************************************************************************************

// Anticorrelated source positions
// theta = 0.1944
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {1092.08, 1920.83, 2181.47, 4255.39}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {3.39, 3.15, 3.62, 2.34}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {662.0,1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {1072.73, 1910.78, 2172.46, 4244.89}; // gauss + erf fit
  double ErecErr_multi[4] = {1.88, 2.21, 3.14, 2.71}; // gauss + erf fit
*/
// theta = 0.1622
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {1004.84, 1792.14, 2037.76, 3994.44}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {1.31, 4.04, 4.07, 5.05}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {662.0,1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {988.83, 1775.68, 2022.11, 3978.12}; // gauss + erf fit
  double ErecErr_multi[4] = {1.63, 3.06, 3.64, 4.33}; // gauss + erf fit
*/
// selected runs, theta = 0.193
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {1077.1, 1912.72, 2172.37, 4240.15};
  double ErecErr_single[4] = {2.86, 4.76, 4.84, 5.4};

  double Etrue_multi[4] = {662.0,1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {1062.44, 1902.49, 2163.29, 4231.53};
  double ErecErr_multi[4] = {3.18, 3.3, 4.6, 5.26};
*/
// selected runs, theta = 0.1956
/*  double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1079.45, 1920.9, 2182.32, 4260.55};
  double ErecErr_single[4] = {3.7, 4.86, 5.13, 5.75};

// selected runs, multi site, theta = 0.2267
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1137.21, 2034.15, 2310.38, 4499.63};
  double ErecErr_multi[4] = {1.79, 4.71, 5.60, 6.62};
*/
// Reprocessed data (01/30/12)
  // selected runs, theta = 0.1916
  /*double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1075.08, 1917.02, 2179.04, 4264.98};
  double ErecErr_single[4] = {2.73, 3.76, 3.87, 5.81};
  
  // selected runs, multi site, theta = 0.2164
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1119.09, 2005.24, 2278.72, 4453.18};
  double ErecErr_multi[4] = {1.27, 3.24, 3.26, 7.03};*/
  
// Reprocessed data (02/07/12)
  // selected runs, theta = 0.1867
  /*double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1070.95, 1909.21, 2169.47, 4245.2};
  double ErecErr_single[4] = {2.49, 4.09, 3.83, 3.78};
  
  // selected runs, multi site, theta = 0.2084
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1104.61, 1980.79, 2250.74, 4398.56};
  double ErecErr_multi[4] = {4.35, 3.74, 3.45, 4.92};*/
  
// Reprocessed data, alpha branch (02/14/12)
  // selected runs, theta = 0.1867
  /*double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1074.94, 1919.66, 2181.43, 4260.96};
  double ErecErr_single[4] = {1.95, 2.37, 2.49, 5.16};
  
  // selected runs, multi site, theta = 0.2084
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1108.55, 1991.98, 2263.64, 4415.11};
  double ErecErr_multi[4] = {1.97, 2.48, 2.96, 6.28};*/
  
// Reprocessed data, alpha branch (02/22/12)
  // selected runs, theta = 0.1860
  double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1073.55, 1917.52, 2179.38, 4258.49};
  double ErecErr_single[4] = {1.93, 3.66, 3.34, 6.28};
  
  // selected runs, multi site, theta = 0.2064
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1103.94, 1984.67, 2255.21, 4402.39};
  double ErecErr_multi[4] = {2.37, 3.76, 4.59, 6.90};

// Anticorrelated source positions, calibrated, individual rotation angles
/*  double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {850.51, 1575.62, 1767.39, 3569.32};
  double ErecErr_single[4] = {3.56, 3.48, 3.67, 5.0};
*/
// Anticorrelated source positions with fScintillationCluster.fAlgorithmUsed == 0 applied
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {1083.73, 1919.8, 2179.91, 4253.74}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {4.02, 3.61, 2.92, 4.06}; // no reclustering (fiducial volume)

  double Etrue_multi[4] = {662.0,1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {1069.26, 1909.84, 2171.32, 4244.4}; // gauss + erf fit
  double ErecErr_multi[4] = {3.45, 1.78, 2.68, 2.78}; // gauss + erf fit
*/
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {964.13, 1782.13, 2026.29, 4253.74};
  double ErecErr_single[4] = {2.68, 3.01, 2.73, 4.06};

  double Etrue_multi[4] = {662.0,1173.2, 1332.5, 2614.0};
  double Erec_multi[4] = {944.48, 1764.79, 2009.49, 4244.4};
  double ErecErr_multi[4] = {2.27, 1.87, 2.38, 2.78};
*/
// Scintillation source positions
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {2509.31, 4180.87, 4682.36, 8608.48}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {17.67, 17.33, 11.27, 28.63}; // no reclustering (fiducial volume)
*/
// Scintillation source positions (12/27/11)
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {2442.89, 4155.65, 4661.13, 8593.81}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {14.8, 21.19, 18, 33.26}; // no reclustering (fiducial volume)*/
  
  // Scintillation source positions, 2D fits, alpha branch (02/22/12)
/*  double Etrue_single[4] = {662.0, 1173.2, 1332.5, 2614.0};
  double Erec_single[4] = {2479.77, 4177.92, 4674.51, 8595.85}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {32.96, 44.90, 44.93, 47.67}; // no reclustering (fiducial volume)
*/
  // ####################################################
  // Calibration with cut on number of u wires. The single site correspond to
  // nU = 1, the multi site corresponds to nU = 2.
  /*double Etrue_single[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_single[4] = {1073.31, 1921.8, 2185.34, 4276.07}; // no reclustering (fiducial volume)
  double ErecErr_single[4] = {1.98, 3.85, 4.04, 6.4974}; // no reclustering (fiducial volume)
  
  double Etrue_multi[4] = {661.657, 1173.228, 1332.492, 2614.511};
  double Erec_multi[4] = {1073.9, 1913.16, 2173.61, 4248.51};
  double ErecErr_multi[4] = {2.21, 4.12, 3.37, 6.40};*/
  
  // ####################################################
  
  TGraphErrors *gr1 = new TGraphErrors(4,Erec_single,Etrue_single,ErecErr_single,0); // E_true vs E_rec (single site)
  TGraphErrors *gr2 = new TGraphErrors(4,Etrue_single,Erec_single,0,ErecErr_single); // E_rec vs E_true (singel site)
  TGraphErrors *gr3 = new TGraphErrors(1,Erec_pp,Etrue_pp,0,ErecErr_pp); // double escape peak
  //TGraphErrors *gr31 = new TGraphErrors(1,Etrue_pp2,Erec_pp2,0,ErecErr_pp2); // double escape peak
  TGraphErrors *gr4 = new TGraphErrors(4,Erec_multi,Etrue_multi,ErecErr_multi,0); // E_true vs E_rec (multi site)
  TGraphErrors *gr5 = new TGraphErrors(4,Etrue_multi,Erec_multi,0,ErecErr_multi); // E_rec vs E_true (multi site)

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.6);

  gr1->GetXaxis()->SetTitle("rec energy");
  gr1->GetYaxis()->SetTitle("true energy");

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.6);

  gr2->GetXaxis()->SetTitle("true energy");
  gr2->GetYaxis()->SetTitle("rec energy");

  gr3->SetMarkerStyle(20);
  gr3->SetMarkerSize(0.6);
  gr3->SetMarkerColor(kRed);

  gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(0.6);

  gr4->GetXaxis()->SetTitle("rec energy");
  gr4->GetYaxis()->SetTitle("true energy");

  gr5->SetMarkerStyle(20);
  gr5->SetMarkerSize(0.6);

  gr5->GetXaxis()->SetTitle("true energy");
  gr5->GetYaxis()->SetTitle("rec energy");

// **** Single Site ***************************************************************************************
// Fit E_rec vs E_true
  TF1 *fit = new TF1("fit","[0]+[1]*x",0,4500);
  fit->SetLineWidth(1);

  gr2->Fit("fit","r");

  double par[2];
  fit->GetParameters(par);

  double a = 1.0/par[1];
  double b = -1.0*a*par[0];

  TF1 *ErrLine1_1 = new TF1("ErrLine1_1",LineFitError,0,4500,11);

  ErrLine1_1->SetParameter(0,par[0]);
  ErrLine1_1->SetParameter(1,par[1]);
  ErrLine1_1->SetParameter(2,Etrue_single[0]);
  ErrLine1_1->SetParameter(3,Etrue_single[1]);
  ErrLine1_1->SetParameter(4,Etrue_single[2]);
  ErrLine1_1->SetParameter(5,Etrue_single[3]);
  ErrLine1_1->SetParameter(6,ErecErr_single[0]);
  ErrLine1_1->SetParameter(7,ErecErr_single[1]);
  ErrLine1_1->SetParameter(8,ErecErr_single[2]);
  ErrLine1_1->SetParameter(9,ErecErr_single[3]);
  ErrLine1_1->SetParameter(10,1);

  ErrLine1_1->SetLineWidth(1);

  TF1 *ErrLine2_1 = new TF1("ErrLine2_1",LineFitError,0,4500,11);

  ErrLine2_1->SetParameter(0,par[0]);
  ErrLine2_1->SetParameter(1,par[1]);
  ErrLine2_1->SetParameter(2,Etrue_single[0]);
  ErrLine2_1->SetParameter(3,Etrue_single[1]);
  ErrLine2_1->SetParameter(4,Etrue_single[2]);
  ErrLine2_1->SetParameter(5,Etrue_single[3]);
  ErrLine2_1->SetParameter(6,ErecErr_single[0]);
  ErrLine2_1->SetParameter(7,ErecErr_single[1]);
  ErrLine2_1->SetParameter(8,ErecErr_single[2]);
  ErrLine2_1->SetParameter(9,ErecErr_single[3]);
  ErrLine2_1->SetParameter(10,-1);

  ErrLine2_1->SetLineWidth(1);

// Fit E_true vs E_rec
  TF1 *LineFit = new TF1("LineFit","[0]+[1]*x",0,9400);

  LineFit->SetLineWidth(1);

  //LineFit->SetParameter(0,b);
  //LineFit->SetParameter(1,a);
  
  LineFit->SetParameter(0,2.0);
  LineFit->SetParameter(1,0.61);

  gr1->Fit("LineFit","r");
  
  fitter = TVirtualFitter::GetFitter();
  TMatrixD* COVM = new TMatrixD(2,2,fitter->GetCovarianceMatrix());
  double cov = COVM[0][1];
  cout << cov << endl;
  
  gMinuit->mnmatu(1);

  double par2[2];
  LineFit->GetParameters(par2);

  TF1 *ErrLine1 = new TF1("ErrLine1",LineFitError,0,9400,11);

  ErrLine1->SetParameter(0,par2[0]);
  ErrLine1->SetParameter(1,par2[1]);
  ErrLine1->SetParameter(2,Etrue_single[0]);
  ErrLine1->SetParameter(3,Etrue_single[1]);
  ErrLine1->SetParameter(4,Etrue_single[2]);
  ErrLine1->SetParameter(5,Etrue_single[3]);
  ErrLine1->SetParameter(6,ErecErr_single[0]*par2[1]);
  ErrLine1->SetParameter(7,ErecErr_single[1]*par2[1]);
  ErrLine1->SetParameter(8,ErecErr_single[2]*par2[1]);
  ErrLine1->SetParameter(9,ErecErr_single[3]*par2[1]);
  ErrLine1->SetParameter(10,1);

  ErrLine1->SetLineWidth(1);

  TF1 *ErrLine2 = new TF1("ErrLine2",LineFitError,0,9400,11);

  ErrLine2->SetParameter(0,par2[0]);
  ErrLine2->SetParameter(1,par2[1]);
  ErrLine2->SetParameter(2,Etrue_single[0]);
  ErrLine2->SetParameter(3,Etrue_single[1]);
  ErrLine2->SetParameter(4,Etrue_single[2]);
  ErrLine2->SetParameter(5,Etrue_single[3]);
  ErrLine2->SetParameter(6,ErecErr_single[0]*par2[1]);
  ErrLine2->SetParameter(7,ErecErr_single[1]*par2[1]);
  ErrLine2->SetParameter(8,ErecErr_single[2]*par2[1]);
  ErrLine2->SetParameter(9,ErecErr_single[3]*par2[1]);
  ErrLine2->SetParameter(10,-1);

  ErrLine2->SetLineWidth(1);

// **** Multi site ****************************************************************************************
// Fit E_rec vs E_true
  TF1 *fitMulti = new TF1("fitMulti","[0]+[1]*x",0,4500);
  fitMulti->SetLineWidth(1);

  gr5->Fit("fitMulti","r");

  double parMulti[2];
  fitMulti->GetParameters(parMulti);

  double aMulti = 1.0/parMulti[1];
  double bMulti = -1.0*aMulti*parMulti[0];

  TF1 *ErrLineMulti1_1 = new TF1("ErrLineMulti1_1",LineFitErrorMulti,0,4500,9);

  ErrLineMulti1_1->SetParameter(0,parMulti[0]);
  ErrLineMulti1_1->SetParameter(1,parMulti[1]);
  ErrLineMulti1_1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti1_1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti1_1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti1_1->SetParameter(5,ErecErr_multi[0]);
  ErrLineMulti1_1->SetParameter(6,ErecErr_multi[1]);
  ErrLineMulti1_1->SetParameter(7,ErecErr_multi[2]);
  ErrLineMulti1_1->SetParameter(8,1);

  ErrLineMulti1_1->SetLineWidth(1);

  TF1 *ErrLineMulti2_1 = new TF1("ErrLineMulti2_1",LineFitErrorMulti,0,4500,9);

  ErrLineMulti2_1->SetParameter(0,parMulti[0]);
  ErrLineMulti2_1->SetParameter(1,parMulti[1]);
  ErrLineMulti2_1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti2_1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti2_1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti2_1->SetParameter(5,ErecErr_multi[0]);
  ErrLineMulti2_1->SetParameter(6,ErecErr_multi[1]);
  ErrLineMulti2_1->SetParameter(7,ErecErr_multi[2]);
  ErrLineMulti2_1->SetParameter(8,-1);

  ErrLineMulti2_1->SetLineWidth(1);

// Fit E_true vs E_rec
  TF1 *LineFitMulti = new TF1("LineFitMulti","[0]+[1]*x",0,4500);

  LineFitMulti->SetLineWidth(1);

  LineFitMulti->SetParameter(0,bMulti);
  LineFitMulti->SetParameter(1,aMulti);

  gr4->Fit("LineFit","r");
  
  gMinuit->mnmatu(1);

  double parMulti2[2];
  LineFitMulti->GetParameters(parMulti2);

  TF1 *ErrLineMulti1 = new TF1("ErrLineMulti1",LineFitErrorMulti,0,4500,9);

  ErrLineMulti1->SetParameter(0,parMulti2[0]);
  ErrLineMulti1->SetParameter(1,parMulti2[1]);
  ErrLineMulti1->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti1->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti1->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti1->SetParameter(5,ErecErr_multi[0]*parMulti2[1]);
  ErrLineMulti1->SetParameter(6,ErecErr_multi[1]*parMulti2[1]);
  ErrLineMulti1->SetParameter(7,ErecErr_multi[2]*parMulti2[1]);
  ErrLineMulti1->SetParameter(8,1);

  ErrLineMulti1->SetLineWidth(1);

  TF1 *ErrLineMulti2 = new TF1("ErrLineMulti2",LineFitErrorMulti,0,4500,11);

  ErrLineMulti2->SetParameter(0,parMulti2[0]);
  ErrLineMulti2->SetParameter(1,parMulti2[1]);
  ErrLineMulti2->SetParameter(2,Etrue_multi[0]);
  ErrLineMulti2->SetParameter(3,Etrue_multi[1]);
  ErrLineMulti2->SetParameter(4,Etrue_multi[2]);
  ErrLineMulti2->SetParameter(5,ErecErr_multi[0]*parMulti2[1]);
  ErrLineMulti2->SetParameter(6,ErecErr_multi[1]*parMulti2[1]);
  ErrLineMulti2->SetParameter(7,ErecErr_multi[2]*parMulti2[1]);
  ErrLineMulti2->SetParameter(8,-1);

  ErrLineMulti2->SetLineWidth(1);

// **** Residuals ***********************************************************************************************
// Single site
  // calculate residual for each data point
  double r1 = (Erec_single[0] - fit->Eval(Etrue_single[0]))/fit->Eval(Etrue_single[0]) * 100.0;
  double r2 = (Erec_single[1] - fit->Eval(Etrue_single[1]))/fit->Eval(Etrue_single[1]) * 100.0;
  double r3 = (Erec_single[2] - fit->Eval(Etrue_single[2]))/fit->Eval(Etrue_single[2]) * 100.0;
  double r4 = (Erec_single[3] - fit->Eval(Etrue_single[3]))/fit->Eval(Etrue_single[3]) * 100.0;

  double rErr1 = ErecErr_single[0]/fit->Eval(Etrue_single[0]) * 100.0;
  double rErr2 = ErecErr_single[1]/fit->Eval(Etrue_single[1]) * 100.0;
  double rErr3 = ErecErr_single[2]/fit->Eval(Etrue_single[2]) * 100.0;
  double rErr4 = ErecErr_single[3]/fit->Eval(Etrue_single[3]) * 100.0;

  double residual[4] = {r1, r2, r3, r4};
  double residualErr[4] = {rErr1, rErr2, rErr3, rErr4};
  TGraphErrors *grResidual_single = new TGraphErrors(4,Etrue_single,residual,0,residualErr);

  grResidual_single->SetMarkerStyle(20);
  grResidual_single->SetMarkerSize(0.8);
  grResidual_single->SetLineStyle(21);

  grResidual_single->GetXaxis()->SetTitle("true energy [keV]");
  grResidual_single->GetYaxis()->SetTitle("residual [%]");

  grResidual_single->GetYaxis()->SetRangeUser(-8.0,8.0);

// Multi site
  // calculate residual for each data point
  double rMulti1 = (Erec_multi[0] - fit->Eval(Etrue_multi[0]))/fit->Eval(Etrue_multi[0]) * 100.0;
  double rMulti2 = (Erec_multi[1] - fit->Eval(Etrue_multi[1]))/fit->Eval(Etrue_multi[1]) * 100.0;
  double rMulti3 = (Erec_multi[2] - fit->Eval(Etrue_multi[2]))/fit->Eval(Etrue_multi[2]) * 100.0;
  double rMulti4 = (Erec_multi[3] - fit->Eval(Etrue_multi[3]))/fit->Eval(Etrue_multi[3]) * 100.0;

  double rErrMulti1 = ErecErr_multi[0]/fit->Eval(Etrue_multi[0]) * 100.0;
  double rErrMulti2 = ErecErr_multi[1]/fit->Eval(Etrue_multi[1]) * 100.0;
  double rErrMulti3 = ErecErr_multi[2]/fit->Eval(Etrue_multi[2]) * 100.0;
  double rErrMulti4 = ErecErr_multi[3]/fit->Eval(Etrue_multi[3]) * 100.0;

  double residualMulti[4] = {rMulti1, rMulti2, rMulti3, rMulti4};
  double residualErrMulti[4] = {rErrMulti1, rErrMulti2, rErrMulti3, rErrMulti4};
  TGraphErrors *grResidual_multi = new TGraphErrors(4,Etrue_multi,residualMulti,0,residualErrMulti);

  grResidual_multi->SetMarkerStyle(22);
  grResidual_multi->SetMarkerSize(0.8);
  grResidual_multi->SetLineStyle(21);

// Double escape Peak
  double r_pp = (Erec_pp[0] - fit->Eval(Etrue_pp[0]))/fit->Eval(Etrue_pp[0]) * 100.0;
  double rErr_pp = ErecErr_pp[0]/fit->Eval(Etrue_pp[0]) * 100.0;

  double residual_pp[1] = {r_pp};
  double residualErr_pp[1] = {rErr_pp};
  TGraphErrors *grResidual_pp = new TGraphErrors(1,Etrue_pp,residual_pp,0,residualErr_pp);

  grResidual_pp->SetMarkerStyle(21);

// Error band single site
  TF1 *ResidualLine1 = new TF1("ResidualLine1",LineFitResidual,0,3000,11);

  ResidualLine1->SetParameter(0,par[0]);
  ResidualLine1->SetParameter(1,par[1]);
  ResidualLine1->SetParameter(2,Etrue_single[0]);
  ResidualLine1->SetParameter(3,Etrue_single[1]);
  ResidualLine1->SetParameter(4,Etrue_single[2]);
  ResidualLine1->SetParameter(5,Etrue_single[3]);
  ResidualLine1->SetParameter(6,ErecErr_single[0]);
  ResidualLine1->SetParameter(7,ErecErr_single[1]);
  ResidualLine1->SetParameter(8,ErecErr_single[2]);
  ResidualLine1->SetParameter(9,ErecErr_single[3]);
  ResidualLine1->SetParameter(10,1);

  ResidualLine1->SetLineWidth(1);

  TF1 *ResidualLine2 = new TF1("ResidualLine2",LineFitResidual,0,3000,11);

  ResidualLine2->SetParameter(0,par[0]);
  ResidualLine2->SetParameter(1,par[1]);
  ResidualLine2->SetParameter(2,Etrue_single[0]);
  ResidualLine2->SetParameter(3,Etrue_single[1]);
  ResidualLine2->SetParameter(4,Etrue_single[2]);
  ResidualLine2->SetParameter(5,Etrue_single[3]);
  ResidualLine2->SetParameter(6,ErecErr_single[0]);
  ResidualLine2->SetParameter(7,ErecErr_single[1]);
  ResidualLine2->SetParameter(8,ErecErr_single[2]);
  ResidualLine2->SetParameter(9,ErecErr_single[3]);
  ResidualLine2->SetParameter(10,-1);

  ResidualLine2->SetLineWidth(1);

// Error band multi site
  TF1 *ResidualLine1Multi = new TF1("ResidualLine1Multi",LineFitResidualMulti,0,3000,11);

  ResidualLine1Multi->SetParameter(0,par[0]);
  ResidualLine1Multi->SetParameter(1,par[1]);
  ResidualLine1Multi->SetParameter(2,parMulti[0]);
  ResidualLine1Multi->SetParameter(3,parMulti[1]);
  ResidualLine1Multi->SetParameter(4,Etrue_multi[0]);
  ResidualLine1Multi->SetParameter(5,Etrue_multi[1]);
  ResidualLine1Multi->SetParameter(6,Etrue_multi[2]);
  ResidualLine1Multi->SetParameter(7,ErecErr_multi[0]);
  ResidualLine1Multi->SetParameter(8,ErecErr_multi[1]);
  ResidualLine1Multi->SetParameter(9,ErecErr_multi[2]);
  ResidualLine1Multi->SetParameter(10,1);

  ResidualLine1Multi->SetLineWidth(1);
  ResidualLine1Multi->SetLineStyle(2);

  TF1 *ResidualLine2Multi = new TF1("ResidualLine2Multi",LineFitResidualMulti,0,3000,11);

  ResidualLine2Multi->SetParameter(0,par[0]);
  ResidualLine2Multi->SetParameter(1,par[1]);
  ResidualLine2Multi->SetParameter(2,parMulti[0]);
  ResidualLine2Multi->SetParameter(3,parMulti[1]);
  ResidualLine2Multi->SetParameter(4,Etrue_multi[0]);
  ResidualLine2Multi->SetParameter(5,Etrue_multi[1]);
  ResidualLine2Multi->SetParameter(6,Etrue_multi[2]);
  ResidualLine2Multi->SetParameter(7,ErecErr_multi[0]);
  ResidualLine2Multi->SetParameter(8,ErecErr_multi[1]);
  ResidualLine2Multi->SetParameter(9,ErecErr_multi[2]);
  ResidualLine2Multi->SetParameter(10,-1);

  ResidualLine2Multi->SetLineWidth(1);
  ResidualLine2Multi->SetLineStyle(2);

  cout << "Scintillation at 511 keV: " << fit->Eval(511) << " + " << ErrLine1_1->Eval(511) - fit->Eval(511) << " - " << fit->Eval(511) - ErrLine2_1->Eval(511) << endl;

  TCanvas *c1 = new TCanvas();
  gr2->Draw("AP");
  //gr3->Draw("Psame");
  //gr31->Draw("Psame");
  //ErrLine1_1->Draw("same");
  //ErrLine2_1->Draw("same");
  
  gr4->SetMarkerColor(kBlue);

  TCanvas *c2 = new TCanvas();
  gr1->Draw("AP");
  //gr3->Draw("Psame");
  gr4->Draw("Psame");
  //LineFit->Draw("same");
  //ErrLine1->Draw("same");
  //ErrLine2->Draw("same");

  TCanvas *c3 = new TCanvas();
  grResidual_single->Draw("AZP");
  grResidual_multi->Draw("ZPsame");
  //grResidual_pp->Draw("Psame");
  ResidualLine1->Draw("same");
  ResidualLine2->Draw("same");
  //ResidualLine1Multi->Draw("same");
  //ResidualLine2Multi->Draw("same");

  grResidual_single->GetXaxis()->SetLimits(0,3000);
  c3->Update();

  TCanvas *c4 = new TCanvas();
  gr5->Draw("AP");
  //ErrLineMulti1_1->Draw("same");
  //ErrLineMulti2_1->Draw("same");

  TCanvas *c5 = new TCanvas();
  gr4->Draw("AP");
  //ErrLineMulti1->Draw("same");
  //ErrLineMulti2->Draw("same");

  return;
}
