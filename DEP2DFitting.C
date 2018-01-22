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

void DEP2DFitting()
{
  gStyle->SetPalette(1);

  TChain *t = new TChain("t","tree");

  t->Add("../../EnergyCalibration/analysis/V6/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V6/2448_noReclustering_Fiducial.root");

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,3500,200,0,15000);

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");

  TF2 *f2 = new TF2("f2",fun2,2000,2240,6400,7500,6);
  f2->SetParameters(200,2100,7000,70,200,-0.15);
  f2->SetParNames("A","E_I","E_S","#sigma_I","#sigma_S","#Theta");

  //f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,2060,2180);
  f2->SetParLimits(2,6700,7300);
  f2->SetParLimits(3,10,200);
  f2->SetParLimits(4,20,600);
  f2->SetParLimits(5,-0.2,-0.08);

  h2DSingle->Fit("f2","rn");

  double par[6];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();

  cout << "End of run " << "  (" << par[1] << " +- " << parErr[1] << ")" << "  (" << par[2] << " +- " << parErr[2] << ")" << endl;

  TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");
  f2->Draw("cont3 same");

  return;
}
