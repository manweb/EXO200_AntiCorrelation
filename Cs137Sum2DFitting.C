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

void Cs137Sum2DFitting()
{
  gStyle->SetPalette(1);
  
  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2450_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2469_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2473_noReclustering_Fiducial.root");

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,1500,200,0,5000);
  TH2F *h2DSingle_calib = new TH2F("h2DSingle_calib","scintillation vs ionization (calibrated)",200,0,1500,200,0,1500);

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");
  //t->Draw("(csc*0.3148-119.5):(ecrec*0.9647+67.4)>>h2DSingle_calib","nsite == 1");
  t->Draw("(csc*0.3194-143.2):(epcrec*0.9369+68.263)>>h2DSingle_calib","nsite == 1");

  TF2 *f2 = new TF2("f2",fun2,500,700,2000,3000,6);
  f2->SetParameters(40,600,2500,70,200,0.0);
  f2->SetParNames("A","E_I","E_S","#sigma_I","#sigma_S","#Theta");

  f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,540,640);
  f2->SetParLimits(2,2200,2700);
  f2->SetParLimits(3,10,100);
  f2->SetParLimits(4,20,600);

  h2DSingle->Fit("f2","r");
  
  TF2 *f2_2 = new TF2("f2_2",fun2,570,750,500,850,6);
  f2_2->SetParameters(40,660,660,70,70,0.3);
  f2_2->SetParNames("A","E_I","E_S","#sigma_I","#sigma_S","#Theta");
  
  f2_2->SetParLimits(0,0,1000);
  f2_2->SetParLimits(1,580,700);
  f2_2->SetParLimits(2,580,700);
  f2_2->SetParLimits(3,10,100);
  f2_2->SetParLimits(4,10,200);
  f2_2->SetParLimits(5,-1.0,-0.1);
  
  h2DSingle_calib->Fit("f2_2","r");

  double par[6];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();

  cout << "End of run (" << par[1] << " +- " << parErr[1] << ")" << "  (" << par[2] << " +- " << parErr[2] << ")" << endl;
  
  TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");
  f2->Draw("cont3 same");
  
  TCanvas *c2 = new TCanvas();
  h2DSingle_calib->Draw("colz");
  f2_2->Draw("cont3 same");

  return;
}
