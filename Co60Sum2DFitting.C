Double_t g2(Double_t *x, Double_t *par) {
   double A = par[0];
   double x_0 = par[1];
   double y_0 = par[2];
   double sigma_x = par[3];
   double sigma_y = par[4];
   double a = TMath::Cos(par[5])**2 / (2*sigma_x**2) + TMath::Sin(par[5])**2 / (2*sigma_y**2);
   double b = -1.0 * TMath::Sin(2*par[5]) / (4*sigma_x**2) + TMath::Sin(2*par[5]) / (4*sigma_y**2);
   double c = TMath::Sin(par[5])**2 / (2*sigma_x**2) + TMath::Cos(par[5])**2 / (2*sigma_y**2);

   double val = A*TMath::Exp(-1.0*(a*(x[0] - x_0)**2 + 2*b*(x[0] - x_0)*(x[1] - y_0) + c*(x[1] - y_0)**2));

   return val;
}

Double_t fun2(Double_t *x, Double_t *par) {
   Double_t *p1 = &par[0];
   Double_t *p2 = &par[6];
   Double_t result = g2(x,p1) + g2(x,p2);
   return result;
}

void Co60Sum2DFitting()
{
  gStyle->SetPalette(1);
  
  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2526_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2538_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2543_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2566_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2578_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2596_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2608_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2620_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2634_noReclustering_Fiducial.root");

  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,2000,200,0,10000);
  TH2F *h2DSingle_calib = new TH2F("h2DSingle_calib","scintillation vs ionization (calibrated)",200,0,2000,200,0,2000);

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");
  //t->Draw("(csc*0.3148-119.5):(ecrec*0.9647+67.4)>>h2DSingle_calib","nsite == 1");
  t->Draw("(csc*0.3194-143.2):(epcrec*0.9369+68.263)>>h2DSingle_calib","nsite == 1");

  TF2 *f2 = new TF2("f2",fun2,1000,1500,3200,6000,12);
  double params[12] = {40,1100,4000,70,200,0.0,40,1300,4600,70,200,0.0};
  f2->SetParameters(params);
  //f2->SetParNames("A1","E_I1","E_S1","#sigma_I1","#sigma_S1","#Theta1","A2","E_I2","E_S2","#sigma_I2","#sigma_S2","#Theta2");

  f2->SetParLimits(0,0,1000);
  f2->SetParLimits(1,1000,1200);
  f2->SetParLimits(2,3600,4400);
  f2->SetParLimits(3,10,100);
  f2->SetParLimits(4,20,600);
  f2->SetParLimits(6,0,1000);
  f2->SetParLimits(7,1250,1400);
  f2->SetParLimits(8,4400,5000);
  f2->SetParLimits(9,10,100);
  f2->SetParLimits(10,20,600);

  h2DSingle->Fit("f2","r");
  
  TF2 *f2_2 = new TF2("f2_2",fun2,1020,1550,950,1550,12);
  double params2[12] = {40,1160,1200,70,70,-0.79,40,1340,1340,70,70,-0.79};
  f2_2->SetParameters(params2);
  //f2->SetParNames("A1","E_I1","E_S1","#sigma_I1","#sigma_S1","#Theta1","A2","E_I2","E_S2","#sigma_I2","#sigma_S2","#Theta2");
  
  f2_2->SetParLimits(0,0,1000);
  f2_2->SetParLimits(1,1100,1240);
  f2_2->SetParLimits(2,1100,1240);
  f2_2->SetParLimits(3,10,100);
  f2_2->SetParLimits(4,10,200);
  f2_2->SetParLimits(5,-0.9,-0.1);
  f2_2->SetParLimits(6,0,1000);
  f2_2->SetParLimits(7,1280,1380);
  f2_2->SetParLimits(8,1250,1400);
  f2_2->SetParLimits(9,10,100);
  f2_2->SetParLimits(10,10,200);
  f2_2->SetParLimits(11,-0.9,-0.1);
  
  h2DSingle_calib->Fit("f2_2","r");

  double par[12];
  double *parErr;
  f2->GetParameters(par);
  parErr = f2->GetParErrors();
  
  cout << "End of run (" << par[1] << " +- " << parErr[1] << ")" << "  (" << par[2] << " +- " << parErr[2] << ")" << endl;
  cout << "End of run (" << par[7] << " +- " << parErr[7] << ")" << "  (" << par[8] << " +- " << parErr[8] << ")" << endl;
  
  TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");
  f2->Draw("cont3 same");
  
  TCanvas *c2 = new TCanvas();
  h2DSingle_calib->Draw("colz");
  f2_2->Draw("cont3 same");

  return;
}
