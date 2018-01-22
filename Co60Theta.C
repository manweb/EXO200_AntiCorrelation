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

  return gauss1 + erf1 + gauss2 + erf2;
}

void Co60Theta()
{
  gStyle->SetPalette(1);

  TChain *t = new TChain("t","tree");
  //t->Add("../../EnergyCalibration/analysis/V6/2496_noReclustering_Fiducial.root");
  /*t->Add("../../EnergyCalibration/analysis/V7/2526_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2538_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2543_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2566_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2578_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2596_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2608_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2620_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2634_noReclustering_Fiducial.root");*/
  
  // Reprocessed data (01/30/12)
  /*t->Add("../../EnergyCalibration/analysis/V8/2526_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2538_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2543_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2566_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2578_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2596_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2608_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2620_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2634_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2635_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2640_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2646_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2653_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2667_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2683_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2689_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2708_noReclustering_Fiducial.root");*/
  
  // Reprocessed data (02/07/12)
  t->Add("../../EnergyCalibration/analysis/V9/2526_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2538_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2543_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2566_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2578_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2596_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2608_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2620_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2634_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2635_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2640_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2646_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2653_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2667_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2683_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2689_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2708_noReclustering_Fiducial.root");

  /*const int nrHistos = 20;

  TH1F *hSingle[nrHistos];
  TH1F *hMulti[nrHistos];
  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,2000,200,0,10000);
  TH2F *h2DMulti = new TH2F("h2DMulti","scintillation vs ionization",200,0,2000,200,0,10000);

  TF1 *fitSingle[nrHistos];
  TF1 *fitMulti[nrHistos];

  double *angle = new double[nrHistos];
  double *res1 = new double[nrHistos];
  double *res2 = new double[nrHistos];
  double *resErr1 = new double[nrHistos];
  double *resErr2 = new double[nrHistos];

  double thetaRange1 = 0.1;
  double thetaRange2 = 1.0; //TMath::Pi()/9.0;
  double dTheta = thetaRange2 - thetaRange1;
  double stp = dTheta/nrHistos;

  for (int i = 0; i < nrHistos; i++) {
     char hSingleName[100];
     char hSingleTitle[100];
     char hMultiName[100];
     char hMultiTitle[100];

     double theta = thetaRange1 + stp*i;

     sprintf(hSingleName,"hSingle_%i",i);
     sprintf(hSingleTitle,"%.0f rotation angle",180.0/TMath::Pi()*theta);
     sprintf(hMultiName,"hMulti_%i",i);
     sprintf(hMultiTitle,"%.0f rotation angle",180.0/TMath::Pi()*theta);

     hSingle[i] = new TH1F(hSingleName,hSingleTitle,200,0,3500);
     hMulti[i] = new TH1F(hMultiName,hMultiTitle,200,0,3500);

     char SingleCMD[200];
     char MultiCMD[200];
     sprintf(SingleCMD,"(csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f) >> %s",theta,theta,hSingleName);
     sprintf(MultiCMD,"csc*TMath::Sin(%.2f) + ecrec*TMath::Cos(%.2f) >> %s",theta,theta,hMultiName);

     t->Draw(SingleCMD,"nsite == 1");
     t->Draw(MultiCMD,"nsite > 1");

     char fitSingleName[100];
     sprintf(fitSingleName,"fit_%i",i);

     double cscOff = 0;

     double fitRange1 = (950-cscOff)*TMath::Sin(theta) + 950*TMath::Cos(theta);
     double fitRange2 = (1500-cscOff)*TMath::Sin(theta) + 1500*TMath::Cos(theta);
     double Peak1 = (1120-cscOff)*TMath::Sin(theta) + 1120*TMath::Cos(theta);
     double Peak2 = (1300-cscOff)*TMath::Sin(theta) + 1300*TMath::Cos(theta);

     fitSingle[i] = new TF1(fitSingleName,fitFunction,fitRange1,fitRange2,7);
     fitSingle[i]->SetParNames("A1","E1","#sigma","R1","R2","E2","#sigma2");

     fitSingle[i]->SetParameters(3000,Peak1,75,0.6,0.8,Peak2,75);
     fitSingle[i]->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
     fitSingle[i]->SetParLimits(2,30,120);
     fitSingle[i]->SetParLimits(3,0.05,1.0);
     fitSingle[i]->SetParLimits(4,0.5,1.5);
     fitSingle[i]->SetParLimits(5,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
     fitSingle[i]->SetParLimits(6,10,130);

     //fitSingle[i] = new TF1(fitSingleName,"gaus(0)+gaus(3)",fitRange1,fitRange2);
     //fitSingle[i]->SetParNames("A1","E1","#sigma1","A2","E2","#sigma2");

     //fitSingle[i]->SetParameters(650,Peak1,75,650,Peak2,75);
     //fitSingle[i]->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
     //fitSingle[i]->SetParLimits(2,10,120);
     //fitSingle[i]->SetParLimits(4,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
     //fitSingle[i]->SetParLimits(5,10,120);

     fitSingle[i]->SetLineWidth(1);
     fitSingle[i]->SetLineColor(kBlue);

     hSingle[i]->Fit(fitSingleName,"r");
     //hMulti[i]->Fit(fitSingleName,"r");

     double par[7];
     fitSingle[i]->GetParameters(par);
     double par1 = fitSingle[i]->GetParameter(1);
     double par2 = fitSingle[i]->GetParameter(2);
     double par3 = fitSingle[i]->GetParameter(5);
     double par4 = fitSingle[i]->GetParameter(6);
     double *parErr = fitSingle[i]->GetParErrors();

     angle[i] = theta;
     res1[i] = fitSingle[i]->GetParameter(2)/fitSingle[i]->GetParameter(1);
     res2[i] = fitSingle[i]->GetParameter(6)/fitSingle[i]->GetParameter(5);
     resErr1[i] = TMath::Sqrt(1.0 / par1**2 * parErr[2]**2 + (par2 / par1**2)**2 * parErr[1]**2);
     resErr2[i] = TMath::Sqrt(1.0 / par3**2 * parErr[6]**2 + (par4 / par3**2)**2 * parErr[5]**2);
  }

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");
  t->Draw("csc:ecrec>>h2DMulti","nsite > 1");

  TGraphErrors *grSingleRes1 = new TGraphErrors(nrHistos,angle,res1,0,resErr1);
  grSingleRes1->SetMarkerStyle(20);
  grSingleRes1->SetMarkerSize(0.8);

  TF1 *fitRes1 = new TF1("fitRes1","pol2",0.15,0.75);
  fitRes1->SetLineWidth(1);
  fitRes1->SetLineColor(kBlue);

  grSingleRes1->Fit("fitRes1","r");

  double par1Res1 = fitRes1->GetParameter(2);
  double par2Res1 = fitRes1->GetParameter(1);

  TGraphErrors *grSingleRes2 = new TGraphErrors(nrHistos,angle,res2,0,resErr2);
  grSingleRes2->SetMarkerStyle(20);
  grSingleRes2->SetMarkerSize(0.8);

  TF1 *fitRes2 = new TF1("fitRes2","pol2",0.15,0.75);
  fitRes2->SetLineWidth(1);
  fitRes2->SetLineColor(kBlue);

  grSingleRes2->Fit("fitRes2","r");

  double par1Res2 = fitRes2->GetParameter(2);
  double par2Res2 = fitRes2->GetParameter(1);
  double *parResErr1 = fitRes1->GetParErrors();
  double *parResErr2 = fitRes2->GetParErrors();

  TF1 *fitRes21 = new TF1("fitRes21","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",0.15,0.75);
  fitRes21->SetParameter(2,-0.5*par2Res1/par1Res1);

  grSingleRes1->Fit("fitRes21","r");

  TF1 *fitRes22 = new TF1("fitRes22","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",0.15,0.75);
  fitRes22->SetParameter(2,-0.5*par2Res2/par1Res2);

  grSingleRes2->Fit("fitRes22","r");

  cout << fitRes21->GetParameter(2) << " +- " << fitRes21->GetParError(2) << "   Res = " << fitRes21->Eval(fitRes21->GetParameter(2)) << " +- " << fitRes21->Eval(fitRes21->GetParameter(2) + fitRes21->GetParError(2)) - fitRes21->Eval(fitRes21->GetParameter(2)) << endl;
  cout << fitRes22->GetParameter(2) << " +- " << fitRes22->GetParError(2) << "   Res = " << fitRes22->Eval(fitRes22->GetParameter(2)) << " +- " << fitRes22->Eval(fitRes22->GetParameter(2) + fitRes22->GetParError(2)) - fitRes22->Eval(fitRes22->GetParameter(2)) << endl;

  cout << "Best resolution1: " << -0.5*par2Res1/par1Res1 << endl;
  cout << "Best resolution2: " << -0.5*par2Res2/par1Res2 << endl;
  cout << "Best resolution1: " << -0.5*par2Res1/par1Res1 << " +- " << TMath::Sqrt(1.0/(2*par1Res1**2) * parResErr1[1]**2 + (par2Res1/(2*par1Res1**2))**2 * parResErr1[2]**2) << endl;
  cout << "Best resolution2: " << -0.5*par2Res2/par1Res2 << " +- " << TMath::Sqrt(1.0/(2*par1Res2**2) * parResErr2[1]**2 + (par2Res2/(2*par1Res2**2))**2 * parResErr2[2]**2) << endl;*/

  //double theta = (fitRes21->GetParameter(2) + fitRes22->GetParameter(2))/2;
  //double theta = 0.47033;
  //double theta = 0.43218;
  //double theta = 0.1956;
  //double theta = 0.2267;
  
  // Reprocessed data (01/30/12)
  //double theta = 0.1916;
  //double theta = 0.2164;
  
  // Reprocessed data (02/07/12)
  double theta = 0.1867;
  //double theta = 0.2084;

  TH1F *hMin = new TH1F("hMin","Co60 (single site, #theta = 0.1916)",200,0,2500);
  char cmd[100];
  //sprintf(cmd,"((csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f))*0.7192 + 49.73 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.6138 - 3.328 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.5791 + 1.925 >> hMin",theta,theta);
  
  // Reprocessed data (01/30/12)
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6117 + 2.443>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5826 + 8.801 >> hMin",theta,theta);
  
  // Reprocessed data (02/07/12)
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598>> hMin",theta,theta);

  t->Draw(cmd,"nsite == 1");

  double cscOff = 0;

/*  double fitRange1 = (3500-cscOff)*TMath::Sin(theta) + 900*TMath::Cos(theta);
  double fitRange2 = (6000-cscOff)*TMath::Sin(theta) + 1500*TMath::Cos(theta);
  double Peak1 = (4300-cscOff)*TMath::Sin(theta) + 1120*TMath::Cos(theta);
  double Peak2 = (4800-cscOff)*TMath::Sin(theta) + 1300*TMath::Cos(theta);
*/
  double fitRange1 = 1000;
  double fitRange2 = 1500;
  double Peak1 = 1173;
  double Peak2 = 1332;

  TF1 *fitMin = new TF1("fitMin",fitFunction,fitRange1,fitRange2,7);
  fitMin->SetParNames("A1","E1","#sigma1","R1","R2","E2","#sigma2");

  fitMin->SetParameters(3000,Peak1,75,0.6,0.8,Peak2,75);
  fitMin->SetParLimits(1,Peak1-Peak1/12.0,Peak1+Peak1/12.0);
  fitMin->SetParLimits(2,10,120);
  fitMin->SetParLimits(3,0.05,1.0);
  fitMin->SetParLimits(4,0.5,1.5);
  fitMin->SetParLimits(5,Peak2-Peak2/12.0,Peak2+Peak2/12.0);
  fitMin->SetParLimits(6,10,130);

  hMin->Fit("fitMin","r");

  double E1 = fitMin->GetParameter(1);
  double sigma1 = fitMin->GetParameter(2);
  double E2 = fitMin->GetParameter(5);
  double sigma2 = fitMin->GetParameter(6);
  double E1_err = fitMin->GetParError(1);
  double sigma1_err = fitMin->GetParError(2);
  double E2_err = fitMin->GetParError(5);
  double sigma2_err = fitMin->GetParError(6);

  cout << "Resolution from best theta: " << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2) << "    " << sigma2 / E2 << " +- " << TMath::Sqrt(sigma2_err**2 / E2**2 + (sigma2/E2**2)**2 * E2_err**2) << endl;

  /*TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");

  TCanvas *c2 = new TCanvas();
  h2DMulti->Draw("colz");

  TCanvas *c3 = new TCanvas();
  c3->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
     c3->cd(i+1);
     hSingle[i]->Draw();
     cout << res1[i] << endl;
     cout << res2[i] << endl;
  }

  TCanvas *c4 = new TCanvas();
  c4->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
     c4->cd(i+1);
     hMulti[i]->Draw();
  }

  TCanvas *c5 = new TCanvas();
  grSingleRes1->Draw("AZP");
  grSingleRes2->Draw("ZPsame");*/

  TCanvas *c6 = new TCanvas();
  hMin->Draw();

  return;
}
