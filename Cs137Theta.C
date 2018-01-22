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

void Cs137Theta()
{
  gStyle->SetPalette(1);

  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/V9/2450_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2469_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2473_noReclustering_Fiducial.root");

  /*const int nrHistos = 20;

  TH1F *hSingle[nrHistos];
  TH1F *hMulti[nrHistos];
  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,1500,200,0,5000);
  TH2F *h2DMulti = new TH2F("h2DMulti","scintillation vs ionization",200,0,1500,200,0,5000);

  TF1 *fitSingle[nrHistos];
  TF1 *fitMulti[nrHistos];

  double *angle = new double[nrHistos];
  double *res = new double[nrHistos];
  double *resErr = new double[nrHistos];

  double thetaRange1 = 0.0;
  double thetaRange2 = 0.8; //TMath::Pi()/9.0;
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

     hSingle[i] = new TH1F(hSingleName,hSingleTitle,200,0,2000);
     hMulti[i] = new TH1F(hMultiName,hMultiTitle,200,0,2000);

     char SingleCMD[200];
     char MultiCMD[200];
     sprintf(SingleCMD,"(csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f) >> %s",theta,theta,hSingleName);
     sprintf(MultiCMD,"csc*TMath::Sin(%.2f) + ecrec*TMath::Cos(%.2f) >> %s",theta,theta,hMultiName);

     t->Draw(SingleCMD,"nsite == 1");
     t->Draw(MultiCMD,"nsite > 1");

     char fitSingleName[100];
     sprintf(fitSingleName,"fit_%i",i);

     double cscOff = 0;

     double fitRange1 = (500-cscOff)*TMath::Sin(theta) + 500*TMath::Cos(theta);
     double fitRange2 = (800-cscOff)*TMath::Sin(theta) + 800*TMath::Cos(theta);
     double Peak = (660-cscOff)*TMath::Sin(theta) + 660*TMath::Cos(theta);

     fitSingle[i] = new TF1(fitSingleName,fitFunction,fitRange1,fitRange2,4);

     fitSingle[i]->SetParNames("A1","E1","#sigma","A2");

     fitSingle[i]->SetParameters(800,Peak,100,0.8);

     fitSingle[i]->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
     fitSingle[i]->SetParLimits(2,10,200);
     fitSingle[i]->SetParLimits(3,0.05,2.0);

     fitSingle[i]->SetLineWidth(1);
     fitSingle[i]->SetLineColor(kBlue);

     hSingle[i]->Fit(fitSingleName,"r");
     //hMulti[i]->Fit(fitSingleName,"r");

     double par[4];
     fitSingle[i]->GetParameters(par);
     double par1 = fitSingle[i]->GetParameter(1);
     double par2 = fitSingle[i]->GetParameter(2);
     double *parErr = fitSingle[i]->GetParErrors();

     angle[i] = theta;
     res[i] = fitSingle[i]->GetParameter(2)/fitSingle[i]->GetParameter(1);
     resErr[i] = TMath::Sqrt(1.0 / par1**2 * parErr[2]**2 + (par[2] / par1**2)**2 * parErr[1]**2);
  }

  t->Draw("csc:ecrec>>h2DSingle","nsite == 1");
  t->Draw("csc:ecrec>>h2DMulti","nsite > 1");

  TGraphErrors *grSingleRes = new TGraphErrors(nrHistos,angle,res,0,resErr);
  //TGraph *grSingleRes = new TGraph(nrHistos,angle,res);
  grSingleRes->SetMarkerStyle(20);
  grSingleRes->SetMarkerSize(0.8);

  TF1 *fitRes = new TF1("fitRes","pol2",0.06,0.7);
  fitRes->SetLineWidth(1);
  fitRes->SetLineColor(kBlue);

  grSingleRes->Fit("fitRes","r");

  double parRes1 = fitRes->GetParameter(2);
  double parRes2 = fitRes->GetParameter(1);
  double *parResErr = fitRes->GetParErrors();

  TF1 *fitRes2 = new TF1("fitRes2","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",0.06,0.7);
  fitRes2->SetParameter(2,-0.5*parRes2/parRes1);

  grSingleRes->Fit("fitRes2","r");
  cout << fitRes2->GetParameter(2) << " +- " << fitRes2->GetParError(2) << "   Res = " << fitRes2->Eval(fitRes2->GetParameter(2)) << " +- " << fitRes2->Eval(fitRes2->GetParameter(2) + fitRes2->GetParError(2)) - fitRes2->Eval(fitRes2->GetParameter(2)) << endl;

  cout << "Best resolution: " << -0.5*parRes2/parRes1 << " +- " << TMath::Sqrt(1.0/(2*parRes1**2) * parResErr[1]**2 + (parRes2/(2*parRes1**2))**2 * parResErr[2]**2) << endl;*/

  //double theta = (fitRes21->GetParameter(2) + fitRes22->GetParameter(2))/2;
  //double theta = 0.37082;
  //double theta = 0.1956;
  //double theta = 0.2267;
  
  // Reprocessed data (01/30/12)
  //double theta = 0.1916;
  //double theta = 0.2164;
  
  // Reprocessed data (02/07/12)
  double theta = 0.1867;
  //double theta = 0.2084;

  TH1F *hMin = new TH1F("hMin","Cs137 (single site, #theta = 0.1916)",200,0,3000);
  char cmd[100];
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.6138 - 3.328 >> hMin",theta,theta);
  //sprintf(cmd,"((csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f))*0.7192 + 49.73 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.5791 + 1.925 >> hMin",theta,theta);
  
  // Reprocessed data (01/30/12)
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6117 + 2.443>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5826 + 8.801 >> hMin",theta,theta);
  
  // Reprocessed data (02/07/12)
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598>> hMin",theta,theta);

  t->Draw(cmd,"nsite == 1");

  double cscOff = 0;

/*  double fitRange1 = (1800-cscOff)*TMath::Sin(theta) + 450*TMath::Cos(theta);
  double fitRange2 = (3300-cscOff)*TMath::Sin(theta) + 750*TMath::Cos(theta);
  double Peak = (2500-cscOff)*TMath::Sin(theta) + 600*TMath::Cos(theta);
*/
  double fitRange1 = 500;
  double fitRange2 = 800;
  double Peak = 662;

  TF1 *fitMin = new TF1("fitMin",fitFunction,fitRange1,fitRange2,4);

  fitMin->SetParNames("A1","E1","#sigma","A2");

  fitMin->SetParameters(10000,Peak,100,0.2);

  fitMin->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
  fitMin->SetParLimits(2,10,200);
  fitMin->SetParLimits(3,0.05,2.0);

  hMin->Fit("fitMin","r");

  double E1 = fitMin->GetParameter(1);
  double sigma1 = fitMin->GetParameter(2);
  double E1_err = fitMin->GetParError(1);
  double sigma1_err = fitMin->GetParError(2);

  cout << "Resolution from best theta: " << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2)  << endl;

  /*TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");

  TCanvas *c2 = new TCanvas();
  h2DMulti->Draw("colz");

  TCanvas *c3 = new TCanvas();
  c3->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
     c3->cd(i+1);
     hSingle[i]->Draw();
     cout << res[i] << endl;
  }

  TCanvas *c4 = new TCanvas();
  c4->Divide(5,4);
  for (int i = 0; i < nrHistos; i++) {
     c4->cd(i+1);
     hMulti[i]->Draw();
  }

  TCanvas *c5 = new TCanvas();
  grSingleRes->Draw("AZP");*/

  TCanvas *c6 = new TCanvas();
  hMin->Draw();

  return;
}
