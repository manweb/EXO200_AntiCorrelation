using namespace RooFit;

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

void Th228RotatedSpectrum()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  TChain *t = new TChain("t","tree");
  t->Add("../../EnergyCalibration/analysis/V11/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V11/2448_noReclustering_Fiducial.root");
  
  TChain *t2 = new TChain("t","tree");
  t2->Add("../../EnergyCalibration/analysis/V12/2424_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2426_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2431_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2432_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2433_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2434_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2447_noReclustering_Fiducial.root");
  t2->Add("../../EnergyCalibration/analysis/V12/2448_noReclustering_Fiducial.root");
  
  /*RooRealVar e_charge("epcrec","epcrec",0,3500);
  RooRealVar e_scint("csc","csc",0,50000);
  RooRealVar mean_ss("mean_ss","mean_ss",2615,2400,2800);
  RooRealVar sigma_ss("sigma_ss","sigma_ss",70,20,200);
  RooRealVar mean_ms("mean_ms","mean_ms",2615,2400,2800);
  RooRealVar sigma_ms("sigma_ms","sigma_ms",70,20,200);
  RooRealVar nWires("nWires","nWires",0,10);
  RooRealVar numSite("nsite","nsite",0,100);
  
  RooRealVar theta_ss("theta_ss","theta_ss",0.1814);
  RooRealVar theta_ms("theta_ms","theta_ms",0.2036);
  
  RooRealVar p0_ss_mean("p0_ss_mean","p0_ss_mean",-1.221);
  RooRealVar p1_ss_mean("p1_ss_mean","p1_ss_mean",0.6188);
  RooRealVar p0_ss_err("p0_ss_err","p0_ss_err",4.803);
  RooRealVar p1_ss_err("p1_ss_err","p1_ss_err",0.002);
  RooRealVar p0_ss("p0_ss","p0_ss",p0_ss_mean.getVal(),p0_ss_mean.getVal()-5*p0_ss_err.getVal(),p0_ss_mean.getVal()+5*p0_ss_err.getVal());
  RooRealVar p1_ss("p1_ss","p1_ss",p1_ss_mean.getVal(),p1_ss_mean.getVal()-5*p1_ss_err.getVal(),p1_ss_mean.getVal()+5*p1_ss_err.getVal());
  
  RooRealVar p0_ms_mean("p0_ms_mean","p0_ms_mean",-0.4439);
  RooRealVar p1_ms_mean("p1_ms_mean","p1_ms_mean",0.5958);
  RooRealVar p0_ms_err("p0_ms_err","p0_ms_err",5.198);
  RooRealVar p1_ms_err("p1_ms_err","p1_ms_err",0.002);
  RooRealVar p0_ms("p0_ms","p0_ms",p0_ms_mean.getVal(),p0_ms_mean.getVal()-5*p0_ms_err.getVal(),p0_ms_mean.getVal()+5*p0_ms_err.getVal());
  RooRealVar p1_ms("p1_ms","p1_ms",p1_ms_mean.getVal(),p1_ms_mean.getVal()-5*p1_ms_err.getVal(),p1_ms_mean.getVal()+5*p1_ms_err.getVal());
  
  RooFormulaVar energy_ss("energy_ss","energy_ss","(@0*sin(@2) + @1*cos(@2))*@3 + @4",RooArgList(e_scint,e_charge,theta_ss,p1_ss,p0_ss));
  RooFormulaVar energy_ms("energy_ms","energy_ms","(@0*sin(@2) + @1*cos(@2))*@3 + @4",RooArgList(e_scint,e_charge,theta_ms,p1_ms,p0_ms));
  
  RooDataSet *data = new RooDataSet("data","data",t,RooArgSet(e_charge,e_scint,nWires,numSite));
  RooRealVar *e_ss = data->addColumn(energy_ss);
  RooRealVar *e_ms = data->addColumn(energy_ms);
  
  RooRealVar ratio_ss("ratio_ss","ratio_ss",0.2,0.05,1.0);
  RooRealVar ratio_ms("ratio_ms","ratio_ms",0.2,0.05,1.0);
  
  RooGaussian gaus_ss("gaus_ss","gaus_ss",*e_ss,mean_ss,sigma_ss);
  RooGaussian gaus_ms("gaus_ms","gaus_ms",*e_ms,mean_ms,sigma_ms);
  
  RooGenericPdf erfc_ss("erfc_ss","erfc_ss","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*e_ss,mean_ss,sigma_ss));
  RooGenericPdf erfc_ms("erfc_ms","erfc_ms","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*e_ms,mean_ms,sigma_ms));
  
  RooAddPdf fit_pdf_ss("fit_pdf_ss","fit_pdf_ss",RooArgList(gaus_ss,erfc_ss),ratio_ss);
  RooAddPdf fit_pdf_ms("fit_pdf_ms","fit_pdf_ms",RooArgList(gaus_ms,erfc_ms),ratio_ms);
  
  e_ss->setRange(2450,2800);
  e_ms->setRange(2450,2800);
  e_ss->setBins(200);
  e_ms->setBins(200);
  
  fit_pdf_ss.fitTo(*data,Cut("nsite == 1 && nWires <= 2"),Range(2450,2800));
  fit_pdf_ms.fitTo(*data,Cut("nsite > 2 || nsite == 1 && nWires > 2"),Range(2450,2800));
  
  RooPlot *frame_ss = e_ss->frame(Range(0,3500));
  data->plotOn(frame_ss,Cut("nsite == 1 && nWires <= 2"),MarkerSize(0.8),DrawOption("ZP"));
  fit_pdf_ss.plotOn(frame_ss,LineWidth(2));
  fit_pdf_ss.plotOn(frame_ss,Components("gaus_ss"),LineStyle(2),LineWidth(2));
  fit_pdf_ss.plotOn(frame_ss,Components("erfc_ss"),LineStyle(2),LineWidth(2));
  
  RooPlot *frame_ms = e_ms->frame(Range(0,3500));
  data->plotOn(frame_ms,Cut("nsite > 2 || nsite == 1 && nWires > 2"),MarkerSize(0.8),DrawOption("ZP"));
  fit_pdf_ms.plotOn(frame_ms,LineWidth(2));
  fit_pdf_ms.plotOn(frame_ms,Components("gaus_ms"),LineStyle(2),LineWidth(2));
  fit_pdf_ms.plotOn(frame_ms,Components("erfc_ms"),LineStyle(2),LineWidth(2));
  
  cout << "Resolution (single site): " << sigma_ss.getVal() / mean_ss.getVal() << endl;
  cout << "Resolution (multi site): " << sigma_ms.getVal() / mean_ms.getVal() << endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  frame_ss->Draw();
  c1->cd(2);
  frame_ms->Draw();*/
  
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
  
  TH1F *h_ss = new TH1F("h_ss",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ms = new TH1F("h_ms",Form("Th228 (multi site, #theta = %.4f)",theta_ms),200,2000,3000);
  char *cmd_ms = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ms",theta_ms,theta_ms,p1_ms,p0_ms);
  
  TH1F *h_ss_1 = new TH1F("h_ss_1",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_1 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_1",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ss_2 = new TH1F("h_ss_2",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f*1.00683 + %.4f >> h_ss_2",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ss_3 = new TH1F("h_ss_3",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_3 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_3",theta_ss,theta_ss,p1_ss,p0_ss);
  
  TH1F *h_ss_2 = new TH1F("h_ss_2",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_2",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ms_2 = new TH1F("h_ms_2",Form("Th228 (multi site, #theta = %.4f)",theta_ms),200,2000,3000);
  char *cmd_ms_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ms_2",theta_ms,theta_ms,p1_ms,p0_ms);
  
  TH1F *h_ss_1_2 = new TH1F("h_ss_1_2",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_1_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f + %.4f >> h_ss_1_2",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  TH1F *h_ss_2_2 = new TH1F("h_ss_2_2",Form("Th228 (single site, #theta = %.4f)",theta_ss),200,2000,3000);
  char *cmd_ss_2_2 = Form("(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*%.4f*1.00683 + %.4f >> h_ss_2_2",theta_ss,theta_ss,p1_1u_ss,p0_1u_ss);
  
  h_ss_1->SetLineColor(kRed);
  h_ss_2->SetLineColor(kGreen);
  
  h_ss_1_2->SetLineColor(kRed);
  h_ss_2_2->SetLineColor(kGreen);
  
  t->Draw(cmd_ss,"nsite == 1 && nWires <= 2");
  t->Draw(cmd_ms,"nsite > 1 || nsite == 1 && nWires > 2");
  t->Draw(cmd_ss_1,"nsite == 1 && nWires == 1");
  t->Draw(cmd_ss_2,"nsite == 1 && nWires == 2");
  
  t2->Draw(cmd_ss_2,"nsite == 1 && nWires <= 2");
  t2->Draw(cmd_ms_2,"nsite > 1 || nsite == 1 && nWires > 2");
  t2->Draw(cmd_ss_1_2,"nsite == 1 && nWires == 1");
  t2->Draw(cmd_ss_2_2,"nsite == 1 && nWires == 2");
  
  TH1F *hSum = h_ss_1->Clone("hSum");
  hSum->Add(h_ss_2);
  
  hSum->SetLineColor(kBlack);
  
  double fitRange1_ss = 2480;
  double fitRange2_ss = 2720;
  double fitRange1_ms = 2480;
  double fitRange2_ms = 2720;
  double Peak = 2615;
  
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
  
  hSum->Fit("fit_ss","r");
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
  //h_ss->Draw("EZP");
  hSum->Draw("EZP");
  h_ss_1->Draw("EZPsame");
  h_ss_2->Draw("EZPsame");
  h_ss->Draw("same");
  
  TCanvas *c2 = new TCanvas();
  h_ms->Draw("EZP");
}