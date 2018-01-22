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

void Th228Theta()
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);

  TChain *t = new TChain("t","tree");
  /*t->Add("../../EnergyCalibration/analysis/V7/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V7/2448_noReclustering_Fiducial.root");*/
  
  // Reprocessed data (01/30/12)
  /*t->Add("../../EnergyCalibration/analysis/V8/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2448_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2817_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2837_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2848_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2865_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2867_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2883_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2897_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2904_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2938_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2943_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2947_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2966_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2981_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/2991_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3007_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3019_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3028_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3034_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3053_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3074_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3091_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3099_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V8/3109_noReclustering_Fiducial.root");*/
  
  // Reprocessed data (02/07/12)
  /*t->Add("../../EnergyCalibration/analysis/V9/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2448_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2817_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2837_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2848_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2865_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2867_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2883_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2897_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2904_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2938_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2943_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2947_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2966_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2981_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/2991_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3007_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3019_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3028_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3034_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3053_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3074_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3091_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3099_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/V9/3109_noReclustering_Fiducial.root");*/
  
  // Reprocessed data, alpha branch (02/18/12)
  /*t->Add("../../EnergyCalibration/analysis/alpha/V1/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V1/2448_noReclustering_Fiducial.root");*/
  
  // Reprocessed data, alpha branch (02/22/12)
  /*t->Add("../../EnergyCalibration/analysis/alpha/V2/2424_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2426_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2431_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2432_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2433_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2434_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2447_noReclustering_Fiducial.root");
  t->Add("../../EnergyCalibration/analysis/alpha/V2/2448_noReclustering_Fiducial.root");*/
  
  // Reprocessed data, alpha branch (03/21/12)
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2424_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2426_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2431_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2432_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2433_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2434_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2447_noReclustering_Fiducial.root");
  t->Add("/nfs/slac/g/exo/maweber/EXO200Analysis/EnergyCalibration/analysis/V12/2448_noReclustering_Fiducial.root");
  
  cout << t->GetEntries() << endl;

  const int nrHistos = 20;

  TH1F *hSingle[nrHistos];
  TH1F *hMulti[nrHistos];
  TH2F *h2DSingle = new TH2F("h2DSingle","scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h2DMulti = new TH2F("h2DMulti","scintillation vs ionization",200,0,3500,200,0,15000);

  TF1 *fitSingle[nrHistos];
  TF1 *fitMulti[nrHistos];

  double *angle_ss = new double[nrHistos];
  double *res_ss = new double[nrHistos];
  double *resErr_ss = new double[nrHistos];
  
  double *angle_ms = new double[nrHistos];
  double *res_ms = new double[nrHistos];
  double *resErr_ms = new double[nrHistos];

  double thetaRange1_ss = 0.10; // 0.13
  double thetaRange2_ss = 0.25; //0.26 //TMath::Pi()/9.0;
  double dTheta_ss = thetaRange2_ss - thetaRange1_ss;
  double stp_ss = dTheta_ss/nrHistos;
  
  double thetaRange1_ms = 0.10; //0.15
  double thetaRange2_ms = 0.25; //0.28 //TMath::Pi()/9.0;
  double dTheta_ms = thetaRange2_ms - thetaRange1_ms;
  double stp_ms = dTheta_ms/nrHistos;
  
  double resFitRange1_ss = 0.14; //0.14
  double resFitRange2_ss = 0.19; //0.22
  
  double resFitRange1_ms = 0.15; //0.16
  double resFitRange2_ms = 0.20; //0.245
  
// #### RooFit ##############################################################################
  
  RooRealVar energy("epcrec","epcrec",0.0,1e6);
  RooRealVar scint("csc","csc",0.0,1e6);
  RooRealVar numSite("nsite","nsite",0,1000);
  RooRealVar nWires("nWires","nWires",-1000,1000);
  RooRealVar sigma_ss("sigma_ss","sigma_ss",70,50,200);
  RooRealVar mean_ss("mean_ss","mean_ss",3000,6000);
  RooRealVar sigma_ms("sigma_ms","sigma_ms",70,50,200);
  RooRealVar mean_ms("mean_ms","mean_ms",3000,6000);
  
  RooRealVar theta_ss("theta_ss","theta_ss",0.0);
  RooRealVar theta_ms("theta_ms","theta_ms",0.0);
  
  RooRealVar ratio_ss("ratio_ss","ratio_ss",0.2,0.05,1.0);
  RooRealVar ratio_ms("ratio_ms","ratio_ms",0.2,0.05,1.0);
  
  //RooAddPdf fit_pdf_ss("fit_pdf_ss","fit_pdf_ss",RooArgList(gaus_ss,erfc_ss),ratio_ss);
  //RooAddPdf fit_pdf_ms("fit_pdf_ms","fit_pdf_ms",RooArgList(gaus_ms,erfc_ms),ratio_ms);
  
  RooDataHist *h_ss;
  RooDataHist *h_ms;
  
  RooDataSet *data = new RooDataSet("data","data",t,RooArgSet(energy,scint,numSite,nWires));
  
  RooFormulaVar *rotated_ss = new RooFormulaVar("rotated_ss","rotated_ss","@0*TMath::Sin(@1) + @2*TMath::Cos(@1)",RooArgList(scint,theta_ss,energy));
  RooFormulaVar *rotated_ms = new RooFormulaVar("rotated_ms","rotated_ms","@0*TMath::Sin(@1) + @2*TMath::Cos(@1)",RooArgList(scint,theta_ms,energy));
  
  //RooRealVar *e_rotated_ss = (RooRealVar*)data->addColumn(*rotated_ss);
  //RooRealVar *e_rotated_ms = (RooRealVar*)data->addColumn(*rotated_ms);
  
  RooDataSet *data_ss = (RooDataSet*)data->reduce("nsite == 1 && nWires <= 2");
  RooDataSet *data_ms = (RooDataSet*)data->reduce("nsite > 1 || nsite == 1 && nWires > 2");
  
  data_ss->SetName("data_ss");
  data_ms->SetName("data_ms");
  
  /*RooGaussian gaus_ss("gaus_ss","gaus_ss",*rotated_ss,mean_ss,sigma_ss);
  RooGaussian gaus_ms("gaus_ms","gaus_ms",*rotated_ms,mean_ms,sigma_ms);
  
  RooGenericPdf erfc_ss("erfc_ss","erfc_ss","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*rotated_ss,mean_ss,sigma_ss));
  RooGenericPdf erfc_ms("erfc_ms","erfc_ms","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*rotated_ms,mean_ms,sigma_ms));*/
  
  RooPlot *frame_ss[nrHistos];
  RooPlot *frame_ms[nrHistos];
  
// ##########################################################################################

  for (int i = 0; i < nrHistos; i++) {
    double theta_ss_c = thetaRange1_ss + stp_ss*i;
    double theta_ms_c = thetaRange1_ms + stp_ms*i;
    
    /*char *hSingleName = Form("hSingle_%i",i);
    char *hSingleTitle = Form("%.4f rotation angle",theta_ss);
    char *hMultiName = Form("hMulti_%i",i);
    char *hMultiTitle = Form("%.4f rotation angle",theta_ms);
    
    hSingle[i] = new TH1F(hSingleName,hSingleTitle,100,3000,5000);
    hMulti[i] = new TH1F(hMultiName,hMultiTitle,100,3000,5000);

    char *SingleCMD = Form("csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f) >> %s",theta_ss,theta_ss,hSingleName);
    char *MultiCMD = Form("csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f) >> %s",theta_ms,theta_ms,hMultiName);
    
    t->Draw(SingleCMD,"nsite == 1 && nWires <= 2","goff");
    t->Draw(MultiCMD,"nsite > 1 || nsite == 1 && nWires > 2","goff");*/
    
// #### RooFit #############################################################################
    
    //h_ss = new RooDataHist("h_ss","h_ss",energy,hSingle[i]);
    //h_ms = new RooDataHist("h_ms","h_ms",energy,hMulti[i]);
    
// #########################################################################################
    
    //char *fitSingleName = Form("fit_%i_ss",i);
    //char *fitMultiName = Form("fit_%i_ms",i);
    
    double fitRange1_ss = 7800*TMath::Sin(theta_ss_c) + 2480*TMath::Cos(theta_ss_c);
    double fitRange2_ss = 9100*TMath::Sin(theta_ss_c) + 2700*TMath::Cos(theta_ss_c);
    double Peak_ss = 8500*TMath::Sin(theta_ss_c) +2610*TMath::Cos(theta_ss_c);
    
    double fitRange1_ms = 8000*TMath::Sin(theta_ms_c) + 2580*TMath::Cos(theta_ms_c);
    double fitRange2_ms = 9280*TMath::Sin(theta_ms_c) + 2740*TMath::Cos(theta_ms_c);
    double Peak_ms = 8500*TMath::Sin(theta_ms_c) +2660*TMath::Cos(theta_ms_c);
    
    /*fitSingle[i] = new TF1(fitSingleName,fitFunction,fitRange1_ss,fitRange2_ss,4);
    fitSingle[i]->SetParNames("A1","E1","#sigma","A2");
    fitSingle[i]->SetParameters(10000,Peak_ss,100,0.1);
    fitSingle[i]->SetParLimits(1,Peak_ss-Peak_ss/12.0,Peak_ss+Peak_ss/12.0);
    fitSingle[i]->SetParLimits(2,50,200);
    fitSingle[i]->SetParLimits(3,0.01,0.2);
    fitSingle[i]->SetLineWidth(2);
    fitSingle[i]->SetLineColor(kBlue);
    
    fitMulti[i] = new TF1(fitMultiName,fitFunction,fitRange1_ms,fitRange2_ms,4);
    fitMulti[i]->SetParNames("A1","E1","#sigma","A2");
    fitMulti[i]->SetParameters(10000,Peak_ms,100,0.07);
    fitMulti[i]->SetParLimits(1,Peak_ms-Peak_ms/12.0,Peak_ms+Peak_ms/12.0);
    fitMulti[i]->SetParLimits(2,50,200);
    fitMulti[i]->SetParLimits(3,0.05,0.1);
    fitMulti[i]->SetLineWidth(2);
    fitMulti[i]->SetLineColor(kBlue);
    
    hSingle[i]->Fit(fitSingleName,"r");
    hMulti[i]->Fit(fitMultiName,"r");*/
    
// #### RooFit #############################################################################
    
    mean_ss.setVal(Peak_ss);
    mean_ms.setVal(Peak_ms);
    
    theta_ss.setVal(theta_ss_c);
    theta_ms.setVal(theta_ms_c);
    
    RooDataSet *data_TMP_ss = (RooDataSet*)data_ss->Clone("data_TMP_ss");
    RooDataSet *data_TMP_ms = (RooDataSet*)data_ms->Clone("data_TMP_ms");
    
    RooRealVar *e_rotated_ss = (RooRealVar*)data_TMP_ss->addColumn(*rotated_ss);
    RooRealVar *e_rotated_ms = (RooRealVar*)data_TMP_ms->addColumn(*rotated_ms);
    
    RooGaussian gaus_ss("gaus_ss","gaus_ss",*e_rotated_ss,mean_ss,sigma_ss);
    RooGaussian gaus_ms("gaus_ms","gaus_ms",*e_rotated_ms,mean_ms,sigma_ms);
    
    RooGenericPdf erfc_ss("erfc_ss","erfc_ss","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*e_rotated_ss,mean_ss,sigma_ss));
    RooGenericPdf erfc_ms("erfc_ms","erfc_ms","0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))",RooArgList(*e_rotated_ms,mean_ms,sigma_ms));
    
    RooAddPdf fit_pdf_ss(Form("fit_pdf_ss_%i",i),Form("fit_pdf_ss_%i",i),RooArgList(gaus_ss,erfc_ss),ratio_ss);
    RooAddPdf fit_pdf_ms(Form("fit_pdf_ms_%i",i),Form("fit_pdf_ms_%i",i),RooArgList(gaus_ms,erfc_ms),ratio_ms);
    
    //fit_pdf_ss.setRange(fitRange1_ss,fitRange2_ss);
    //fit_pdf_ms.setRange(fitRange1_ms,fitRange2_ms);
    
    //energy.setRange(fitRange1_ss,fitRange2_ss);
    //energy.setBins(200);
    
    //RooRealVar *e_rotated_ss = (RooRealVar*)data_ss->addColumn(*rotated_ss);
    //RooRealVar *e_rotated_ms = (RooRealVar*)data_ms->addColumn(*rotated_ms);
    
    e_rotated_ss->setRange(fitRange1_ss,fitRange2_ss);
    e_rotated_ms->setRange(fitRange1_ms,fitRange2_ms);
    e_rotated_ss->setBins(200);
    e_rotated_ms->setBins(200);
    
    fit_pdf_ss.fitTo(*data_TMP_ss,Range(fitRange1_ss,fitRange2_ss),NumCPU(8));
    fit_pdf_ms.fitTo(*data_TMP_ms,Range(fitRange1_ms,fitRange2_ms),NumCPU(8));
    
    frame_ss[i] = e_rotated_ss->frame(Bins(200),Range(3000,6000));
    data_TMP_ss->plotOn(frame_ss[i],DrawOption("ZP"),MarkerStyle(20),MarkerSize(0.5));
    fit_pdf_ss.plotOn(frame_ss[i]);
    fit_pdf_ss.plotOn(frame_ss[i],Components("gaus_ss"),LineWidth(2),LineStyle(2));
    fit_pdf_ss.plotOn(frame_ss[i],Components("erfc_ss"),LineWidth(2),LineStyle(2));
    fit_pdf_ss.paramOn(frame_ss[i]);
    frame_ss[i]->SetTitle(Form("single site, #Theta = %.4f",theta_ss_c));
    
    frame_ms[i] = e_rotated_ms->frame(Bins(200),Range(3000,6000));
    data_TMP_ms->plotOn(frame_ms[i],DrawOption("ZP"),MarkerStyle(20),MarkerSize(0.5));
    fit_pdf_ms.plotOn(frame_ms[i]);
    fit_pdf_ms.plotOn(frame_ms[i],Components("gaus_ms"),LineWidth(2),LineStyle(2));
    fit_pdf_ms.plotOn(frame_ms[i],Components("erfc_ms"),LineWidth(2),LineStyle(2));
    fit_pdf_ms.paramOn(frame_ms[i]);
    frame_ms[i]->SetTitle(Form("multi site, #Theta = %.4f",theta_ms_c));
    
// #########################################################################################
    
    //double par_ss[4];
    //fitSingle[i]->GetParameters(par_ss);
    //double par1_ss = fitSingle[i]->GetParameter(1);
    //double par2_ss = fitSingle[i]->GetParameter(2);
    double par1_ss = mean_ss.getVal();
    double par2_ss = sigma_ss.getVal();
    double par1Err_ss = mean_ss.getError();
    double par2Err_ss = sigma_ss.getError();
    //double *parErr_ss = fitSingle[i]->GetParErrors();
    //double par1Err_ss = parErr_ss[1];
    //double par2Err_ss = parErr_ss[2];
    
    angle_ss[i] = theta_ss_c;
    //res_ss[i] = fitSingle[i]->GetParameter(2)/fitSingle[i]->GetParameter(1);
    //resErr_ss[i] = TMath::Sqrt(1.0 / par1_ss**2 * parErr_ss[2]**2 + (par_ss[2] / par1_ss**2)**2 * parErr_ss[1]**2);
    
    res_ss[i] = par2_ss/par1_ss;
    resErr_ss[i] = TMath::Sqrt(1.0 / par1_ss**2 * par2Err_ss**2 + (par2_ss / par1_ss**2)**2 * par1Err_ss**2);
    
    cout << "theta = " << angle_ss[i] << "\t res = " << res_ss[i] << endl;
    
    //double par_ms[4];
    //fitMulti[i]->GetParameters(par_ms);
    //double par1_ms = fitMulti[i]->GetParameter(1);
    //double par2_ms = fitMulti[i]->GetParameter(2);
    double par1_ms = mean_ss.getVal();
    double par2_ms = sigma_ss.getVal();
    double par1Err_ms = mean_ss.getError();
    double par2Err_ms = sigma_ss.getError();
    //double *parErr_ms = fitMulti[i]->GetParErrors();
    //double par1Err_ms = parErr_ms[1];
    //double par2Err_ms = parErr_ms[2];
    
    angle_ms[i] = theta_ms_c;
    //res_ms[i] = fitMulti[i]->GetParameter(2)/fitMulti[i]->GetParameter(1);
    //resErr_ms[i] = TMath::Sqrt(1.0 / par1_ms**2 * parErr_ms[2]**2 + (par_ms[2] / par1_ms**2)**2 * parErr_ms[1]**2);
    
    res_ms[i] = par2_ms/par1_ms;
    resErr_ms[i] = TMath::Sqrt(1.0 / par1_ms**2 * par2Err_ms**2 + (par2_ms / par1_ms**2)**2 * par1Err_ms**2);
    
    cout << "theta = " << angle_ms[i] << "\t res = " << res_ms[i] << endl;
  }

  t->Draw("csc:epcrec>>h2DSingle","nsite == 1 && nWires <= 2","goff");
  t->Draw("csc:epcrec>>h2DMulti","nsite > 1 || nsite == 1 && nWires > 2","goff");
  
// #### Single site minimization ############################################################################
  TGraphErrors *grSingleRes = new TGraphErrors(nrHistos,angle_ss,res_ss,0,resErr_ss);
  grSingleRes->SetMarkerStyle(20);
  grSingleRes->SetMarkerSize(0.8);

  TF1 *fitRes_ss = new TF1("fitRes_ss","pol2",resFitRange1_ss,resFitRange2_ss);
  fitRes_ss->SetLineWidth(2);
  fitRes_ss->SetLineColor(kRed);

  grSingleRes->Fit("fitRes_ss","rn");

  double parRes1_ss = fitRes_ss->GetParameter(2);
  double parRes2_ss = fitRes_ss->GetParameter(1);
  double *parResErr_ss = fitRes_ss->GetParErrors();

  TF1 *fitRes2_ss = new TF1("fitRes2_ss","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",resFitRange1_ss,resFitRange2_ss);
  fitRes2_ss->SetParameter(2,-0.5*parRes2_ss/parRes1_ss);
  fitRes2_ss->SetLineWidth(2);
  fitRes2_ss->SetLineColor(kBlue);

  grSingleRes->Fit("fitRes2_ss","r");
  cout << fitRes2_ss->GetParameter(2) << " +- " << fitRes2_ss->GetParError(2) << "   Res = " << fitRes2_ss->Eval(fitRes2_ss->GetParameter(2)) << " +- " << fitRes2_ss->Eval(fitRes2_ss->GetParameter(2) + fitRes2_ss->GetParError(2)) - fitRes2_ss->Eval(fitRes2_ss->GetParameter(2)) << endl;
// ##########################################################################################################
  
// #### Multi site minimization ############################################################################
  TGraphErrors *grMultiRes = new TGraphErrors(nrHistos,angle_ms,res_ms,0,resErr_ms);
  grMultiRes->SetMarkerStyle(20);
  grMultiRes->SetMarkerSize(0.8);
  
  TF1 *fitRes_ms = new TF1("fitRes_ms","pol2",resFitRange1_ms,resFitRange2_ms);
  fitRes_ms->SetLineWidth(2);
  fitRes_ms->SetLineColor(kRed);
  
  grMultiRes->Fit("fitRes_ms","rn");
  
  double parRes1_ms = fitRes_ms->GetParameter(2);
  double parRes2_ms = fitRes_ms->GetParameter(1);
  double *parResErr_ms = fitRes_ms->GetParErrors();
  
  TF1 *fitRes2_ms = new TF1("fitRes2_ms","[0]*x*x - 2*[0]*[2]*x + [0]*[2]*[2] + [1]",resFitRange1_ms,resFitRange2_ms);
  fitRes2_ms->SetParameter(2,-0.5*parRes2_ms/parRes1_ms);
  fitRes2_ms->SetLineWidth(2);
  fitRes2_ms->SetLineColor(kBlue);
  
  grMultiRes->Fit("fitRes2_ms","r");
  cout << fitRes2_ms->GetParameter(2) << " +- " << fitRes2_ms->GetParError(2) << "   Res = " << fitRes2_ms->Eval(fitRes2_ms->GetParameter(2)) << " +- " << fitRes2_ms->Eval(fitRes2_ms->GetParameter(2) + fitRes2_ms->GetParError(2)) - fitRes2_ms->Eval(fitRes2_ms->GetParameter(2)) << endl;
// ##########################################################################################################

  //cout << "Best resolution: " << -0.5*parRes2/parRes1 << " +- " << TMath::Sqrt(1.0/(2*parRes1**2) * parResErr[1]**2 + (parRes2/(2*parRes1**2))**2 * parResErr[2]**2) << endl;

  //double theta = (fitRes21->GetParameter(2) + fitRes22->GetParameter(2))/2;
  //double theta = 0.5369;
  //double theta = 0.1956;
  //double theta = 0.2267;
  
  // Reprocessed data (01/30/12)
  //double theta = 0.1916;
  //double theta = 0.2164;
  
  // Reprocessed data (02/07/12)
  /*double theta = 0.1867;
  double theta = 0.2084;

  TH1F *hMin = new TH1F("hMin","Th228 (single site, #theta = 0.1867)",200,0,3500);
  char cmd[100];
  //sprintf(cmd,"((csc*0.3148-119.5)*TMath::Sin(%.2f) + (ecrec*0.9647+67.4)*TMath::Cos(%.2f))*0.7192 + 49.73 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.6138 - 3.328 >> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + ecrec*TMath::Cos(%.4f))*0.5791 + 1.925 >> hMin",theta,theta);
  
  // Reprocessed data (01/30/12)
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6117 + 2.443>> hMin",theta,theta);
  //sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5826 + 8.801 >> hMin",theta,theta);
  
  // Reprocessed data (02/07/12)
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.6152 + 1.196>> hMin",theta,theta);
  sprintf(cmd,"(csc*TMath::Sin(%.4f) + epcrec*TMath::Cos(%.4f))*0.5936 - 0.1598>> hMin",theta,theta);

  t->Draw(cmd,"nsite == 1");

  double cscOff = 0;*/

/*  double fitRange1 = (7000-cscOff)*TMath::Sin(theta) + 2300*TMath::Cos(theta);
  double fitRange2 = (10500-cscOff)*TMath::Sin(theta) + 2950*TMath::Cos(theta);
  double Peak = (8600-cscOff)*TMath::Sin(theta) + 2650*TMath::Cos(theta);
*/
  //double fitRange1 = 2400;
  //double fitRange2 = 2800;
  //double Peak = 2615;

/*  double fitRange1 = 3800;
  double fitRange2 = 4500;
  double Peak = 4215;
*/
  /*TF1 *fitMin = new TF1("fitMin",fitFunction,fitRange1,fitRange2,4);

  fitMin->SetParNames("A1","E1","#sigma","A2");

  fitMin->SetParameters(10000,Peak,100,0.2);

  fitMin->SetParLimits(1,Peak-Peak/12.0,Peak+Peak/12.0);
  fitMin->SetParLimits(2,20,200);
  fitMin->SetParLimits(3,0.05,2.0);

  hMin->Fit("fitMin","r");

  double E1 = fitMin->GetParameter(1);
  double sigma1 = fitMin->GetParameter(2);
  double E1_err = fitMin->GetParError(1);
  double sigma1_err = fitMin->GetParError(2);

  cout << "Resolution from best theta: " << sigma1 / E1 << " +- " << TMath::Sqrt(sigma1_err**2 / E1**2 + (sigma1/E1**2)**2 * E1_err**2)  << endl;*/

  TCanvas *c1 = new TCanvas();
  h2DSingle->Draw("colz");

  TCanvas *c2 = new TCanvas();
  h2DMulti->Draw("colz");

  TCanvas *c3 = new TCanvas();
  c3->Divide(5,4);
  TCanvas *c3_2 = new TCanvas();
  c3_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ss.ps[");
  for (int i = 0; i < nrHistos; i++) {
     //c3->cd(i+1);
     //hSingle[i]->Draw("EZP");
     c3_2->cd();
     //hSingle[i]->Draw("EZP");
     frame_ss[i]->Draw();
     c3_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ss.ps");
  }
  c3_2->cd();
  grSingleRes->Draw("AZP");
  c3_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ss.ps");
  c3_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ss.ps]");

  TCanvas *c4 = new TCanvas();
  c4->Divide(5,4);
  TCanvas *c4_2 = new TCanvas();
  c4_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ms.ps[");
  for (int i = 0; i < nrHistos; i++) {
     //c4->cd(i+1);
     //hMulti[i]->Draw("EZP");
     c4_2->cd();
     //hMulti[i]->Draw("EZP");
     frame_ms[i]->Draw();
     c4_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ms.ps");
  }
  c4_2->cd();
  grMultiRes->Draw("AZP");
  c4_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ms.ps");
  c4_2->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/AntiCorralation/scripts/ThetaMinimization_ms.ps]");

  /*TCanvas *c5 = new TCanvas();
  grSingleRes->Draw("AZP");

  TCanvas *c6 = new TCanvas();
  grMultiRes->Draw("AZP");*/
  //hMin->Draw();

  return;
}
