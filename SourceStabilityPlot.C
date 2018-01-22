void SourceStabilityPlot()
{
  gROOT->SetStyle("Plain");
  
  TTree* tCs137 = new TTree;
  TTree* tCo60 = new TTree;
  TTree* tTh228 = new TTree;
  
  tCs137->ReadFile("Cs137.dat");
  tCo60->ReadFile("Co60.dat");
  tTh228->ReadFile("Th228.dat");
  
  int runID_Cs137;
  int time_Cs137;
  double E_I_Cs137;
  double E_S_Cs137;
  double sigma_I_Cs137;
  double sigma_S_Cs137;
  double E_I_err_Cs137;
  double E_S_err_Cs137;
  double sigma_I_err_Cs137;
  double sigma_S_err_Cs137;
  double A_Cs137;
  double R_Cs137;
  double A_err_Cs137;
  double R_err_Cs137;
  
  int runID_Co60;
  int time_Co60;
  double E_I1_Co60;
  double E_S1_Co60;
  double E_I2_Co60;
  double E_S2_Co60;
  double sigma_I1_Co60;
  double sigma_S1_Co60;
  double sigma_I2_Co60;
  double sigma_S2_Co60;
  double E_I1_err_Co60;
  double E_S1_err_Co60;
  double E_I2_err_Co60;
  double E_S2_err_Co60;
  double sigma_I1_err_Co60;
  double sigma_S1_err_Co60;
  double sigma_I2_err_Co60;
  double sigma_S2_err_Co60;
  double A1_Co60;
  double R1_Co60;
  double A2_Co60;
  double R2_Co60;
  double A1_err_Co60;
  double R1_err_Co60;
  double A2_err_Co60;
  double R2_err_Co60;
  
  int runID_Th228;
  int time_Th228;
  double E_I_Th228;
  double E_S_Th228;
  double sigma_I_Th228;
  double sigma_S_Th228;
  double E_I_err_Th228;
  double E_S_err_Th228;
  double sigma_I_err_Th228;
  double sigma_S_err_Th228;
  double A_Th228;
  double R_Th228;
  double A_err_Th228;
  double R_err_Th228;
  
  int maxRunID = 0;
  
  tCs137->SetBranchAddress("RunNumber",&runID_Cs137);
  tCs137->SetBranchAddress("Time",&time_Cs137);
  tCs137->SetBranchAddress("E_I",&E_I_Cs137);
  tCs137->SetBranchAddress("E_S",&E_S_Cs137);
  tCs137->SetBranchAddress("sigma_I",&sigma_I_Cs137);
  tCs137->SetBranchAddress("sigma_S",&sigma_S_Cs137);
  tCs137->SetBranchAddress("E_I_err",&E_I_err_Cs137);
  tCs137->SetBranchAddress("E_S_err",&E_S_err_Cs137);
  tCs137->SetBranchAddress("sigma_I_err",&sigma_I_err_Cs137);
  tCs137->SetBranchAddress("sigma_S_err",&sigma_S_err_Cs137);
  tCs137->SetBranchAddress("A",&A_Cs137);
  tCs137->SetBranchAddress("R",&R_Cs137);
  tCs137->SetBranchAddress("A_err",&A_err_Cs137);
  tCs137->SetBranchAddress("R_err",&R_err_Cs137);
  
  tCo60->SetBranchAddress("RunNumber",&runID_Co60);
  tCo60->SetBranchAddress("Time",&time_Co60);
  tCo60->SetBranchAddress("E_I1",&E_I1_Co60);
  tCo60->SetBranchAddress("E_S1",&E_S1_Co60);
  tCo60->SetBranchAddress("E_I2",&E_I2_Co60);
  tCo60->SetBranchAddress("E_S2",&E_S2_Co60);
  tCo60->SetBranchAddress("sigma_I1",&sigma_I1_Co60);
  tCo60->SetBranchAddress("sigma_S1",&sigma_S1_Co60);
  tCo60->SetBranchAddress("sigma_I2",&sigma_I2_Co60);
  tCo60->SetBranchAddress("sigma_S2",&sigma_S2_Co60);
  tCo60->SetBranchAddress("E_I1_err",&E_I1_err_Co60);
  tCo60->SetBranchAddress("E_S1_err",&E_S1_err_Co60);
  tCo60->SetBranchAddress("E_I2_err",&E_I2_err_Co60);
  tCo60->SetBranchAddress("E_S2_err",&E_S2_err_Co60);
  tCo60->SetBranchAddress("sigma_I1_err",&sigma_I1_err_Co60);
  tCo60->SetBranchAddress("sigma_S1_err",&sigma_S1_err_Co60);
  tCo60->SetBranchAddress("sigma_I2_err",&sigma_I2_err_Co60);
  tCo60->SetBranchAddress("sigma_S2_err",&sigma_S2_err_Co60);
  tCo60->SetBranchAddress("A1",&A1_Co60);
  tCo60->SetBranchAddress("R1",&R1_Co60);
  tCo60->SetBranchAddress("A2",&A2_Co60);
  tCo60->SetBranchAddress("R2",&R2_Co60);
  tCo60->SetBranchAddress("A1_err",&A1_err_Co60);
  tCo60->SetBranchAddress("R1_err",&R1_err_Co60);
  tCo60->SetBranchAddress("A2_err",&A2_err_Co60);
  tCo60->SetBranchAddress("R2_err",&R2_err_Co60);
  
  tTh228->SetBranchAddress("RunNumber",&runID_Th228);
  tTh228->SetBranchAddress("Time",&time_Th228);
  tTh228->SetBranchAddress("E_I",&E_I_Th228);
  tTh228->SetBranchAddress("E_S",&E_S_Th228);
  tTh228->SetBranchAddress("sigma_I",&sigma_I_Th228);
  tTh228->SetBranchAddress("sigma_S",&sigma_S_Th228);
  tTh228->SetBranchAddress("E_I_err",&E_I_err_Th228);
  tTh228->SetBranchAddress("E_S_err",&E_S_err_Th228);
  tTh228->SetBranchAddress("sigma_I_err",&sigma_I_err_Th228);
  tTh228->SetBranchAddress("sigma_S_err",&sigma_S_err_Th228);
  tTh228->SetBranchAddress("A",&A_Th228);
  tTh228->SetBranchAddress("R",&R_Th228);
  tTh228->SetBranchAddress("A_err",&A_err_Th228);
  tTh228->SetBranchAddress("R_err",&R_err_Th228);
  
  int nCs137 = tCs137->GetEntries();
  int nCo60 = tCo60->GetEntries();
  int nTh228 = tTh228->GetEntries();
  
  double *xCs137 = new double[nCs137];
  double *ICs137 = new double[nCs137];
  double *SCs137 = new double[nCs137];
  double *ISCs137 = new double[nCs137];
  double *SSCs137 = new double[nCs137];
  double *IECs137 = new double[nCs137];
  double *SECs137 = new double[nCs137];
  double *ISECs137 = new double[nCs137];
  double *SSECs137 = new double[nCs137];
  double *ACs137 = new double[nCs137];
  double *AECs137 = new double[nCs137];
  double *RCs137 = new double[nCs137];
  double *RECs137 = new double[nCs137];
  double *xCo60_1 = new double[nCo60];
  double *ICo60_1 = new double[nCo60];
  double *SCo60_1 = new double[nCo60];
  double *ISCo60_1 = new double[nCo60];
  double *SSCo60_1 = new double[nCo60];
  double *IECo60_1 = new double[nCo60];
  double *SECo60_1 = new double[nCo60];
  double *ISECo60_1 = new double[nCo60];
  double *SSECo60_1 = new double[nCo60];
  double *xCo60_2 = new double[nCo60];
  double *ICo60_2 = new double[nCo60];
  double *SCo60_2 = new double[nCo60];
  double *ISCo60_2 = new double[nCo60];
  double *SSCo60_2 = new double[nCo60];
  double *IECo60_2 = new double[nCo60];
  double *SECo60_2 = new double[nCo60];
  double *ISECo60_2 = new double[nCo60];
  double *SSECo60_2 = new double[nCo60];
  double *ACo60_1 = new double[nCo60];
  double *AECo60_1 = new double[nCo60];
  double *RCo60_1 = new double[nCo60];
  double *RECo60_1 = new double[nCo60];
  double *ACo60_2 = new double[nCo60];
  double *AECo60_2 = new double[nCo60];
  double *RCo60_2 = new double[nCo60];
  double *RECo60_2 = new double[nCo60];
  double *xTh228 = new double[nTh228];
  double *ITh228 = new double[nTh228];
  double *STh228 = new double[nTh228];
  double *ISTh228 = new double[nTh228];
  double *SSTh228 = new double[nTh228];
  double *IETh228 = new double[nTh228];
  double *SETh228 = new double[nTh228];
  double *ISETh228 = new double[nTh228];
  double *SSETh228 = new double[nTh228];
  double *ATh228 = new double[nTh228];
  double *AETh228 = new double[nTh228];
  double *RTh228 = new double[nTh228];
  double *RETh228 = new double[nTh228];
  
  for (int i = 0; i < nCs137; i++) {
    tCs137->GetEntry(i);
    
    xCs137[i] = runID_Cs137;
    ICs137[i] = E_I_Cs137;
    SCs137[i] = E_S_Cs137;
    ISCs137[i] = sigma_I_Cs137;
    SSCs137[i] = sigma_S_Cs137;
    IECs137[i] = E_I_err_Cs137;
    SECs137[i] = E_S_err_Cs137;
    ISECs137[i] = sigma_I_err_Cs137;
    SSECs137[i] = sigma_S_err_Cs137;
    ACs137[i] = A_Cs137;
    AECs137[i] = A_err_Cs137;
    RCs137[i] = R_Cs137;
    RECs137[i] = R_err_Cs137;
    
    if (runID_Cs137 > maxRunID) {maxRunID = runID_Cs137;}
  }
  
  for (int i = 0; i < nCo60; i++) {
    tCo60->GetEntry(i);
    
    xCo60_1[i] = runID_Co60;
    ICo60_1[i] = E_I1_Co60;
    SCo60_1[i] = E_S1_Co60;
    ISCo60_1[i] = sigma_I1_Co60;
    SSCo60_1[i] = sigma_S1_Co60;
    IECo60_1[i] = E_I1_err_Co60;
    SECo60_1[i] = E_S1_err_Co60;
    ISECo60_1[i] = sigma_I1_err_Co60;
    SSECo60_1[i] = sigma_S1_err_Co60;
    ICo60_2[i] = E_I2_Co60;
    SCo60_2[i] = E_S2_Co60;
    ISCo60_2[i] = sigma_I2_Co60;
    SSCo60_2[i] = sigma_S2_Co60;
    IECo60_2[i] = E_I2_err_Co60;
    SECo60_2[i] = E_S2_err_Co60;
    ISECo60_2[i] = sigma_I2_err_Co60;
    SSECo60_2[i] = sigma_S2_err_Co60;
    ACo60_1[i] = A1_Co60;
    AECo60_1[i] = A1_err_Co60;
    RCo60_1[i] = R1_Co60;
    RECo60_1[i] = R1_err_Co60;
    ACo60_2[i] = A2_Co60;
    AECo60_2[i] = A2_err_Co60;
    RCo60_2[i] = R2_Co60;
    RECo60_2[i] = R2_err_Co60;
    
    if (runID_Co60 > maxRunID) {maxRunID = runID_Co60;}
  }
  
  for (int i = 0; i < nTh228; i++) {
    tTh228->GetEntry(i);
    
    xTh228[i] = runID_Th228;
    ITh228[i] = E_I_Th228;
    STh228[i] = E_S_Th228;
    ISTh228[i] = sigma_I_Th228;
    SSTh228[i] = sigma_S_Th228;
    IETh228[i] = E_I_err_Th228;
    SETh228[i] = E_S_err_Th228;
    ISETh228[i] = sigma_I_err_Th228;
    SSETh228[i] = sigma_S_err_Th228;
    ATh228[i] = A_Th228;
    AETh228[i] = A_err_Th228;
    RTh228[i] = R_Th228;
    RETh228[i] = R_err_Th228;
    
    if (runID_Th228 > maxRunID) {maxRunID = runID_Th228;}
  }
  
  TGraphErrors *grCs137_I = new TGraphErrors(nCs137,xCs137,ICs137,0,IECs137);
  TGraphErrors *grCs137_S = new TGraphErrors(nCs137,xCs137,SCs137,0,SECs137);
  TGraphErrors *grCs137_IS = new TGraphErrors(nCs137,xCs137,ISCs137,0,ISECs137);
  TGraphErrors *grCs137_SS = new TGraphErrors(nCs137,xCs137,SSCs137,0,SSECs137);
  TGraphErrors *grCs137_A = new TGraphErrors(nCs137,xCs137,ACs137,0,AECs137);
  TGraphErrors *grCs137_R = new TGraphErrors(nCs137,xCs137,RCs137,0,RECs137);
  TGraphErrors *grCo60_1_I = new TGraphErrors(nCo60,xCo60_1,ICo60_1,0,IECo60_1);
  TGraphErrors *grCo60_1_S = new TGraphErrors(nCo60,xCo60_1,SCo60_1,0,SECo60_1);
  TGraphErrors *grCo60_1_IS = new TGraphErrors(nCo60,xCo60_1,ISCo60_1,0,ISECo60_1);
  TGraphErrors *grCo60_1_SS = new TGraphErrors(nCo60,xCo60_1,SSCo60_1,0,SSECo60_1);
  TGraphErrors *grCo60_2_I = new TGraphErrors(nCo60,xCo60_1,ICo60_2,0,IECo60_2);
  TGraphErrors *grCo60_2_S = new TGraphErrors(nCo60,xCo60_1,SCo60_2,0,SECo60_2);
  TGraphErrors *grCo60_2_IS = new TGraphErrors(nCo60,xCo60_1,ISCo60_2,0,ISECo60_2);
  TGraphErrors *grCo60_2_SS = new TGraphErrors(nCo60,xCo60_1,SSCo60_2,0,SSECo60_2);
  TGraphErrors *grCo60_1_A = new TGraphErrors(nCo60,xCo60_1,ACo60_1,0,AECo60_1);
  TGraphErrors *grCo60_1_R = new TGraphErrors(nCo60,xCo60_1,RCo60_1,0,RECo60_1);
  TGraphErrors *grCo60_2_A = new TGraphErrors(nCo60,xCo60_1,ACo60_2,0,AECo60_2);
  TGraphErrors *grCo60_2_R = new TGraphErrors(nCo60,xCo60_1,RCo60_2,0,RECo60_2);
  TGraphErrors *grTh228_I = new TGraphErrors(nTh228,xTh228,ITh228,0,IETh228);
  TGraphErrors *grTh228_S = new TGraphErrors(nTh228,xTh228,STh228,0,SETh228);
  TGraphErrors *grTh228_IS = new TGraphErrors(nTh228,xTh228,ISTh228,0,ISETh228);
  TGraphErrors *grTh228_SS = new TGraphErrors(nTh228,xTh228,SSTh228,0,SSETh228);
  TGraphErrors *grTh228_A = new TGraphErrors(nTh228,xTh228,ATh228,0,AETh228);
  TGraphErrors *grTh228_R = new TGraphErrors(nTh228,xTh228,RTh228,0,RETh228);
  
  grCs137_I->SetTitle("Ionization");
  grCs137_S->SetTitle("Scintillation");
  grCs137_IS->SetTitle("sigma ionization");
  grCs137_SS->SetTitle("sigma scintillation");
  
  grCs137_I->GetXaxis()->SetTitle("run number");
  grCs137_S->GetXaxis()->SetTitle("run number");
  grCs137_IS->GetXaxis()->SetTitle("run number");
  grCs137_SS->GetXaxis()->SetTitle("run number");
  grCs137_A->GetXaxis()->SetTitle("run number");
  grCs137_R->GetXaxis()->SetTitle("run number");
  
  grCs137_I->GetYaxis()->SetTitle("peak position");
  grCs137_S->GetYaxis()->SetTitle("peak position");
  grCs137_IS->GetYaxis()->SetTitle("peak sigma");
  grCs137_SS->GetYaxis()->SetTitle("peak sigma");
  grCs137_A->GetYaxis()->SetTitle("peak position");
  grCs137_R->GetYaxis()->SetTitle("resolution");
  
  grCs137_I->GetYaxis()->SetTitleOffset(1.2);
  grCs137_S->GetYaxis()->SetTitleOffset(1.2);
  grCs137_IS->GetYaxis()->SetTitleOffset(1.2);
  grCs137_SS->GetYaxis()->SetTitleOffset(1.2);
  grCs137_A->GetYaxis()->SetTitleOffset(1.2);
  grCs137_R->GetYaxis()->SetTitleOffset(1.2);
  
  grCo60_1_I->GetXaxis()->SetTitle("run number");
  grCo60_1_S->GetXaxis()->SetTitle("run number");
  grCo60_1_IS->GetXaxis()->SetTitle("run number");
  grCo60_1_SS->GetXaxis()->SetTitle("run number");
  grCo60_1_A->GetXaxis()->SetTitle("run number");
  grCo60_1_R->GetXaxis()->SetTitle("run number");
  
  grCo60_1_I->GetYaxis()->SetTitle("peak position");
  grCo60_1_S->GetYaxis()->SetTitle("peak position");
  grCo60_1_IS->GetYaxis()->SetTitle("peak sigma");
  grCo60_1_SS->GetYaxis()->SetTitle("peak sigma");
  grCo60_1_A->GetYaxis()->SetTitle("peak position");
  grCo60_1_R->GetYaxis()->SetTitle("resolution");
  
  grCo60_1_I->GetYaxis()->SetTitleOffset(1.2);
  grCo60_1_S->GetYaxis()->SetTitleOffset(1.2);
  grCo60_1_IS->GetYaxis()->SetTitleOffset(1.2);
  grCo60_1_SS->GetYaxis()->SetTitleOffset(1.2);
  grCo60_1_A->GetYaxis()->SetTitleOffset(1.2);
  grCo60_1_R->GetYaxis()->SetTitleOffset(1.2);
  
  grCo60_2_I->GetXaxis()->SetTitle("run number");
  grCo60_2_S->GetXaxis()->SetTitle("run number");
  grCo60_2_IS->GetXaxis()->SetTitle("run number");
  grCo60_2_SS->GetXaxis()->SetTitle("run number");
  grCo60_2_A->GetXaxis()->SetTitle("run number");
  grCo60_2_R->GetXaxis()->SetTitle("run number");
  
  grCo60_2_I->GetYaxis()->SetTitle("peak position");
  grCo60_2_S->GetYaxis()->SetTitle("peak position");
  grCo60_2_IS->GetYaxis()->SetTitle("peak sigma");
  grCo60_2_SS->GetYaxis()->SetTitle("peak sigma");
  grCo60_2_A->GetYaxis()->SetTitle("peak position");
  grCo60_2_R->GetYaxis()->SetTitle("resolution");
  
  grCo60_2_I->GetYaxis()->SetTitleOffset(1.2);
  grCo60_2_S->GetYaxis()->SetTitleOffset(1.2);
  grCo60_2_IS->GetYaxis()->SetTitleOffset(1.2);
  grCo60_2_SS->GetYaxis()->SetTitleOffset(1.2);
  grCo60_2_A->GetYaxis()->SetTitleOffset(1.2);
  grCo60_2_R->GetYaxis()->SetTitleOffset(1.2);
  
  grTh228_I->GetXaxis()->SetTitle("run number");
  grTh228_S->GetXaxis()->SetTitle("run number");
  grTh228_IS->GetXaxis()->SetTitle("run number");
  grTh228_SS->GetXaxis()->SetTitle("run number");
  grTh228_A->GetXaxis()->SetTitle("run number");
  grTh228_R->GetXaxis()->SetTitle("run number");
  
  grTh228_I->GetYaxis()->SetTitle("peak position");
  grTh228_S->GetYaxis()->SetTitle("peak position");
  grTh228_IS->GetYaxis()->SetTitle("peak sigma");
  grTh228_SS->GetYaxis()->SetTitle("peak sigma");
  grTh228_A->GetYaxis()->SetTitle("peak position");
  grTh228_R->GetYaxis()->SetTitle("resolution");
  
  grTh228_I->GetYaxis()->SetTitleOffset(1.2);
  grTh228_S->GetYaxis()->SetTitleOffset(1.2);
  grTh228_IS->GetYaxis()->SetTitleOffset(1.2);
  grTh228_SS->GetYaxis()->SetTitleOffset(1.2);
  grTh228_A->GetYaxis()->SetTitleOffset(1.2);
  grTh228_A->GetYaxis()->SetTitleOffset(1.2);
  
  grCs137_I->GetYaxis()->SetRangeUser(500,3000);
  grCs137_S->GetYaxis()->SetRangeUser(1000,10000);
  grCs137_IS->GetYaxis()->SetRangeUser(40,150);
  grCs137_SS->GetYaxis()->SetRangeUser(300,700);
  
  grCs137_I->GetXaxis()->SetLimits(2400,maxRunID+100);
  grCs137_S->GetXaxis()->SetLimits(2400,maxRunID+100);
  grCs137_IS->GetXaxis()->SetLimits(2400,maxRunID+100);
  grCs137_SS->GetXaxis()->SetLimits(2400,maxRunID+100);
  
  grCs137_I->SetMarkerStyle(20);
  grCs137_S->SetMarkerStyle(20);
  grCs137_IS->SetMarkerStyle(20);
  grCs137_SS->SetMarkerStyle(20);
  grCs137_A->SetMarkerStyle(20);
  grCs137_R->SetMarkerStyle(20);
  grCo60_1_I->SetMarkerStyle(22);
  grCo60_1_S->SetMarkerStyle(22);
  grCo60_1_IS->SetMarkerStyle(22);
  grCo60_1_SS->SetMarkerStyle(22);
  grCo60_1_A->SetMarkerStyle(22);
  grCo60_1_R->SetMarkerStyle(22);
  grCo60_2_I->SetMarkerStyle(23);
  grCo60_2_S->SetMarkerStyle(23);
  grCo60_2_IS->SetMarkerStyle(23);
  grCo60_2_SS->SetMarkerStyle(23);
  grCo60_2_A->SetMarkerStyle(23);
  grCo60_2_R->SetMarkerStyle(23);
  grTh228_I->SetMarkerStyle(21);
  grTh228_S->SetMarkerStyle(21);
  grTh228_IS->SetMarkerStyle(21);
  grTh228_SS->SetMarkerStyle(21);
  grTh228_A->SetMarkerStyle(21);
  grTh228_R->SetMarkerStyle(21);
  
  grCs137_I->SetMarkerSize(0.8);
  grCs137_S->SetMarkerSize(0.8);
  grCs137_IS->SetMarkerSize(0.8);
  grCs137_SS->SetMarkerSize(0.8);
  grCs137_A->SetMarkerSize(0.8);
  grCs137_R->SetMarkerSize(0.8);
  grCo60_1_I->SetMarkerSize(0.8);
  grCo60_1_S->SetMarkerSize(0.8);
  grCo60_1_IS->SetMarkerSize(0.8);
  grCo60_1_SS->SetMarkerSize(0.8);
  grCo60_1_A->SetMarkerSize(0.8);
  grCo60_1_R->SetMarkerSize(0.8);
  grCo60_2_I->SetMarkerSize(0.8);
  grCo60_2_S->SetMarkerSize(0.8);
  grCo60_2_IS->SetMarkerSize(0.8);
  grCo60_2_SS->SetMarkerSize(0.8);
  grCo60_2_A->SetMarkerSize(0.8);
  grCo60_2_R->SetMarkerSize(0.8);
  grTh228_I->SetMarkerSize(0.8);
  grTh228_S->SetMarkerSize(0.8);
  grTh228_IS->SetMarkerSize(0.8);
  grTh228_SS->SetMarkerSize(0.8);
  grTh228_A->SetMarkerSize(0.8);
  grTh228_R->SetMarkerSize(0.8);
  
  char Cs137_I_label[50];
  char Cs137_S_label[50];
  char Co60_I1_label[50];
  char Co60_S1_label[50];
  char Co60_I2_label[50];
  char Co60_S2_label[50];
  char Th228_I_label[50];
  char Th228_S_label[50];
  sprintf(Cs137_I_label,"Cs137 (%.2f +- %.2f)",TMath::Mean(nCs137,ICs137),TMath::RMS(nCs137,ICs137));
  sprintf(Cs137_S_label,"Cs137 (%.2f +- %.2f)",TMath::Mean(nCs137,SCs137),TMath::RMS(nCs137,SCs137));
  sprintf(Co60_I1_label,"Co60 (1) (%.2f +- %.2f)",TMath::Mean(nCo60,ICo60_1),TMath::RMS(nCo60,ICo60_1));
  sprintf(Co60_S1_label,"Co60 (1) (%.2f +- %.2f)",TMath::Mean(nCo60,SCo60_1),TMath::RMS(nCo60,SCo60_1));
  sprintf(Co60_I2_label,"Co60 (2) (%.2f +- %.2f)",TMath::Mean(nCo60,ICo60_2),TMath::RMS(nCo60,ICo60_2));
  sprintf(Co60_S2_label,"Co60 (2) (%.2f +- %.2f)",TMath::Mean(nCo60,SCo60_2),TMath::RMS(nCo60,SCo60_2));
  sprintf(Th228_I_label,"Th228 (%.2f +- %.2f)",TMath::Mean(nTh228,ITh228),TMath::RMS(nTh228,ITh228));
  sprintf(Th228_S_label,"Th228 (%.2f +- %.2f)",TMath::Mean(nTh228,STh228),TMath::RMS(nTh228,STh228));
  
  TLegend *l1 = new TLegend(0.65,0.15,0.88,0.4);
  l1->AddEntry(grCs137_I,Cs137_I_label,"lp");
  l1->AddEntry(grCo60_1_I,Co60_I1_label,"lp");
  l1->AddEntry(grCo60_2_I,Co60_I2_label,"lp");
  l1->AddEntry(grTh228_I,Th228_I_label,"lp");
  
  TLegend *l2 = new TLegend(0.65,0.15,0.88,0.4);
  l2->AddEntry(grCs137_S,Cs137_S_label,"lp");
  l2->AddEntry(grCo60_1_S,Co60_S1_label,"lp");
  l2->AddEntry(grCo60_2_S,Co60_S2_label,"lp");
  l2->AddEntry(grTh228_S,Th228_S_label,"lp");
  
  l1->SetFillColor(0);
  l2->SetFillColor(0);
  
  TCanvas *c1 = new TCanvas();
  grCs137_I->Draw("AZP");
  grCo60_1_I->Draw("ZPsame");
  grCo60_2_I->Draw("ZPsame");
  grTh228_I->Draw("ZPsame");
  l1->Draw("same");
  c1->SaveAs("Ionization.png","RECREATE");
  
  TCanvas *c2 = new TCanvas();
  grCs137_S->Draw("AZP");
  grCo60_1_S->Draw("ZPsame");
  grCo60_2_S->Draw("ZPsame");
  grTh228_S->Draw("ZPsame");
  l2->Draw("same");
  c2->SaveAs("Scintillation.png","RECREATE");
  
  TCanvas *c3 = new TCanvas();
  grCs137_IS->Draw("AZP");
  grCo60_1_IS->Draw("ZPsame");
  grCo60_2_IS->Draw("ZPsame");
  grTh228_IS->Draw("ZPsame");
  c3->SaveAs("SigmaIonization.png","RECREATE");
  
  TCanvas *c4 = new TCanvas();
  grCs137_SS->Draw("AZP");
  grCo60_1_SS->Draw("ZPsame");
  grCo60_2_SS->Draw("ZPsame");
  grTh228_SS->Draw("ZPsame");
  c4->SaveAs("SigmaScintillation.png","RECREATE");
  
  // individual plots
  grCs137_I->GetYaxis()->SetRangeUser(TMath::MinElement(nCs137,ICs137) - 20, TMath::MaxElement(nCs137,ICs137) + 20);
  grCs137_S->GetYaxis()->SetRangeUser(TMath::MinElement(nCs137,SCs137) - 50, TMath::MaxElement(nCs137,SCs137) + 50);
  grCs137_I->GetXaxis()->SetRangeUser(TMath::MinElement(nCs137,xCs137) - 50, TMath::MaxElement(nCs137,xCs137) + 50);
  grCs137_S->GetXaxis()->SetRangeUser(TMath::MinElement(nCs137,xCs137) - 50, TMath::MaxElement(nCs137,xCs137) + 50);
  
  grCs137_I->SetTitle("Cs137");
  grCs137_S->SetTitle("Cs137");
  grCs137_A->SetTitle("Cs137");
  grCs137_R->SetTitle("Cs137");
  grCo60_1_I->SetTitle("Co60 (1)");
  grCo60_1_S->SetTitle("Co60 (1)");
  grCo60_1_A->SetTitle("Co60 (1)");
  grCo60_1_R->SetTitle("Co60 (1)");
  grCo60_2_I->SetTitle("Co60 (2)");
  grCo60_2_S->SetTitle("Co60 (2)");
  grCo60_2_A->SetTitle("Co60 (2)");
  grCo60_2_R->SetTitle("Co60 (2)");
  grTh228_I->SetTitle("Th228");
  grTh228_S->SetTitle("Th228");
  grTh228_A->SetTitle("Th228");
  grTh228_R->SetTitle("Th228");
  
  TCanvas *c5 = new TCanvas();
  c5->Divide(2,2);
  c5->cd(1);
  grCs137_I->Draw("AZP");
  c5->cd(2);
  grCo60_1_I->Draw("AZP");
  c5->cd(3);
  grCo60_2_I->Draw("AZP");
  c5->cd(4);
  grTh228_I->Draw("AZP");
  c5->SaveAs("IonizationIndividual.png");
  
  TCanvas *c6 = new TCanvas();
  c6->Divide(2,2);
  c6->cd(1);
  grCs137_S->Draw("AZP");
  c6->cd(2);
  grCo60_1_S->Draw("AZP");
  c6->cd(3);
  grCo60_2_S->Draw("AZP");
  c6->cd(4);
  grTh228_S->Draw("AZP");
  c6->SaveAs("ScintillationIndividual.png");
  
  TCanvas *c7 = new TCanvas();
  c7->Divide(2,2);
  c7->cd(1);
  grCs137_A->Draw("AZP");
  c7->cd(2);
  grCo60_1_A->Draw("AZP");
  c7->cd(3);
  grCo60_2_A->Draw("AZP");
  c7->cd(4);
  grTh228_A->Draw("AZP");
  c7->SaveAs("AnticorrelationIndividual.png");
  
  TCanvas *c8 = new TCanvas();
  c8->Divide(2,2);
  c8->cd(1);
  grCs137_R->Draw("AZP");
  c8->cd(2);
  grCo60_1_R->Draw("AZP");
  c8->cd(3);
  grCo60_2_R->Draw("AZP");
  c8->cd(4);
  grTh228_R->Draw("AZP");
  c8->SaveAs("ResolutionIndividual.png");
  
  return;
}