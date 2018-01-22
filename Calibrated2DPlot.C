{
  double Charge_p0_ss = 68.263;
  double Charge_p1_ss = 0.9369;
  double Charge_p0_ms = 113.448;
  double Charge_p1_ms = 0.9382;

  double Scint_p0_ss = -143.2;
  double Scint_p1_ss = 0.3194;
  double Scint_p0_ms = -143.2;
  double Scint_p1_ms = 0.3194;
  
  double Rotated_p0_ss = 3.81;
  double Rotated_p1_ss = 0.6114;
  double Rotated_p0_ms = 7.6951;
  double Rotated_p1_ms = 0.59;

  double Rotated_theta_ss = 0.1860;

  double p0_DiagonalCut_ss = 2717.72;
  double p1_DiagonalCut_ss = 3.4844;
  double p0_DiagonalCut_ms = 2894.05;
  double p1_DiagonalCut_ms = 3.7024;

  gStyle->SetPalette(100);
  gStyle->SetOptLogz(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TFile *f = new TFile("../../EnergyCalibration/analysis/alpha/V2/2424_noReclustering_Fiducial.root","READ");
  TTree *t = (TTree*)f->Get("t");
  
  TH2F *h = new TH2F("h","h",400,500,3500,400,500,3500);
  h->GetXaxis()->SetTitle("ionization (keV)");
  h->GetXaxis()->SetTitleSize(0.03);
  h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->SetTitle("scintillation (keV)");
  h->GetYaxis()->SetTitleSize(0.03);
  h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetYaxis()->SetTitleOffset(1.6);
  
  double Calib_p1 = p1_DiagonalCut_ss*Scint_p1_ss/Charge_p1_ss;
  double Calib_p0 = Scint_p1_ss*(p0_DiagonalCut_ss - p1_DiagonalCut_ss*Charge_p0_ss/Charge_p1_ss);
  
  char *cut = Form("(csc*%.4f + %.4f) <= %.4f * (epcrec*%.4f + %.4f) + %.4f && nsite == 1",Scint_p1_ss, Scint_p0_ss,Calib_p1,Charge_p1_ss,Charge_p0_ss,Calib_p0);
  
  char *form = Form("x*%.4f + %.4f",Calib_p1,Calib_p0);
  
  char *draw = Form("(csc*%.4f+%.4f):(epcrec*%.4f+%.4f)>>h",Scint_p1_ss,Scint_p0_ss,Charge_p1_ss,Charge_p0_ss);
  
  cout << cut << endl;
  cout << form << endl;
  
  t->Draw(draw,"nsite == 1","col goff");
  
  TF1 *DiagCut = new TF1("DiagCut",form,0,3500);
  DiagCut->SetLineWidth(2);
  DiagCut->SetLineStyle(2);
  DiagCut->SetLineColor(kRed);
  
  double xMaxCalib = 3500.0;
  double xMinCalib = 1000.0;
  double tanCalib = Charge_p1_ss / Scint_p1_ss * TMath::Tan(Rotated_theta_ss);
  double thetaCalib = TMath::ATan(tanCalib);
  
  double xMaxCalibRotated = xMaxCalib/TMath::Cos(thetaCalib);
  double xMinCalibRotated = xMinCalib/TMath::Cos(thetaCalib);
  double xMaxRotated = (xMaxCalibRotated - Charge_p0_ss) / Charge_p1_ss;
  double RotatedMax = xMaxRotated * Rotated_p1_ss + Rotated_p0_ss;
  double calib = TMath::Cos(thetaCalib) + TMath::Sin(thetaCalib);
  
  TGaxis *axis = new TGaxis(xMinCalib,xMinCalib*tanCalib,xMaxCalib,xMaxCalib*tanCalib,xMinCalib/calib,xMaxCalibRotated/calib,510);
  axis->SetTitle("rotated axis (keV)");
  axis->SetLineColor(kRed);
  axis->SetTitleColor(kRed);
  axis->CenterTitle(true);
  axis->SetLabelSize(0.03);
  axis->SetTitleSize(0.03);
  axis->SetTitleOffset(2.0);
  axis->SetLabelColor(kRed);
  axis->SetLabelOffset(-0.01);
  
  TLine *l = new TLine(2614.5*calib*TMath::Cos(thetaCalib),2614.5*calib*TMath::Sin(thetaCalib),2614.5*calib*TMath::Cos(thetaCalib) - (3200.0 - 2614.5*calib*TMath::Sin(thetaCalib))*tanCalib,3200);
  l->SetLineColor(kRed);
  l->SetLineWidth(2);
  
  TArc *arc = new TArc(2614.5*calib*TMath::Cos(thetaCalib),2614.5*calib*TMath::Sin(thetaCalib),120.0,90.0+180.0/TMath::Pi()*thetaCalib,180.0+180.0/TMath::Pi()*thetaCalib);
  arc->SetLineColor(kRed);
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  h->Draw("col");
  DiagCut->Draw("same");
  axis->Draw("same");
  l->Draw("same");
  arc->Draw("same");
}