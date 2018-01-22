void Superimposed()
{
  gStyle->SetPalette(1);

  TFile *f1 = new TFile("../../EnergyCalibration/analysis/V6/2448_noReclustering_Fiducial.root","READ");
  TFile *f2 = new TFile("../../EnergyCalibration/analysis/V6/2526_noReclustering_Fiducial.root","READ");
  TFile *f3 = new TFile("../../EnergyCalibration/analysis/V6/2449_noReclustering_Fiducial.root","READ");
  TFile *f4 = new TFile("2417_2448.root","READ");
  TFile *f5 = new TFile("../../LowBackgroundData/analysis/LowBackgroundDataMasked.root","READ");

  TTree *t1 = (TTree*)f1->Get("t");
  TTree *t2 = (TTree*)f2->Get("t");
  TTree *t3 = (TTree*)f3->Get("t");
  TTree *t4 = (TTree*)f4->Get("treeUP");
  TTree *t5 = (TTree*)f5->Get("t");

  TH2F *h1 = new TH2F("h1","Scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h2 = new TH2F("h2","Scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h3 = new TH2F("h3","Scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h4 = new TH2F("h4","Scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h5 = new TH2F("h5","Scintillation vs ionization",200,0,3500,200,0,15000);
  TH2F *h6 = new TH2F("h6","Scintillation vs ionization",200,0,3500,200,0,15000);

  t1->Draw("csc:ecrec>>h1","nsite == 1","colz");
  //t4->Draw("(fcsc-2*2087.84):fepclUP>>h4","","colzsame");
  t4->Draw("fcsc:fesum>>h6","","same");
  t2->Draw("csc:ecrec>>h2","nsite == 1","colzsame");
  t3->Draw("csc:ecrec>>h3","nsite == 1","colzsame");
  t5->Draw("csc:ecrec>>h5","nsite == 1","same");

  double x[4] = {585.435, 1121.55, 1301.35, 2632.44};
  double y[4] = {2509.31, 4180.87, 4682.36, 8608.48};

  //double x[4] = {562.51, 1087.96, 1261.99, 2587.56};
  //double y[4] = {2628.6, 4329.76, 4821.81, 8795.56};

  TGraph *gr = new TGraph(4,x,y);

  TF1 *fit = new TF1("fit","[0]+[1]*x",0,3500);
  gr->Fit("fit","r");

  gr->Draw("Psame");

}
