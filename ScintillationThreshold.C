ScintillationThreshold()
{
  gSystem->Load("libEXOUtilities");
  gROOT->SetStyle("Plain");
  
  TChain *t = new TChain("tree");
  
  t->Add("/nfs/slac/g/exo_data3/exo_data/data/WIPP/processed/ateam/1605/recon00001605-*.root");
  
  EXOEventData *ED = 0;
  
  t->SetBranchAddress("EventBranch",&ED);
  
  int nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;
  
  TH1F *h1 = new TH1F("h1","h1",200,0,3500);
  TH1F *h2 = new TH1F("h2","h2",200,0,3500);
  
  bool GoodCluster = true;
  for (int i = 0; i < nentries; i++) {
    t->GetEntry(i);

    int nCL = ED->GetNumChargeClusters();
    
    GoodCluster = true;
    for (int clID = 0; clID < nCL; clID++) {
      EXOChargeCluster *cl = ED->GetChargeCluster(clID);
      for (int clIDTMP = 0; clIDTMP < nCL; clIDTMP++) {
        if (clIDTMP == clID) {continue;}
        EXOChargeCluster *clTMP = ED->GetChargeCluster(clIDTMP);
        if (TMath::Abs(cl->fCollectionTime - clTMP->fCollectionTime) < 120000.0) {GoodCluster = false; break;}
      }
      
      if (!GoodCluster) {continue;}
      
      if (cl->fX == -999 || cl->fY == -999) {continue;}
      if (TMath::Sqrt(cl->fX * cl->fX + cl->fY * cl->fY) > 163.0) {continue;}

      double energy = cl->fCorrectedEnergy*TMath::Exp(120.0/239.0);
      
      EXOScintillationCluster *sc = cl->GetScintillationCluster();
      if (sc) {h1->Fill(energy);}
      h2->Fill(energy);
    }
  }
  
  TH1F *h = h1->Clone("h");
  h->Sumw2();
  h->Divide(h2);
  
  TCanvas *c1 = new TCanvas();
  h->Draw("EZP");
}
