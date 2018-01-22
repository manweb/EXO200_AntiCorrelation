void ThetaPlot()
{
  double true_energy[4] = {662.0, 1173.0, 1333.0, 2615.0};
  double res_single[4] = {0.1927, 0.1666, 0.1584, 0.1937};
  double res_multi[4] = {0.2387, 0.2012, 0.1778, 0.2264};
  double res_single_err[4] = {0.0213, 0.0229, 0.0133, 0.0155};
  double res_multi_err[4] = {0.025, 0.0664, 0.0305, 0.01823};

  TGraphErrors *gr_res_single = new TGraphErrors(4,true_energy,res_single,0,res_single_err);
  TGraphErrors *gr_res_multi = new TGraphErrors(4,true_energy,res_multi,0,res_multi_err);

  gr_res_single->SetMarkerStyle(20);
  gr_res_single->SetMarkerSize(0.8);

  gr_res_multi->SetMarkerStyle(20);
  gr_res_multi->SetMarkerSize(0.8);
  gr_res_multi->SetMarkerColor(kRed);

  gr_res_single->GetXaxis()->SetTitle("energy [keV]");
  gr_res_single->GetYaxis()->SetTitle("rotation angle [rad]");

  TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
  l->AddEntry(gr_res_single,"single single","lp");
  l->AddEntry(gr_res_multi,"multi site","lp");

  TCanvas *c1 = new TCanvas();
  gr_res_single->Draw("AP");
  gr_res_multi->Draw("Psame");
  l->Draw("same");

  return;
}
