void ThetaPlot()
{
  double true_energy[4] = {662.0, 1173.0, 1333.0, 2615.0};
/*  double res_single[4] = {0.1927, 0.1666, 0.1584, 0.1937};
  double res_multi[4] = {0.2387, 0.2012, 0.1778, 0.2264};
  double res_single_err[4] = {0.0213, 0.0229, 0.0133, 0.0155};
  double res_multi_err[4] = {0.025, 0.0664, 0.0305, 0.01823};
*/
/*  double res_single[4] = {0.1192, 0.1549, 0.1504, 0.19};
  double res_multi[4] = {0.153, 0.1816, 0.1672, 0.2143};
  double res_single_err[4] = {0.0359, 0.0241, 0.0137, 0.0161};
  double res_multi_err[4] = {0.0649, 0.053, 0.026, 0.0097};
*/

// coarse scan
/*  double res_single[4] = {0.1434, 0.1657, 0.1568, 0.195};
  double res_multi[4] = {0.1976, 0.1957, 0.1791, 0.2288};
  double res_single_err[4] = {0.026, 0.0178, 0.0118, 0.0129};
  double res_multi_err[4] = {0.05, 0.0611, 0.0328, 0.01};
*/
// fine scan
/*  double res_single[4] = {0.1422, 0.1677, 0.156, 0.193};
  double res_multi[4] = {0.196, 0.1931, 0.1781, 0.226};
  double res_single_err[4] = {0.0143, 0.032, 0.019, 0.0153};
  double res_multi_err[4] = {0.0331, 0.051, 0.026, 0.0136};
*/
// new fit function
/*  double res_single[4] = {0.1355, 0.1652, 0.1553, 0.1929};
  double res_multi[4] = {0.1828, 0.1931, 0.1803, 0.2233};
  double res_single_err[4] = {0.0022, 0.0027, 0.0015, 0.0008};
  double res_multi_err[4] = {0.0042, 0.0035, 0.0016, 0.0007};
*/
//run by run minimization to include systematics
  double res_single[4] = {0.137573, 0.168237111, 0.154925111, 0.195558};
  double res_multi[4] = {0.185147, 0.204344667, 0.175362556, 0.2267355};
  double res_single_err[4] = {0.007100061, 0.009654485, 0.006004525, 0.00347403};
  double res_multi_err[4] = {0.023103319, 0.016186326, 0.008287013, 0.001243259};
  
//run by run minimization to include systematics, ionization and scintillation calibrated
  double res_single[4] = {0.370816333, 0.470329444, 0.432183222, 0.5369125};
  double res_multi[4] = {0.017532569, 0.018490219, 0.014659815, 0.008002388};

  TGraphErrors *gr_res_single = new TGraphErrors(4,true_energy,res_single,0,res_single_err);
  //TGraphErrors *gr_res_multi = new TGraphErrors(4,true_energy,res_multi,0,res_multi_err);

  gr_res_single->SetMarkerStyle(20);
  gr_res_single->SetMarkerSize(0.8);

/*  gr_res_multi->SetMarkerStyle(20);
  gr_res_multi->SetMarkerSize(0.8);
  gr_res_multi->SetMarkerColor(kRed);
*/
  gr_res_single->GetXaxis()->SetTitle("energy [keV]");
  gr_res_single->GetYaxis()->SetTitle("rotation angle [rad]");

  cout << gr_res_single->GetMean(2) << endl;

/*  TLegend *l = new TLegend(0.8,0.8,0.9,0.9);
  l->AddEntry(gr_res_single,"single site","lp");
  l->AddEntry(gr_res_multi,"multi site","lp");
*/
  TCanvas *c1 = new TCanvas();
  gr_res_single->Draw("AP");
  //gr_res_multi->Draw("Psame");
  //l->Draw("same");

  return;
}
