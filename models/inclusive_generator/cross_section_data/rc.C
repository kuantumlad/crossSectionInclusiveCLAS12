{
  ifstream a("data_5_35_0.5_7.5_7.5GeVInclusiveElastic.dat"); //"data_5_35_1.5_10.6_10.6GeVInclusiveElastic.dat");
  float zhopa;
  float w, q2, cs, csRC;
  float c = 0;
  int t1825 = 0;
  int t2075 = 0;

  const int nBins = 1734;
  float wA1825[nBins];
  float csA1825[nBins];
  float csARC1825[nBins];

  float wA2075[nBins];
  float csA2075[nBins];
  float csARC2075[nBins];

  float RC1825[nBins];
  float RC2075[nBins];


  for (int i = 0; i < 303114*8; i++){
    c++;
    a>>zhopa;
    if (c == 5)
      q2 = zhopa;
    if (c == 6)
      w = zhopa;
    if (c == 7)
      cs = zhopa;
    if (c == 8) {
      csRC = zhopa;
      c = 0;
      if (q2 > 2.5011 && q2 < 2.5015){
	if (t1825 < nBins && cs > 0 && csRC > 0){
	  wA1825[t1825] = w;
	  csA1825[t1825] = cs;
	  csARC1825[t1825] = csRC;
	}
	t1825++;
      }
      if (q2 > 3.333 && q2 < 3.3335){
        if (t2075 < nBins && cs > 0 && csRC > 0){
          wA2075[t2075] = w;
          csA2075[t2075] = cs;
          csARC2075[t2075] = csRC;
        }
        t2075++;
      }

    }
  }

  for (int t = 0; t < nBins;t++){
    RC1825[t] = csA1825[t]/csARC1825[t];
    RC2075[t] = csA2075[t]/csARC2075[t];
  }


  TGraph *csG1825 = new TGraph (t1825, wA1825, csA1825);
  TGraph *csG2075 = new TGraph (t2075, wA2075, csA2075);

  TGraph *csGRC1825 = new TGraph (t1825, wA1825, csARC1825);
  TGraph *csGRC2075 = new TGraph (t2075, wA2075, csARC2075);
  csGRC1825->SetMarkerColor(2);
  csGRC2075->SetMarkerColor(2);


  TGraph *rcG1825 = new TGraph (t1825, wA1825, RC1825);
  TGraph *rcG2075 = new TGraph (t2075, wA2075, RC2075);


  TCanvas *aC = new TCanvas("aC", "aC", 10, 10, 1000, 600);
  aC->Divide(2, 1);
  aC->cd(1);
  csG1825->Draw("A*");
  csGRC1825->Draw("A* same");

  aC->cd(2);
  csG2075->Draw("A*");
  csGRC2075->Draw("A* same");

  aC->SaveAs("sameCS.pdf");
  TCanvas *rcC = new TCanvas("rcC", "rcC", 10, 10, 1000, 600);
  rcC->Divide(2, 1);
  rcC->cd(1);
  rcG1825->Draw("A*");
  rcC->cd(2);
  rcG2075->Draw("A*");
  rcC->SaveAs("rcTest.pdf");
}
