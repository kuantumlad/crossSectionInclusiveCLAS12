{
  float beamE = 6.5;
  
  float wLow = 0.8;
  float wHigh = 2.3;

  float q2Low = 0.05;
  float q2High = 2;
  

  if (beamE > 7.4 && beamE < 7.6){
    TFile *f = new TFile("data_5_35_0.5_7.5_7.5GeVInclusiveElastic.root");
    wHigh = 4;
    q2High = 6;

  }

  if (beamE > 6.4 && beamE < 6.6){
    TFile *f = new TFile("data_5_35_0.5_6.5_6.5GeVInclusiveElastic.root");
    wHigh = 4;
    q2High = 6;

  }

  if (beamE > 10.5 && beamE < 10.7){
    TFile *f = new TFile("data_5_35_1.0_10.6_10.6GeVInclusiveElastic.root");
    wHigh = 4.5;
    q2High = 10;
  }

  if (beamE > 2.1 && beamE < 2.3){
    TFile *f = new TFile("data_5_35_0.1_2.2_2.2GeVInclusiveElastic.root");
  }
  hQ2vsWGen->SetAxisRange(wLow, wHigh, "x");
  hQ2vsWGen->SetAxisRange(q2Low, q2High, "y");
  hWGen->SetAxisRange(wLow, wHigh, "x");
  
  TCanvas *lala = new TCanvas("lala", "lala", 10, 10, 800, 800);
  lala->Divide(2, 2);
  lala->cd(1);
  gPad->SetLogz();
  hQ2vsWGen->Draw("COLZ");
  lala->cd(2);
  hWGen->Draw("COLZ");
  lala->cd(3);
  hThetaGen->Draw("");
  lala->cd(4);
  hEnergyvsThetaGen->Draw("COLZ");
  gPad->SetLogz();
  lala->SaveAs(Form("monitoringE%3.1f.pdf", beamE));
}
