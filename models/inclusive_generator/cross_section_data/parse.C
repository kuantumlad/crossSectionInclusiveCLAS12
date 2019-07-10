{
  ifstream a("fort.26");
  float zhopa;
  TH1F *cs[19];
  TH1F *csRC[19];
  TH1F *RC[19];
  for (int q = 0; q < 19; q++){
    cs[q] = new TH1F(Form("csQ%d", q + 1), Form("csQ%d", q + 1), 203, 0.952500, 4.992500);
    csRC[q] = new TH1F(Form("csRCQ%d", q + 1), Form("csRCQ%d", q + 1), 203, 0.952500, 4.992500);
    RC[q] = new TH1F(Form("RCQ%d", q + 1), Form("RCQ%d", q + 1), 203, 0.952500, 4.992500);
  }
  for (int q = 0; q < 19; q++){
    for (int w = 0; w < 203; w++){
      for (int c = 0; c < 8; c++){
	a>>zhopa;
	if (c == 6 && zhopa > 0){
	  cs[q]->SetBinContent(w+1, zhopa);
	  }
	if (c == 7 && zhopa > 0) {
	  csRC[q]->SetBinContent(w+1, zhopa);
	}
      }
    }
  }

  for (int q = 0; q < 19; q++){
    for (int w = 0; w < 203; w++){
      if (csRC[q]->GetBinContent(w+1) > 0){
	RC[q]->SetBinContent(w + 1, cs[q]->GetBinContent(w+1)/csRC[q]->GetBinContent(w+1));
      }
    }
  }

  TCanvas *aC = new TCanvas("aC", "aC", 10, 10, 1000, 600);
  
  aC->Divide(5, 4);
  for (int q = 0; q < 19; q++){
    aC->cd(q+1);
    cs[q]->Draw();
    csRC[q]->SetLineColor(2);
    csRC[q]->Draw("same");
  }
  aC->SaveAs("sameCS.pdf");
  
  TCanvas *rcC = new TCanvas("rcC", "rcC", 10, 10, 1000, 600);
  rcC->Divide(5, 4);
  for (int q = 0; q < 19; q++){
    rcC->cd(q+1);
    RC[q]->Draw();
  }
  rcC->SaveAs("rcTest.pdf");
}
