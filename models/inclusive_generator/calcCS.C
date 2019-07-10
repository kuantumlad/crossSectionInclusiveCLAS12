void calcCS()
{
  float c1, c2, c3, c4, c5, c6, c7, c8;
  vector<float> cs;
  vector<float> theta;
  vector<float> energy;
  ifstream inFile("cross_section_data/data_5_35_0.5_7.5_7.5GeVInclusiveElastic.dat");//data_2_35_0.5_10.6_10.6GeV.dat");
  float thetaStep = 0.5;
  float energyStep = 0.01;
  float csTot;
  float norm = 0;
  while(inFile.good()){
    inFile>>c1>>c2>>c3>>c4>>c5>>c6>>c7>>c8;
    theta.push_back(c1);
    energy.push_back(c2);
    cs.push_back(c4);
    norm = norm + c4*TMath::Sin(c1*TMath::DegToRad())*thetaStep*TMath::DegToRad()*energyStep*2*TMath::Pi();
  }
  cout <<norm <<endl;
}
