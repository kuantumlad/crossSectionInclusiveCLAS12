/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
{

  gSystem->CompileMacro("DisFunctions.C", "kf");
  gSystem->CompileMacro("GenFunctions.C", "kf");
  gSystem->CompileMacro("InclusiveRateAnalysis.C");
  t=new InclusiveRateAnalysis();

  float beamE = 7.5;

  //2 GeV 
  if (beamE > 2.1 && beamE < 2.3){
    t->SetKinematicLimits(5, 35, 0.005, 0.1, 2.222, 0.002, 5472378, 2.222);
    t->SetFile("cross_section_data/data_5_35_0.1_2.2_2.2GeVInclusiveElastic.dat");
  }
  if (beamE > 10.5 && beamE < 10.7){
    //10 GeV
    t->SetKinematicLimits( 5, 35, 0.1, 0.5, 10.6, 0.005, 333214, 10.604);
    t->SetFile("cross_section_data/data_5_35_1.0_10.6_10.6GeVInclusiveElastic.dat");
  }
  //7.5 GeV
  if (beamE > 7.4 && beamE < 7.6){
    t->SetKinematicLimits( 5, 35, 0.1, 0.5, 7.546, 0.005,  279779, 7.546);
    t->SetFile("cross_section_data/data_5_35_0.5_7.5_7.5GeVInclusiveElastic.dat");
  }

  //6.5 GeV
  if (beamE > 6.4 && beamE < 6.6){
    t->SetKinematicLimits( 5, 35, 0.1, 0.5, 6.535, 0.005, 248648, 6.535);
    t->SetFile("cross_section_data/data_5_35_0.5_6.5_6.5GeVInclusiveElastic.dat");
  }


  cout << beamE <<endl;
  //cout <<  Form("/work/clas12/bclary/CLAS12/inclusive_studies/lunds/incl_%3.1/rad", beamE) << endl;
  
  t->PrepareOutput("");
  t->BookHistos();
  t->SetLuminosity( pow( 10, 35 ) );
  t->ReadDataFile();
  t->PlotCrossSectionDependences();
  t->GeneratePseudoData(10000000, 10000, "/work/clas12/bclary/CLAS12/inclusive_studies/lunds/incl_7p5/norad/");
  t->CompareGenEventsToOriginal();
  t->SaveOutput();
  delete t;
  
}
