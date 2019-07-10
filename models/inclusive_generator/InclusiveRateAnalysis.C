/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#include<InclusiveRateAnalysis.h>

#ifndef _INCLUSIVERATEANALYSIS_C_
#define _INCLUSIVERATEANALYSIS_C_

ClassImp(InclusiveRateAnalysis)

InclusiveRateAnalysis::InclusiveRateAnalysis()
{

  wCut  = 0.0 ;

  q2Cut = 0.0 ;

  integratedCrossSection    = 0.0 ;
			          
  integratedCrossSectionRad = 0.0 ;

  mapLundIdToNumber.clear();

  MapLundIdToParticleNumber( mapLundIdToNumber );

  mapLundIdToMass.clear();

  MapLundIdToParticleMass( mapLundIdToMass );

  /* here the defaul value of L=10^35 cm-2 s-1 is set. It can be set from script and commented here, if needed */
  // SetLuminosity( pow( 10, 35 ) );

}


void InclusiveRateAnalysis::SetFile( const char* _file_name )
{

  file_name = _file_name ;

  cout << " Analyzing file " << file_name << endl ;

}

/* -------------------------- */
/*  Read data file from Java  */
/* -------------------------- */
void InclusiveRateAnalysis::ReadJavaDataFiles( const char* _java_file_name, Int_t _n_java_files )
{

  /* clear vectors to store values */
  ClearVectors();

  integratedCrossSectionRad = 1.01242; /* this is for 11 GeV: 0.370401 ; */

  Int_t _n_gen_events = 0 ;

  ifstream _input_java_file ;

  for( Int_t _if = 0 ; _if < _n_java_files ; _if++ ){

    cout << " Now reading file " << Form( "%s%d.dat", _java_file_name, _if) << endl ;

    _input_java_file.open( Form( "%s%d.dat", _java_file_name, _if) );

    if (_input_java_file.is_open())
      {
	while ( _input_java_file.good() )
	  {

	    TString _event_type = "";
	    Double_t _px, _py, _pz, _theta, _phi ;

	    _input_java_file >> _event_type >> _px >> _py >> _pz >> _theta >> _phi ;

	    Double_t _energy = TMath::Sqrt( _px*_px + _py*_py + _pz*_pz + ElectronMass*ElectronMass );

	    if(      !_event_type.CompareTo( "gen" )                ) { energyEleGen.push_back( _energy ); thetaEleGen.push_back( _theta );  phiEleGen.push_back( _phi ); hEnergyGen -> Fill( _energy ); hEnergyvsThetaGen -> Fill( _theta, _energy ); hThetavsPhiGen -> Fill( _phi, _theta ); _n_gen_events++ ; }
	    else if( !_event_type.CompareTo( "rec" ) && _pz > -99.9 ) { energyEleRec.push_back( _energy ); thetaEleRec.push_back( _theta );  phiEleRec.push_back( _phi ); hEnergy    -> Fill( _energy ); hEnergyvsTheta    -> Fill( _theta, _energy ); hThetavsPhi    -> Fill( _phi, _theta ); }
	    // else       cout << " WARNING! Something is wrong with the format of the java dat file! " << endl ;

	  }
      }

    _input_java_file.close();

  } /* it closes the loop on the files */

  cout << " Analyzed " << _n_gen_events << " generated events " << endl ;

  /* check the distribution of the cross-section */
  TCanvas* _c_kin = new TCanvas("c_java_kin_check", "", 1000, 700 );

  _c_kin -> Divide( 2, 2 );

  _c_kin -> cd( 1 ) ;

  gPad -> SetLogz( 1 ) ;

  hEnergyvsThetaGen -> Draw( "colz" );
   
  _c_kin -> cd( 2 ) ;

  gPad -> SetLogz( 1 ) ;

  hEnergyvsTheta    -> Draw( "colz" );
   
  _c_kin -> cd( 3 ) ;

  gPad -> SetLogz( 1 ) ;

  hThetavsPhiGen    -> Draw( "colz" );
   
  _c_kin -> cd( 4 ) ;

  gPad -> SetLogz( 1 ) ;

  hThetavsPhi       -> Draw( "colz" );

  _c_kin -> SaveAs( Form( "%s.pdf", _c_kin -> GetName() ) );

  output_file -> WriteTObject( _c_kin );

  // _c_gen -> Divide( 2 );

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  TLegend* _leg = new TLegend( 0.91, 0.16, 0.99, 0.87 );

  _leg -> SetFillStyle( 0 );

  _leg -> SetLineColor( 0 );

  /* check the distribution of the cross-section */
  TCanvas* _c_gen = new TCanvas("c_java_gen_check", "", 1000, 700 );

  _c_gen -> cd();

  /* generated cross-sections */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_temp_py = ( hEnergyvsThetaGen -> ProjectionY( Form("h_java_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    //  if( !_it ){

    //   _h_temp_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

    //   _h_temp_py -> GetXaxis() -> SetLabelFont(   132     );
	 	      
    //   _h_temp_py -> GetYaxis() -> SetLabelFont(   132     );
	 	      
    //   _h_temp_py -> GetZaxis() -> SetLabelFont(   132     );
	 	      
    //   _h_temp_py -> GetXaxis() -> SetTitleFont(   132     );
	 	      
    //   _h_temp_py -> GetYaxis() -> SetTitleFont(   132     );
	 	      
    //   _h_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    // }

    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    // cout << " Color set to " << _color << "\t" << _it << "\t" << _i_multiple << endl ;

    _h_temp_py -> SetLineColor(   _color );

    _h_temp_py -> SetMarkerColor( _color );

    // cout << _h_temp_py -> GetEntries() << " \t " << _it << endl ;

    _leg -> AddEntry( _h_temp_py, Form( "#theta = %g", thetaMin + thetaStep*_it ), "L" );
 
    Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

    Double_t _weight = 1.0/(energyStep*thetaStep*TMath::DegToRad()*_sin_theta*2*TMath::Pi()*_n_gen_events/integratedCrossSectionRad ) ;

    _h_temp_py -> Scale( _weight );

    // _h_temp_py -> Rebin( 5 );

    _h_temp_py -> SetMinimum( 0.00001  );

    _h_temp_py -> SetMaximum( 100   );

    // _c_gen -> cd( 1 );

    gPad -> SetLogy() ;

    if( !_it ) _h_temp_py -> Draw();
    else       _h_temp_py -> Draw("same");
  
  } /* generation check stops here */

  _leg -> Draw( );

  /* check the distribution of the cross-section */
  TCanvas* _c_rec = new TCanvas("c_java_rec_check", "", 1734, 844 );

  _c_rec -> cd();

  /* reconstructed cross-sections */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_rec_temp_py = ( hEnergyvsTheta -> ProjectionY( Form("h_rec_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    if( !_it ){

      _h_rec_temp_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

      _h_rec_temp_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_rec_temp_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_rec_temp_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_rec_temp_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_rec_temp_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_rec_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    }
    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    // cout << " Color set to " << _color << "\t" << _it << "\t" << _i_multiple << endl ;

    _h_rec_temp_py -> SetLineColor(   _color );

    _h_rec_temp_py -> SetMarkerColor( _color );

   // cout << _h_rec_temp_py -> GetEntries() << " \t " << _it << endl ;

    Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

    Double_t _weight = 1.0/(thetaStep*TMath::DegToRad()*energyStep*_sin_theta*2*TMath::Pi()*_n_gen_events/integratedCrossSectionRad ) ;

    _h_rec_temp_py -> Scale( _weight );

    // _h_rec_temp_py -> Rebin( 5 );

    _h_rec_temp_py -> SetMinimum( 0.0001  );

    _h_rec_temp_py -> SetMaximum( 100   );

    // _c_rec -> cd( 1 );

    gPad -> SetLogy() ;

    if( !_it ) _h_rec_temp_py -> Draw();
    else       _h_rec_temp_py -> Draw("same");

  } /* generation check stops here */

  /* check the distribution of the rates */
  TCanvas* _c_rates = new TCanvas("c_java_rates", "", 1734, 844 );

  _c_rates -> cd();

  /* reconstructed cross-sections */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_rate_temp_py = ( hEnergyvsTheta -> ProjectionY( Form("h_rate_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    if( !_it ){

      _h_rate_temp_py -> GetYaxis() -> SetTitle( Form ("Inclusive Electron Rate/%g GeV/%g deg/sec", energyStep, thetaStep ) );

      _h_rate_temp_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_rate_temp_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_rate_temp_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_rate_temp_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_rate_temp_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_rate_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    }

    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    _h_rate_temp_py -> SetLineColor(   _color );

    _h_rate_temp_py -> SetMarkerColor( _color );

    /* it comes from L=10^35 cm-2 s-1 x 10-30 = 10^5, that is the conversion of the cross-section from microbarn to cm^2*/
    Double_t _lumi_factor  = Luminosity*pow( 10, -30 ) ;

    Double_t _scale_factor = integratedCrossSectionRad/_n_gen_events*_lumi_factor ;

    _h_rate_temp_py -> Scale( _scale_factor ) ;

    _h_rate_temp_py -> SetMinimum( 0.0001  );

    _h_rate_temp_py -> SetMaximum( 100    );

    gPad -> SetLogy() ;

    if( !_it ) _h_rate_temp_py -> Draw();
    else       _h_rate_temp_py -> Draw("same");
  
  } /* generation check stops here */

  _leg -> Draw( );

  output_file -> WriteTObject( _c_rec   );

  output_file -> WriteTObject( _c_rates );

  output_file -> WriteTObject( _c_gen   );

}

/* ------------------------------------------- */
/*  Read data from Java - including the fastMC */
/* ------------------------------------------- */
void InclusiveRateAnalysis::ReadHadronJavaFiles( const char* _java_file_name, Int_t _n_java_files, Double_t _int_dis_cross_section )
{

  /* clear vectors to store values */
  ClearVectors();

  /* if this files are produced through clasDIS, the cross-section coming from it has to be used */
  // integratedDISCrossSection = 0.02733 /* tot cross-section@6GeV = 1.01242; */ /* tot cross-section @11 GeV: 0.370401 ; */
  integratedDISCrossSection = _int_dis_cross_section ; /* tot cross-section@6GeV = 1.01242; */ /* tot cross-section @11 GeV: 0.370401 ; */

  /* set parameter - cross-section - for PiP */
  // Double_t _PiPFraction = 0.840248 ;

  // Double_t _intCrossSectionPiP = _PiPFraction*integratedCrossSectionRad ;

  Int_t _n_gen_events = 0 ;

  Int_t _n_rec_hadrons[ n_hadron_types ] ;

  Int_t _n_gen_hadrons[ n_hadron_types ] ;

  // Int_t _n_rec_pip_events = 0 ;

  // Int_t _n_gen_pip_events = 0 ;

  ifstream _input_java_file ;

  for( Int_t _if = 0 ; _if < _n_java_files ; _if++ ){

    cout << " Now reading file " << Form( "%s%d.dat", _java_file_name, _if) << endl ;

    _input_java_file.open( Form( "%s%d.dat", _java_file_name, _if) );

    if (_input_java_file.is_open())
      {
	while ( _input_java_file.good() )
	  {

	    TString _event_type = "";
	    Int_t   _lund_id ;
	    Double_t _px, _py, _pz, _theta, _phi ;

	    _input_java_file >> _event_type >> _lund_id >> _px >> _py >> _pz >> _theta >> _phi ;

	    // cout << " event_type is " << _event_type.Data() << endl ;

	    // if( !_event_type.CompareTo( "gen" ) && _lund_id == 211 ) cout << _lund_id << "\t" << _px << "\t" << _py << "\t" << _pz << "\t" << _theta << "\t" << _phi << endl ;

	    /* check if the lund id is mapped - otherwise, attribute the generic number n_hadron_types + 1 */

	    if ( mapLundIdToNumber.find( _lund_id ) != mapLundIdToNumber.end() ) {

	      Int_t    _hadron_id   = mapLundIdToNumber[             _lund_id   ];

	      Double_t _hadron_mass = mapLundIdToMass[   TMath::Abs( _lund_id ) ];

	      // cout << " hadron found - mass is " << _hadron_mass << endl ;

	      // Int_t    _hadron_id   = 0 ;
	      // Double_t _hadron_mass = 0;

	      // if( _lund_id == 211 ){

	      //   _hadron_id = 0 ;
	      //   _hadron_mass = PiChMass ;

	      // } else if( _lund_id == -211 ) {

	      //   _hadron_id = 1 ;
	      //   _hadron_mass = PiChMass ;

	      // } else if( _lund_id == 321 ) {

	      //   _hadron_id = 2 ;
	      //   _hadron_mass = KaonChMass ;

	      // } else if( _lund_id == -321 ) {

	      //   _hadron_id = 3 ;
	      //   _hadron_mass = KaonChMass ;

	      // }

	      Double_t _energy = TMath::Sqrt( _px*_px + _py*_py + _pz*_pz + _hadron_mass*_hadron_mass );

	      /* check first the generated electrons */
	      if(      !_event_type.CompareTo( "gen" ) && _lund_id ==  11 ) { hEnergyvsThetaGen -> Fill( _theta, _energy ); hThetavsPhiGen -> Fill( _phi, _theta ); _n_gen_events++ ; /* cout << _event_type.Data() << "\t" << _lund_id << endl ; */ }
	      /* then the hadrons */
	      else if( !_event_type.CompareTo( "gen" )) { hGenHadronEnergyvsTheta[ _hadron_id ] -> Fill( _theta, _energy ); hGenHadronThetavsPhi[ _hadron_id ] -> Fill( _phi, _theta ); _n_gen_hadrons[ _hadron_id ]++ ; }
	      else if( !_event_type.CompareTo( "rec" )) { hRecHadronEnergyvsTheta[ _hadron_id ] -> Fill( _theta, _energy ); hRecHadronThetavsPhi[ _hadron_id ] -> Fill( _phi, _theta ); _n_rec_hadrons[ _hadron_id ]++ ; }

	      // if(      !_event_type.CompareTo( "gen" ) && _lund_id ==  11 ) { hEnergyvsThetaGen -> Fill( _theta, _energy ); hThetavsPhiGen -> Fill( _phi, _theta ); _n_gen_events++ ; /* cout << _event_type.Data() << "\t" << _lund_id << endl ; */ }
	      // else if( !_event_type.CompareTo( "gen" ) && _lund_id == 211 ) { hGenHadronEnergyvsTheta[ _hadron_id ] -> Fill( _theta, _energy ); hGenHadronThetavsPhi[ _hadron_id ] -> Fill( _phi, _theta ); _n_gen_pip_events++ ; }
	      // else if( !_event_type.CompareTo( "rec" ) && _lund_id == 211 ) { hRecHadronEnergyvsTheta[ _hadron_id ] -> Fill( _theta, _energy ); hRecHadronThetavsPhi[ _hadron_id ] -> Fill( _phi, _theta ); _n_rec_pip_events++ ; }

	      /* closes the condition on the lund id being mapped */
	    } else cout << " Hadron not found - lund id is " << _lund_id << endl ;

	  }
      }

    _input_java_file.close();

  } /* it closes the loop on the files */

  /* calculate the pip fraction */
  // Double_t _PiPFraction = ( Double_t )( (double)_n_gen_pip_events/(double)_n_gen_events ) ;

  // Double_t _intCrossSectionPiP = _PiPFraction*integratedCrossSectionRad ;

  // cout << " PiP Fraction is " << _PiPFraction << "\t" << _n_gen_events << "\t" << _n_gen_pip_events << endl ;

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  /* generated cross-sections for any hadron species */
  for( Int_t _ih = 0 ; _ih < n_hadron_types ; _ih++ ){
  // for( Int_t _ih = 0 ; _ih < 1 ; _ih++ ){

    TCanvas* _c_hadrons = new TCanvas( Form( "c_java_hadrons_%d", _ih ), "", 1734, 844 );

    _c_hadrons -> cd();

    for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

      Double_t _theta_now = thetaMin + thetaStep*_it ;

      // TH1D* _h_ele_py = ( hEnergyvsThetaGen -> ProjectionY( Form("h_ele_%d_gen_projY_theta_%g",    _ih, _theta_now ), _it + 1, _it + 1 ) ) ;

      TH1D* _h_gen_py = ( hGenHadronEnergyvsTheta[ _ih ] -> ProjectionY( Form("h_hadron_%d_gen_projY_theta_%g", _ih, _theta_now ), _it + 1, _it + 1 ) ) ;

      TH1D* _h_rec_py = ( hRecHadronEnergyvsTheta[ _ih ] -> ProjectionY( Form("h_hadron_%d_rec_projY_theta_%g", _ih, _theta_now ), _it + 1, _it + 1 ) ) ;

      //  if( !_it ){

        _h_gen_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

        // _h_gen_py -> GetXaxis() -> SetTitle( "E_{h} (GeV)" );
        _h_gen_py -> GetXaxis() -> SetTitle( "E_{#pi^{+}} (GeV)" );

        _h_gen_py -> GetXaxis() -> SetLabelFont(   132     );
	 	      
        _h_gen_py -> GetYaxis() -> SetLabelFont(   132     );
	 	      
        _h_gen_py -> GetXaxis() -> SetTitleFont(   132     );
	 	      
        _h_gen_py -> GetYaxis() -> SetTitleFont(   132     );
	 	      
      //   _h_gen_py -> GetZaxis() -> SetLabelFont(   132     );
	 	      
      //   _h_gen_py -> GetZaxis() -> SetTitleFont(   132     );

      // }

      Int_t _color = 1 + _it ;

      Int_t _i_multiple = (int)( _it/99 ) ;

      if( _color > 99 ) _color -= 99*_i_multiple ;

      _h_gen_py -> SetLineColor(   _color );

      _h_gen_py -> SetMarkerColor( _color );

      _h_rec_py -> SetLineColor(   _color );

      _h_rec_py -> SetMarkerColor( _color );

      // _leg -> AddEntry( _h_gen_py, Form( "#theta = %g", thetaMin + thetaStep*_it ), "L" );
 
      Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

      Double_t _weight = 1.0/(energyStep*thetaStep*TMath::DegToRad()*_sin_theta*2*TMath::Pi()*_n_gen_events/integratedDISCrossSection ) ;

      // _h_ele_py -> Scale( _weight );

      // _h_ele_py -> SetMinimum( 0.00001  );

      // _h_ele_py -> SetMaximum( 100   );

      _h_gen_py -> Scale( _weight );

      _h_gen_py -> SetMinimum( 0.0000001  );

      _h_gen_py -> SetMaximum( 1   );

      _h_rec_py -> Scale( _weight );

      // _h_rec_py -> SetMinimum( 0.00001  );

      // _h_rec_py -> SetMaximum( 100   );

      gPad -> SetLogy() ;

      if( !_it ) _h_gen_py -> Draw();
      else       _h_gen_py -> Draw("same");
  
      TCanvas* _c_comp_hadrons = new TCanvas( Form( "c_comp_%d_theta%g", _ih, _theta_now ) );

      gPad ->SetLogy( 1 );

      // _h_ele_py -> Draw();

      TLegend _l( 0.5, 0.6, 0.7, 0.8 );

      _l.SetFillStyle( 0 );

      _l.SetLineColor( 0 );

      _l.AddEntry( _h_gen_py, "generated",     "L" );

      _l.AddEntry( _h_rec_py, "reconstructed", "L" );

      _h_gen_py -> SetLineColor( 1 );

      _h_rec_py -> SetLineColor( kTeal + 1 );

      _h_gen_py -> SetStats(0);

      _h_rec_py -> SetStats(0);

      _h_gen_py -> Draw();

      _h_rec_py -> Draw("same");

      _l.Draw();

      output_file -> WriteTObject( _c_comp_hadrons );


    } /* generation check stops here */

  // _leg -> Draw( );

  output_file -> WriteTObject( _c_hadrons );

  } /* it closes the loop on the hadron types */

}

/* ------------------------------------------- */
/*  Read data from Java - including the fastMC */
/* ------------------------------------------- */
void InclusiveRateAnalysis::ReadJavaRecAndFastMC( const char* _java_file_name, Int_t _n_java_files )
{

  /* clear vectors to store values */
  ClearVectors();

  integratedCrossSectionRad = 1.01242; /* this is for 11 GeV: 0.370401 ; */

  Int_t _n_gen_events = 0 ;

  Int_t _n_rec_events = 0 ;

  Int_t _n_fMC_events = 0 ;

  ifstream _input_java_file ;

  for( Int_t _if = 0 ; _if < _n_java_files ; _if++ ){

    cout << " Now reading file " << Form( "%s%d.dat", _java_file_name, _if) << endl ;

    _input_java_file.open( Form( "%s%d.dat", _java_file_name, _if) );

    if (_input_java_file.is_open())
      {
	while ( _input_java_file.good() )
	  {

	    TString _event_type = "";
	    Double_t _px, _py, _pz, _theta, _phi ;

	    _input_java_file >> _event_type >> _px >> _py >> _pz >> _theta >> _phi ;

	    // cout << " event_type is " << _event_type.Data() << endl ;

	    Double_t _energy = TMath::Sqrt( _px*_px + _py*_py + _pz*_pz + ElectronMass*ElectronMass );

	    if(      !_event_type.CompareTo( "gen" )                   ) { energyEleGen.push_back( _energy ); thetaEleGen.push_back( _theta );  phiEleGen.push_back( _phi ); hEnergyGen -> Fill( _energy ); hEnergyvsThetaGen -> Fill( _theta, _energy ); hThetavsPhiGen -> Fill( _phi, _theta ); _n_gen_events++ ; }
	    else if( !_event_type.CompareTo( "fastmc" ) && _pz > -99.9 ) { hEnergyvsThetaFastMC    -> Fill( _theta, _energy ); hThetavsPhiFastMC    -> Fill( _phi, _theta ); _n_fMC_events++ ; }
	    else if( !_event_type.CompareTo( "rec" )    && _pz > -99.9 ) { hEnergyvsTheta          -> Fill( _theta, _energy ); hThetavsPhi          -> Fill( _phi, _theta ); 
	    _n_rec_events++ ; }
	    // else       cout << " WARNING! Something is wrong with the format of the java dat file! " << endl ;

	  }
      }

    _input_java_file.close();

  } /* it closes the loop on the files */


  // cout << Form( "Found %d gen events, %d in the new fastMC and %d in rec", _n_gen_events, _n_fMC_events, _n_rec_events ) << endl ;

  cout << Form( "Analyzed %d gen events; found %d rec events and %d fastMC events", _n_gen_events, _n_rec_events, _n_fMC_events ) << endl ;

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  // TLegend* _leg = new TLegend( 0.91, 0.16, 0.99, 0.87 );

  // _leg -> SetFillStyle( 0 );

  // _leg -> SetLineColor( 0 );

  // _leg -> Draw( );

  /* check the distribution of the cross-section */
  TCanvas* _c_fastMC = new TCanvas("c_java_fastMC", "", 1734, 844 );

  _c_fastMC -> cd();

  /* reconstructed cross-sections */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_fastMC_temp_py = ( hEnergyvsThetaFastMC -> ProjectionY( Form("h_fastMC_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    if( !_it ){

      _h_fastMC_temp_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

      _h_fastMC_temp_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_fastMC_temp_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_fastMC_temp_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_fastMC_temp_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_fastMC_temp_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_fastMC_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    }
    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    // cout << " Color set to " << _color << "\t" << _it << "\t" << _i_multiple << endl ;

    _h_fastMC_temp_py -> SetLineColor(   _color );

    _h_fastMC_temp_py -> SetMarkerColor( _color );

   // cout << _h_fastMC_temp_py -> GetEntries() << " \t " << _it << endl ;

    Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

    Double_t _weight = 1.0/(thetaStep*TMath::DegToRad()*energyStep*_sin_theta*2*TMath::Pi()*_n_gen_events/integratedCrossSectionRad ) ;

    _h_fastMC_temp_py -> Scale( _weight );

    // _h_fastMC_temp_py -> Rebin( 5 );

    _h_fastMC_temp_py -> SetMinimum( 0.0001  );

    _h_fastMC_temp_py -> SetMaximum( 100   );

    // _c_fastMC -> cd( 1 );

    gPad -> SetLogy() ;

    if( !_it ) _h_fastMC_temp_py -> Draw();
    else       _h_fastMC_temp_py -> Draw("same");

  } /* generation check stops here */

  /* check the distribution of the fastMC_rates */
  TCanvas* _c_fastMC_rates = new TCanvas("c_java_fastMC_rates", "", 1734, 844 );

  _c_fastMC_rates -> cd();

  /* reconstructed cross-sections */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_rate_fastmc_py = ( hEnergyvsThetaFastMC -> ProjectionY( Form("h_fastMC_rate_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    if( !_it ){

      _h_rate_fastmc_py -> GetYaxis() -> SetTitle( Form( "Inclusive Electron Rate/%g GeV/%g deg/sec", energyStep, thetaStep ) );

      _h_rate_fastmc_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_rate_fastmc_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_rate_fastmc_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_rate_fastmc_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_rate_fastmc_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_rate_fastmc_py -> GetZaxis() -> SetTitleFont(   132     );

    }

    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    _h_rate_fastmc_py -> SetLineColor(   _color );

    _h_rate_fastmc_py -> SetMarkerColor( _color );

    Double_t _lumi_factor  = Luminosity*pow( 10, -30 ) ;

    Double_t _scale_factor = integratedCrossSectionRad/_n_gen_events*_lumi_factor ;

    _h_rate_fastmc_py -> Scale( _scale_factor ) ;

    // _h_rate_fastmc_py -> SetMinimum( 0.00001  );

    // _h_rate_fastmc_py -> SetMaximum( 1    );

    gPad -> SetLogy() ;

    if( !_it ) _h_rate_fastmc_py -> Draw();
    else       _h_rate_fastmc_py -> Draw("same");
  
  } /* generation check stops here */

  output_file -> WriteTObject( _c_fastMC   );

  output_file -> WriteTObject( _c_fastMC_rates );

}

/* ------------------------------------------- */
/*  Compare Generated Events to original ones  */
/* ------------------------------------------- */
void InclusiveRateAnalysis::CompareJavaEventsToGen( const char* _gen_root_file, const char* _java_root_file )
{

  TFile* _gen_file  = new TFile( _gen_root_file  );  gDirectory -> cd(); 

  TFile* _java_file = new TFile( _java_root_file );  gDirectory -> cd(); 

  TH2D* _h_e_theta  = (TH2D*)( _gen_file -> Get( "hEnergyvsThetaGen" ) );

  Int_t _n_2d_bins = _h_e_theta -> GetNbinsX();

  /* get cross-section histograms */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _actualTheta = thetaMin + thetaStep*_it ;

    /* generator curve */
    TGraph* _ge_gen  = (TGraph*)( _gen_file -> Get( Form( "ge_sigma_rad_theta_%g", _actualTheta ) ) );

    /* generated events from java bank */
    TH1D*   _h_java  = (TH1D*)(   _java_file -> Get( Form( "h_java_projY_theta_%g", _actualTheta ) ) );

    /* reconstructed events from java bank */
    TH1D*   _h_rec   = (TH1D*)(   _java_file -> Get( Form( "h_rec_projY_theta_%g",  _actualTheta ) ) );

    if( _ge_gen == NULL || _h_java == NULL || _h_rec == NULL ) cout << " WARNING! Histo(s) not found!" << endl ; 

    // if( !_it ){

    //   _ge_gen -> SetLineColor(   2 );

    //   _ge_gen -> SetMarkerColor( 2 );

    // } else {

    _ge_gen -> SetLineColor(   1 );

    _ge_gen -> SetMarkerColor( 1 );

    _h_java -> SetLineColor(   kViolet + 1 );

    _h_java -> SetMarkerColor( kViolet + 1 );

    _h_rec  -> SetLineColor(   kTeal   + 1 );

    _h_rec  -> SetMarkerColor( kTeal   + 1 );

    // }

    /* draw the two superimposed in a canvas */
    TCanvas* _c_temp = new TCanvas( Form( "c_compare_java_theta_%g", _actualTheta ), "", 600, 600 );

    _c_temp -> cd();

    gPad    -> SetLogy( 1 );

    _h_java -> SetStats( 0 ) ;

    _h_java -> Draw();

    _h_rec  -> Draw("same");

    _ge_gen -> Draw("same");

    TLegend _leg( 0.3, 0.65, 0.5, 0.85 );

    _leg.SetFillStyle( 0 );

    _leg.SetLineColor( 0 );

    _leg.SetTextFont( 132 );

    _leg.AddEntry( _ge_gen, "cross-section calculation", "L" );

    _leg.AddEntry( _h_java, "generated",                 "L" );

    _leg.AddEntry( _h_rec , "reconstructed",             "L" );

    _leg.Draw();

    TPaveText* _pt_temp = new TPaveText( 6, 10, 9, 70 );
    // pt->AddLine(.0,.5,1.,.5);

    _pt_temp -> SetFillColor(  0 );

    _pt_temp -> SetTextSize( 0.032 );

    _pt_temp -> SetTextFont( 132 );

    _pt_temp -> AddText( Form( "Comparison for theta = %g (deg)", _actualTheta ) );

    _pt_temp -> Draw( "same" );

    output_file -> WriteTObject( _h_java );

    output_file -> WriteTObject( _h_rec  );

    output_file -> WriteTObject( _ge_gen );

    output_file -> WriteTObject( _c_temp );

  } /* it closes the loop on the theta bins */

}

/* ------------------------------------------- */
/*  Compare rates from the fastMC to JAVA      */
/* ------------------------------------------- */
void InclusiveRateAnalysis::CompareFastMCToJava( const char* _fastmc_root_file, const char* _java_root_file )
{

  TFile* _fastmc_file  = new TFile( _fastmc_root_file  );  gDirectory -> cd(); 

  TFile* _java_file = new TFile( _java_root_file );  gDirectory -> cd(); 

  TH2D* _h_e_theta  = (TH2D*)( _fastmc_file -> Get( "hEnergyvsThetaGen" ) );

  Int_t _n_2d_bins = _h_e_theta -> GetNbinsX();

  /* get cross-section histograms */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _actualTheta = thetaMin + thetaStep*_it ;

    /* generator curve */
    TGraph* _ge_fastmc  = (TGraph*)( _fastmc_file -> Get( Form( "ge_rate_theta_%g", _actualTheta ) ) );

    /* reconstructed events from java bank */
    TH1D*   _h_rec   = (TH1D*)(   _java_file -> Get( Form( "h_rate_projY_theta_%g",  _actualTheta ) ) );

    if( _ge_fastmc == NULL || _h_rec == NULL ) cout << " WARNING! Histo(s) not found!" << endl ; 

    // if( !_it ){

    //   _ge_fastmc -> SetLineColor(   2 );

    //   _ge_fastmc -> SetMarkerColor( 2 );

    // } else {

    _ge_fastmc -> SetLineColor(   1 );

    _ge_fastmc -> SetMarkerColor( 1 );

    _h_rec  -> SetLineColor(   kTeal   + 1 );

    _h_rec  -> SetMarkerColor( kTeal   + 1 );

    // }

    /* draw the two superimposed in a canvas */
    TCanvas* _c_temp = new TCanvas( Form( "c_compare_rate_theta_%g", _actualTheta ), "", 600, 600 );

    _c_temp -> cd();

    gPad    -> SetLogy( 1 );

    _h_rec -> SetStats( 0 ) ;

    _h_rec -> GetYaxis() -> SetTitle( Form( "Inclusive Electron Rate/%g GeV/%g deg/sec", energyStep, thetaStep ) );

    _h_rec -> GetXaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetYaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetZaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetXaxis() -> SetTitleFont(   132     );
		      
    _h_rec -> GetYaxis() -> SetTitleFont(   132     );
		      
    _h_rec -> GetZaxis() -> SetTitleFont(   132     );

    _h_rec -> Draw();

    _ge_fastmc -> Draw("same");

    TLegend _leg( 0.3, 0.65, 0.5, 0.85 );

    _leg.SetFillStyle( 0 );

    _leg.SetLineColor( 0 );

    _leg.SetTextFont( 132 );

    _leg.AddEntry( _ge_fastmc, "rate from fastMC",    "L" );

    _leg.AddEntry( _h_rec,     "rate from gemc+java", "L" );

    _leg.Draw();

    TPaveText* _pt_temp = new TPaveText( 6, 10, 9, 70 );
    // pt->AddLine(.0,.5,1.,.5);

    _pt_temp -> SetFillColor(  0 );

    _pt_temp -> SetTextSize( 0.032 );

    _pt_temp -> SetTextFont( 132 );

    _pt_temp -> AddText( Form( "Comparison for theta = %g (deg)", _actualTheta ) );

    _pt_temp -> Draw( "same" );

    output_file -> WriteTObject( _h_rec  );

    output_file -> WriteTObject( _ge_fastmc );

    output_file -> WriteTObject( _c_temp );

  } /* it closes the loop on the theta bins */

}

/* ------------------------------------------------------------- */
/*  Compare rates from the fastMC to JAVA and to the NEW fastMC  */
/* ------------------------------------------------------------- */
void InclusiveRateAnalysis::CompareNewFastMCToJava( const char* _fastmc_root_file, const char* _java_root_file, const char* _new_fastmc_root_file )
{

  TFile* _fastmc_file      = new TFile( _fastmc_root_file      );  gDirectory -> cd(); 

  TFile* _new_fastmc_file  = new TFile( _new_fastmc_root_file  );  gDirectory -> cd(); 

  TFile* _java_file        = new TFile( _java_root_file        );  gDirectory -> cd(); 

  TH2D* _h_e_theta  = (TH2D*)( _fastmc_file -> Get( "hEnergyvsThetaGen" ) );

  Int_t _n_2d_bins = _h_e_theta -> GetNbinsX();

  /* get cross-section histograms */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _actualTheta = thetaMin + thetaStep*_it ;

    /* generator curve */
    TGraph* _ge_fastmc  = (TGraph*)( _fastmc_file -> Get( Form( "ge_rate_theta_%g", _actualTheta ) ) );

    /* reconstructed events from java fastMC */
    TH1D*   _h_fastmc   = (TH1D*)(   _new_fastmc_file -> Get( Form( "h_fastMC_rate_projY_theta_%g",  _actualTheta ) ) );

    /* reconstructed events from java bank */
    TH1D*   _h_rec      = (TH1D*)(   _java_file -> Get( Form( "h_rate_projY_theta_%g",  _actualTheta ) ) );

    if( _ge_fastmc == NULL || _h_fastmc == NULL || _h_rec == NULL ) cout << " WARNING! Histo(s) not found!" << endl ; 

    // if( !_it ){

    //   _ge_fastmc -> SetLineColor(   2 );

    //   _ge_fastmc -> SetMarkerColor( 2 );

    // } else {

    _ge_fastmc -> SetLineColor(   1 );

    _ge_fastmc -> SetMarkerColor( 1 );

    _h_rec     -> SetLineColor(   kTeal   + 1 );

    _h_rec     -> SetMarkerColor( kTeal   + 1 );

    _h_fastmc  -> SetLineColor(   kViolet   + 1 );

    _h_fastmc  -> SetMarkerColor( kViolet   + 1 );

    // }

    /* draw the two superimposed in a canvas */
    TCanvas* _c_temp = new TCanvas( Form( "c_compare_rate_theta_%g", _actualTheta ), "", 600, 600 );

    _c_temp -> cd();

    gPad    -> SetLogy( 1 );

    _h_rec -> SetStats( 0 ) ;

    _h_rec -> SetMinimum( 0.0001 );

    _h_rec -> SetMaximum( 100    );

    _h_rec -> GetYaxis() -> SetTitle( Form( "Inclusive Electron Rate/%g GeV/%g deg/sec", energyStep, thetaStep ) );

    _h_rec -> GetXaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetYaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetZaxis() -> SetLabelFont(   132     );
		      
    _h_rec -> GetXaxis() -> SetTitleFont(   132     );
		      
    _h_rec -> GetYaxis() -> SetTitleFont(   132     );
		      
    _h_rec -> GetZaxis() -> SetTitleFont(   132     );

    _h_rec -> Draw();

    _ge_fastmc -> Draw("same");

    _h_fastmc  -> Draw("same");

    TLegend _leg( 0.3, 0.65, 0.5, 0.85 );

    _leg.SetFillStyle( 0 );

    _leg.SetLineColor( 0 );

    _leg.SetTextFont( 132 );

    _leg.AddEntry( _ge_fastmc, "rate from fastMC (clasev)", "L" );

    _leg.AddEntry( _h_fastmc,  "rate from fastMC (java)",   "L" );

    _leg.AddEntry( _h_rec,     "rate from gemc+java",       "L" );

    _leg.Draw();

    TPaveText* _pt_temp = new TPaveText( 6, 10, 9, 70 );
    // pt->AddLine(.0,.5,1.,.5);

    _pt_temp -> SetFillColor(  0 );

    _pt_temp -> SetTextSize( 0.032 );

    _pt_temp -> SetTextFont( 132 );

    _pt_temp -> AddText( Form( "Comparison for theta = %g (deg)", _actualTheta ) );

    _pt_temp -> Draw( "same" );

    output_file -> WriteTObject( _h_rec  );

    output_file -> WriteTObject( _h_fastmc );

    output_file -> WriteTObject( _ge_fastmc );

    output_file -> WriteTObject( _c_temp );

  } /* it closes the loop on the theta bins */

}


/* ------------------------------------------- */
/*  Compare Generated Events to original ones  */
/* ------------------------------------------- */
void InclusiveRateAnalysis::CompareGenEventsToOriginal()
{

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  /* get cross-section histograms */
  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _actualTheta = thetaMin + thetaStep*_it ;

    TGraph* _ge_temp = (TGraph*)( output_file -> Get( Form( "ge_sigma_rad_theta_%g", _actualTheta ) ) );

    TH1D*   _h_temp  = (TH1D*)(   output_file -> Get( Form("h_projY_theta_%g",       _actualTheta ) ) );

    if( _ge_temp == NULL || _h_temp == NULL ) cout << Form( " WARNING! %s not found! Check if you are properly storing it in the output file \n", _h_temp -> GetName() );

    else {

      if( !_it ){

	_ge_temp -> SetLineColor(   2 );

	_ge_temp -> SetMarkerColor( 2 );

      } else {

	_ge_temp -> SetLineColor(   1 );

	_ge_temp -> SetMarkerColor( 1 );

      }

      /* draw the two superimposed in a canvas */
      TCanvas* _c_temp = new TCanvas( Form( "c_compare_theta_%g", _actualTheta ), "", 1200, 800 );

      _c_temp  -> cd();

      _h_temp  -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

      // Int_t _max_bin = (int)( _h_temp -> GetMaximumBin() );

      // Int_t _min_bin = (int)( _h_temp -> GetMinimumBin() );

      // Double_t _y_max = _h_temp  -> GetBinContent( _max_bin ) ; 

      // Double_t _y_min = _h_temp  -> GetBinContent( _min_bin ) ; 

      // _h_temp  -> GetYaxis() -> SetRangeUser( _y_min*1.1, _y_max*1.1 ) ;

      _h_temp  -> SetStats( 0 ) ;

      _h_temp  -> Draw();

      _ge_temp -> Draw("same");

      gPad     -> SetLogy( 1 );

      TPaveText* _pt_temp = new TPaveText( 3, 10, 5, 80 );
      // pt->AddLine(.0,.5,1.,.5);

      _pt_temp -> SetFillColor(  0 );

      _pt_temp -> SetTextSize( 0.032 );

      _pt_temp -> SetTextFont( 132 );

      _pt_temp -> AddText( Form( "Comparison for theta = %g (deg)", _actualTheta ) );

      _pt_temp -> Draw( "same" );

      output_file -> WriteTObject( _c_temp );

    } /* it closes the condition on having the ge and the histo */

  } /* it closes the loop on the theta bins */

}

// ----------------------
//  Generate pseudo data
// ----------------------

void InclusiveRateAnalysis::GeneratePseudoData( Int_t _n_events, Int_t _n_events_per_lund_file, const char* _dir_for_lund_files )
{

  if( _n_events < 0 ) _n_events = fileNLines ;

  Int_t _n_events_per_lundfile = _n_events_per_lund_file ;

  /* only store 20k events per lundfile */
  Int_t _n_lund_files = TMath::Floor( _n_events/_n_events_per_lundfile ) ;

  if( _n_events % _n_events_per_lundfile ) _n_lund_files += 1 ; /* to account for the residual events */

  cout << " _n_lund_files is " << _n_lund_files << "\n" << endl ;

  /* out lundfile */
  ofstream _out_lund_file[ _n_lund_files ] ;

  for( Int_t _ilf = 0 ; _ilf < _n_lund_files; _ilf++ ){ _out_lund_file[ _ilf ].open( Form( "%s/out_lund_files_n%d.dat", _dir_for_lund_files, _ilf ) );  cout <<"file: " << Form( "%s/out_lund_files_n%d.dat", _dir_for_lund_files, _ilf ) <<endl;
    cout  <<"################### file directory: " << Form( "%s", _dir_for_lund_files) <<endl;
  }

  /* produce _n_events simulated data */
  for( Int_t _iev = 0 ; _iev < _n_events ; _iev++  ){

    if( _iev && _iev % _n_events_per_lundfile == 0 ) cout << " Generated " << _iev << " events " << endl ; 
    
    Int_t _file_number = TMath::Floor( _iev/_n_events_per_lundfile ) ;
    if( _iev && _iev % _n_events_per_lundfile == 0 ) cout << "file number : " <<  _file_number <<endl;
    Double_t _theta  = 0 ;
    Double_t _energy = 0 ;

    Double_t _phi ;

    hEnergyvsTheta      -> GetRandom2( _theta, _energy );

    /* assign a spread in theta/energy */
    TRandom3 _rGenEnergySpread ;

    TRandom3 _rGenThetaSpread  ;

    _rGenEnergySpread.SetSeed( 0 ) ;
                       
    _rGenThetaSpread.SetSeed(  0 ) ;

    /* fill the histo with the generated data */
    hEnergyvsThetaGen      -> Fill(       _theta, _energy );

    hEnergyGen             -> Fill(               _energy );

    hThetaGen              -> Fill(               _theta  );

    TRandom3 _rGen;

    _rGen.SetSeed(0);

    /* phi generated in degrees */
    _phi = _rGen.Uniform( 0.0, 360.0 );

    Double_t _p = TMath::Sqrt( _energy*_energy - ElectronMass*ElectronMass );

    /* set four-momentum components */
    Double_t _px = _p*TMath::Sin( _theta*TMath::DegToRad() )*TMath::Cos( _phi*TMath::DegToRad() ) ;
    Double_t _py = _p*TMath::Sin( _theta*TMath::DegToRad() )*TMath::Sin( _phi*TMath::DegToRad() ) ;
    Double_t _pz = _p*TMath::Cos( _theta*TMath::DegToRad() ) ;

    TLorentzVector Pbeam, Pout ;
    Pbeam.SetPxPyPzE(   0.0,  0.0, energyBeam, energyBeam );
    Pout.SetPxPyPzE(    _px,  _py,        _pz,    _energy );

    Double_t _q2, _xb, _w, _y ;

    DefDISVar(Pbeam, Pout, _q2, _xb, _w, _y );
    

    Double_t _nu = energyBeam - _energy ;

    hQ2vsWGen      -> Fill(   _w,    _q2 );

    hWGen          -> Fill(   _w         );

    hThetaGen      -> Fill(   _theta     );

    hThetavsPhiGen -> Fill( _phi, _theta );
    
    // hQ2Gen         -> Fill(   _q2        );

    _out_lund_file[ _file_number ] << "         1" << "  " <<  "1.0" << "  " <<  "1.0" << "   " <<  "0" << "   " <<  "1" << "   " <<  _xb << "   " <<  _y << "   " <<  _w*_w << "   " <<  _q2 << "   " <<  _nu <<  endl ;

    _out_lund_file[ _file_number ] << "  " << 1 << "   " <<  "-1.0" << "   " <<  "1" << "   " <<  11 << "     " <<  0 << "     " <<  0 << "   " <<   _px << "   " <<   _py << "   " <<  _pz << "   " <<  _energy << "   " <<  ElectronMass << "   " <<  "0.0000" << "   " <<   "0.0000" << "   " <<   "0.0000" <<  endl ;
    
  } /* it closes the loop on the number of generated events */

  for( Int_t _ilf = 0 ; _ilf < _n_lund_files; _ilf++ ){ _out_lund_file[ _ilf ].close() ; }

  /* check the distribution of the generated energy vs theta dependence */
  TCanvas* _c_theta_energy = new TCanvas("c_theta_energy", "", 1734, 844 );

  _c_theta_energy   -> Divide( 2 );

  _c_theta_energy   -> cd( 1 );

  hEnergyvsTheta    -> Draw( "colz" );

  _c_theta_energy   -> cd( 2 );

  hEnergyvsThetaGen -> Draw( "colz" );

  /* check the distribution of the generated cross-section */
  TCanvas* _c_gen = new TCanvas("c_gen_check", "", 1000, 700 );

  gPad -> SetLogy() ;

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  TLegend* _leg = new TLegend( 0.91, 0.06, 0.99, 0.97 );

  _leg -> SetFillStyle( 0 );

  _leg -> SetLineColor( 0 );

  /* check the total number of events */
  Int_t _n_theta_event_tot = 0 ;

  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    Double_t _theta_now = thetaMin + thetaStep*_it ;

    TH1D* _h_temp_py = ( hEnergyvsThetaGen -> ProjectionY( Form("h_projY_theta_%g", _theta_now ), _it + 1, _it + 1 ) ) ;

    if( !_it ){

      _h_temp_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

      _h_temp_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_temp_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    }

    _n_theta_event_tot += _h_temp_py -> GetEntries();

    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    _h_temp_py -> SetStats(            0 );

    _h_temp_py -> SetLineColor(   _color );

    _h_temp_py -> SetMarkerColor( _color );

    _leg -> AddEntry( _h_temp_py, Form( "#theta = %g", thetaMin + thetaStep*_it ), "L" );
 
    Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

    Double_t _weight = 1.0/(energyStep*thetaStep*TMath::DegToRad()*_sin_theta*2*TMath::Pi()*_n_events/integratedCrossSectionRad ) ;

    _h_temp_py -> Scale( _weight );

    _h_temp_py -> SetMinimum( 0.00001  );

    _h_temp_py -> SetMaximum( 1000   );

    if( !_it ) _h_temp_py -> Draw();
    else       _h_temp_py -> Draw("same");

  } /* generation check stops here */

  _leg -> Draw( );

  output_file -> WriteTObject( _c_gen          );

  output_file -> WriteTObject( _c_theta_energy );

}


// ----------------------
//  Generate pseudo data
// ----------------------

void InclusiveRateAnalysis::TestGeneratedLund( Int_t _n_lund_files, const char* _lund_file_body )
{

  Int_t _n_events = 0 ;

  /* in lundfile */
  ifstream _input_lund_file[ _n_lund_files ] ;

  for( Int_t _ilf = 0 ; _ilf < _n_lund_files; _ilf++ ){

    cout << Form( "Analyzing lund file %s%d.dat", _lund_file_body, _ilf ) << endl ;

    _input_lund_file[ _ilf ].open( Form( "%s%d.dat", _lund_file_body, _ilf ) );

    if ( _input_lund_file[ _ilf ].is_open())
      {

	Int_t _n_line = 0 ;

	while ( _input_lund_file[ _ilf ].good() )
	  {

	    Double_t _c1, _c2, _c3, _c4, _c5, _c6, _c7, _c8, _c9, _c10, _c11;

	    Double_t _px, _py, _pz, _energy, _m, _vx, _vy, _vz, _xb, _y, _w2, _q2, _nu  ;

	    // const char* _temp ;

	    if( _n_line % 2 == 0) _input_lund_file[ _ilf ] >> _c1 >> _c2 >> _c3 >> _c4 >> _c5 >> _xb  >>  _y >>  _w2 >>  _q2 >>  _nu ;
	    else                  _input_lund_file[ _ilf ] >> _c6 >> _c7 >> _c8 >> _c9 >> _c10 >> _c11 >> _px  >>  _py >>  _pz >>  _energy >> _m >> _vx >> _vy >> _vz ;

	    //        1  1.0  1.0   0   1   1.09761   0.0431907   0.793329   0.978567   0.475098
	    // 1   1.0   1   11     0     0   -0.870304   -0.420589   10.4804   10.5249   0.000510999   0.0000   0.0000   0.0000

	    // if( _n_line < 20 ){

	    //   cout << _n_line << endl ;

	    //   if( _n_line % 2 == 0) cout << "\t" << _c1 << "\t" << _c2 << "\t" << _c3 << "\t" << _c4 << "\t" << _c5 << "\t" << _xb  << "\t" <<  _y << "\t" <<  _w2 << "\t" <<  _q2 << "\t" <<  _nu << "\t" <<  endl;
	    //   else                  cout << "\t" << _c6 << "\t" << _c7 << "\t" << _c8 << "\t" << _c9 << "\t" << _c10 << "\t" << _c11 << "\t" << _px  << "\t" <<  _py << "\t" <<  _pz << "\t" <<  _energy << "\t" << _m << "\t" << _vx << "\t" << _vy << "\t" << _vz << endl;

	    // }

	    /* fill vectors once per event, i.e. only for odd lines */
	    if( _n_line % 2 != 0){

	      TLorentzVector _temp_4v ;

	      _temp_4v.SetPxPyPzE( _px, _py, _pz, _energy );

	      hEnergyvsTheta -> Fill( _temp_4v.Theta()*TMath::RadToDeg(), _temp_4v.E() );

	      hEnergy        -> Fill( _energy ) ;

	      _n_events++ ;

	    }

	    _n_line++ ;

	  }

	_n_line = 0 ;

      } /* close the opening of the file */

    _input_lund_file[ _ilf ].close();

  } /* it closes the loop on the files */

  // for( Int_t _ilf = 0 ; _ilf < _n_lund_files; _ilf++ ){ _out_lund_file[ _ilf ].close() ; }

  /* check the distribution of the generated cross-section */
  TCanvas* _c_gen = new TCanvas("c_lund_check", "", 1000, 700 );

  integratedCrossSectionRad = 1.01242; /* this is for 11 GeV: 0.370401 ; */

  gPad -> SetLogy() ;

  Int_t _n_2d_bins = hEnergyvsThetaGen -> GetNbinsX();

  TLegend* _leg = new TLegend( 0.91, 0.06, 0.99, 0.97 );

  _leg -> SetFillStyle( 0 );

  _leg -> SetLineColor( 0 );

  for( Int_t _it = 0 ; _it < _n_2d_bins ; _it++ ){

    TH1D* _h_temp_py = ( hEnergyvsTheta -> ProjectionY( Form("h_lund_projY_thetaBin_%d", _it ), _it + 1, _it + 1 ) ) ;

    _h_temp_py -> SetLineColor( _it + 1 );

    if( !_it ){

      _h_temp_py -> GetYaxis() -> SetTitle( "d#sigma/dEd#Omega #mu b/GeV/str" );

      _h_temp_py -> GetXaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetYaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetZaxis() -> SetLabelFont(   132     );
		      
      _h_temp_py -> GetXaxis() -> SetTitleFont(   132     );
		      
      _h_temp_py -> GetYaxis() -> SetTitleFont(   132     );
		      
      _h_temp_py -> GetZaxis() -> SetTitleFont(   132     );

    }

    Int_t _color = 1 + _it ;

    Int_t _i_multiple = (int)( _it/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    // cout << " Color set to " << _color << "\t" << _it << "\t" << _i_multiple << endl ;

    _h_temp_py -> SetLineColor(   _color );

    _h_temp_py -> SetMarkerColor( _color );


    // cout << _h_temp_py -> GetEntries() << " \t " << _it << endl ;

    _leg -> AddEntry( _h_temp_py, Form( "#theta = %g", thetaMin + thetaStep*_it ), "L" );
 
    Double_t _theta_now = thetaMin + thetaStep*_it ;

    Double_t _sin_theta = TMath::Sin( _theta_now*TMath::DegToRad() ) ;

    Double_t _weight = 1.0/(thetaStep*TMath::DegToRad()*energyStep*_sin_theta*2*TMath::Pi()*_n_events/integratedCrossSectionRad ) ;

    _h_temp_py -> Scale( _weight );

    // _h_temp_py -> Rebin( 5 );

    _h_temp_py -> SetMinimum( 0.00001  );

    _h_temp_py -> SetMaximum( 1000   );

    if( !_it ) _h_temp_py -> Draw();
    else       _h_temp_py -> Draw("same");

  } /* generation check stops here */

  _leg -> Draw( );

  output_file -> WriteTObject( _c_gen          );

}


/* -------------------------- */
/*    Set kinematical limits  */
/* -------------------------- */

void InclusiveRateAnalysis::SetKinematicLimits( Double_t _theta_min, Double_t _theta_max, Double_t _theta_step, Double_t _energy_min, Double_t _energy_max, Double_t _energy_step, Int_t _file_n_lines, Double_t _e_beam, Int_t _histos_n_bins  )
{

  thetaMin  = _theta_min  ;
  thetaMax  = _theta_max  ;
  thetaStep = _theta_step ;

  thetaPointN  = (int)( ( _theta_max - _theta_min )/( _theta_step ) ) + 1 ;

  energyMin  = _energy_min  ;
  energyMax  = _energy_max  ;
  energyStep = _energy_step ;

  energyPointN  = (int)( ( _energy_max - _energy_min )/( _energy_step ) ) + 1 ;

  sigmaMin   = 0.000000001 ;

  sigmaMax   = 1.0 ;

  fileNLines = _file_n_lines ;

  histosNBins = _histos_n_bins ;

  energyBeam = _e_beam ;

  cout << Form( "\n Analyzing the range %g - %g in step of %g, with Npoints=%d and E_{beam} = %g \n", thetaMin,  thetaMax,  thetaStep,  thetaPointN, energyBeam  ) << endl ; 

  cout << Form( "\n Analyzing the range %g - %g in step of %g, with Npoints=%d and E_{beam} = %g \n", energyMin, energyMax, energyStep, energyPointN, energyBeam ) << endl ; 

  return;

}


/* --------------------------------------------- */
/*   read file data - format of the probe system */
/* --------------------------------------------- */

void InclusiveRateAnalysis::ReadDataFile(){

  ifstream _input_file ;

  _input_file.open( file_name );

  /* clear vectors to store values */

  ClearVectors();

  Int_t _n_line = 0 ;

  if (_input_file.is_open())
    {
      while ( _input_file.good() )
  	{

	  Double_t _c1, _c2, _c3, _c4, _c5, _c6, _c7, _c8 ;

	  _input_file >> _c1 >> _c2 >> _c3 >> _c4 >> _c5 >> _c6 >> _c7 >> _c8 ;

	  thetaEle.push_back(      _c1 ) ;
	  energyEle.push_back(     _c2 ) ;
	  dsdThetadE.push_back(    _c3 ) ;   
	  dsdThetadErad.push_back( _c4 ) ;
	  Q2.push_back(            _c5 ) ;            
	  W.push_back(             _c6 ) ;           
	  dsdWdQ2.push_back(       _c7 ) ;      
	  dsdWdQ2rad.push_back(    _c8 ) ;   

	  /* fill this histo with a weight given by the cross-section times dE dOmega = dE 2 pi sinTheta dTheta - being the other terms (dE 2 pi dTheta) constant, we only consider sinTheta, since it is the only one affecting the distribution - we just want to have the relative weight, not the absolute cross-section value  */
	  Double_t _weight      = _c4*TMath::Sin( _c1*TMath::DegToRad() );

	  Double_t _weight_norm = _weight*thetaStep*TMath::DegToRad()*energyStep*2*TMath::Pi() ;

	  dsdThetadEradVolume.push_back( _weight_norm );

	  hEnergyvsTheta      -> Fill( _c1, _c2, _weight_norm );

	  hQ2vsW              -> Fill( _c6, _c5      );

	  hW                  -> Fill( _c6           );

	  hTheta              -> Fill( _c1           );

	  hEnergy             -> Fill( _c2           );

	  // if( _c1 == 10 && _c2 == 0.505 ) cout << _c1 << "\t" << _c2 << "\t" << _c4 << "\t" <<_weight_norm << " histo = " << hEnergyvsTheta -> GetBinContent( 11, 2 ) << endl ;
	  // if( _c1 == 7.05 || _c1 == 15.2 || _c1 == 15.7 ) cout << _c1 << "\t" << _c2 << "\t" << _c4 << "\t" <<_weight_norm << " histo = " << hEnergyvsTheta -> GetBinContent( 11, 2 ) << endl ;

	  _n_line++;

	}

      _input_file.close();

    }

  // cout << " Tot xs and xs vol = " << 

}

void InclusiveRateAnalysis::BookHistos()
{

  output_file -> cd();

  hEnergyvsThetaGen         = new TH2D( "hEnergyvsThetaGen"        , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax ) ;
  hEnergyvsTheta            = new TH2D( "hEnergyvsTheta"           , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax ) ;
  hEnergyvsThetaFastMC      = new TH2D( "hEnergyvsThetaFastMC"     , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax ) ;

  // hEnergyvsThetaGen         = new TH2D( "hEnergyvsThetaGen"        , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin - thetaStep*0.5, thetaMax + thetaStep*0.5, energyPointN, energyMin - energyStep*0.5, energyMax + energyStep*0.5 ) ;
  // hEnergyvsTheta            = new TH2D( "hEnergyvsTheta"           , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin - thetaStep*0.5, thetaMax + thetaStep*0.5, energyPointN, energyMin - energyStep*0.5, energyMax + energyStep*0.5 ) ;

  hEnergyvsTheta    -> Sumw2();

  hEnergyvsThetaGen -> Sumw2();

  hEnergyvsThetaKinCut      = new TH2D( "hEnergyvsThetaKinCut"     , ";#theta_{e'} (deg);E_{e'}",                                 thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax                                 ) ;

  hThetavsPhiGen            = new TH2D( "hThetavsPhiGen"           , ";#phi_{e'};#theta_{e'} (deg)",                                      360,   -180.0,    180.0, thetaPointN,   thetaMin, thetaMax                                  ) ;
 
  hThetavsPhi               = new TH2D( "hThetavsPhi"              , ";#phi_{e'};#theta_{e'} (deg)",                                      360,   -180.0,    180.0, thetaPointN,   thetaMin, thetaMax                                  ) ;

  hThetavsPhiFastMC         = new TH2D( "hThetavsPhiFastMC"        , ";#phi_{e'};#theta_{e'} (deg)",                                      360,   -180.0,    180.0, thetaPointN,   thetaMin, thetaMax                                  ) ;
 
  hQ2vsW                    = new TH2D( "hQ2vsW"                   , "Q^{2} (GeV^{2}); W (GeV)",                                  20*histosNBins,      0.0,      5.0,  20*histosNBins,        0,       15.0                                 ) ;
  hQ2vsWGen                 = new TH2D( "hQ2vsWGen"                , "Q^{2} (GeV^{2}); W (GeV)",                                  20*histosNBins,      0.0,      5.0,  20*histosNBins,        0,       15.0                                 ) ;

  hW                        = new TH1D( "hW"                       , "W (GeV)"                 ,                                  histosNBins,      0.8,      5.5                                                                     ) ;
  hWGen                     = new TH1D( "hWGen"                    , "W (GeV)"                 ,                                2*histosNBins,      0.8,      5.5                                                                     ) ;
  //hQ2Gen                    = new TH1D( "hQ2Gen"                   , "Q2(GeV2)"                ,                                10*histosNBins,     0.0,      1.5                                                                     ) ;

  hTheta                    = new TH1D( "hTheta"                   , "#theta_{e'} (deg)"       ,                                  thetaPointN, thetaMin,   thetaMax                                                                     ) ;
  hEnergy                   = new TH1D( "hEnergy"                  , "E_{e'} (GeV)"            ,                                  energyPointN, energyMin, 2*energyMax                                                                     ) ;
  hEnergyGen                = new TH1D( "hEnergyGen"               , "E_{e'} (GeV)"            ,                                  energyPointN, energyMin, 2*energyMax                                                                     ) ;
  hThetaGen                 = new TH1D( "hThetaGen"                , "#theta_{e'} (deg)"       ,                                  histosNBins, thetaMin,   thetaMax                                                                     ) ;

  /* hadron histos */

  /* Create the directory to store histos for hadrons */
  TDirectory* _dir_for_hadrons_plots = (TDirectory*)( output_file -> mkdir( "hadron_plots" ) );

  for( Int_t _ih = 0 ; _ih < n_hadron_types ; _ih++ ){

    hGenHadronEnergyvsTheta[ _ih ] = new TH2D( Form("hGenHadronEnergyvsTheta_h%d", _ih ), ";#theta_{e'} (deg);E_{e'}",                  thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax                                 ) ;

    hGenHadronThetavsPhi[    _ih ] = new TH2D( Form("hGenHadronThetacsPhi_h%d",    _ih ), ";#phi (deg);#theta_{e'} (deg)}",             thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax                                 ) ;

    hRecHadronEnergyvsTheta[ _ih ] = new TH2D( Form("hRecHadronEnergyvsTheta_h%d", _ih ), ";#theta_{e'} (deg);E_{e'}",                  thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax                                 ) ;

    hRecHadronThetavsPhi[    _ih ] = new TH2D( Form("hRecHadronThetacsPhi_h%d",    _ih ), ";#phi (deg);#theta_{e'} (deg)}",             thetaPointN, thetaMin, thetaMax, energyPointN, energyMin, energyMax                                 ) ;

    hGenHadronEnergyvsTheta[ _ih ] -> SetDirectory( _dir_for_hadrons_plots );
                                           
    hGenHadronThetavsPhi[    _ih ] -> SetDirectory( _dir_for_hadrons_plots );
                                           
    hRecHadronEnergyvsTheta[ _ih ] -> SetDirectory( _dir_for_hadrons_plots );
                                           
    hRecHadronThetavsPhi[    _ih ] -> SetDirectory( _dir_for_hadrons_plots );


  } /* it closes the loop on the hadorn types */

  gDirectory -> cd();

}

void InclusiveRateAnalysis::ClearVectors()
{

  thetaEle.clear(      ) ;
  energyEle.clear(     ) ;
  dsdThetadE.clear(    ) ;   
  dsdThetadErad.clear( ) ;
  W.clear(             ) ;            
  Q2.clear(            ) ;           
  dsdWdQ2.clear(       ) ;      
  dsdWdQ2rad.clear(    ) ;   

  dsdThetadEradVolume.clear();

  return ;

}

void InclusiveRateAnalysis::PlotCrossSectionDependences()
{

  TGraph* _ge_sigma_rad[ thetaPointN ];

  TMultiGraph *_mg_sigma_rad = new TMultiGraph();

  _mg_sigma_rad -> SetName( "mg_sigma_rad" );

  TGraph* _ge_rate[ thetaPointN ];

  TMultiGraph *_mg_rate = new TMultiGraph();

  _mg_rate -> SetName( "mg_rate" );

  TLegend* _leg = new TLegend( 0.87, 0.16, 0.97, 0.87 );

  _leg -> SetFillStyle( 1 );

  _leg -> SetFillColor( 0 );

  _leg -> SetLineColor( 0 );

  /* fill a tgraph for E dependencies at a fixed theta */

  vector<Int_t> vectorNPointE   ;

  vector<Int_t> vectorNPointTot ;

  vectorNPointE.clear() ;
                   
  vectorNPointTot.clear() ;

  for( Int_t _ith = 0; _ith < thetaPointN ; _ith++ ){

    Double_t _actualTheta = thetaMin + thetaStep*_ith ;

    Int_t _nPointE   = 0 ;

    Int_t _nPointTot = 0 ;

    Bool_t _thetaFound = kFALSE ;

    for( Int_t _ip = 0 ; _ip < fileNLines ; _ip++ ){

      if( TMath::Abs( thetaEle.at( _ip ) - _actualTheta ) < 0.00001 ){

	_nPointTot ++ ;

	_thetaFound = kTRUE ;

	if(  W.at( _ip ) > wCut && Q2.at( _ip ) > q2Cut  ) 	_nPointE++ ;

      } else if( _thetaFound ) {

	break ;

      }
    }

    // cout << Form( "There are n = %d energy points for theta = %g", _nPointE, _actualTheta ) << endl ;

    vectorNPointE.push_back(   _nPointE );

    vectorNPointTot.push_back( _nPointTot );

  }

  Int_t _n_check     = 0 ;
  Int_t _n_check_tot = 0 ;

  for( Int_t _i = 0; _i < (int)(vectorNPointE.size());   _i++ ) { _n_check     += vectorNPointE.at(   _i ); /* cout << Form( "theta = %g has %d points", thetaMin + thetaStep*_i, vectorNPointE.at(   _i ) ) << endl ; */ }
  for( Int_t _i = 0; _i < (int)(vectorNPointTot.size()); _i++ )   _n_check_tot += vectorNPointTot.at( _i );

  // cout << " tot lines = " << _n_check << "\t" << _n_check_tot << endl ;

  Int_t _istart = 0 ;

  Double_t _total_rate = 0.0 ;

  for( Int_t _ith = 0; _ith < thetaPointN ; _ith++ ){

    if( vectorNPointE.at( _ith ) ){

    Double_t _actualTheta = thetaMin + thetaStep*_ith ;

    _ge_sigma_rad[ _ith ] =  new TGraph( vectorNPointE.at( _ith ) );
    _ge_sigma_rad[ _ith ] -> SetName( Form( "ge_sigma_rad_theta_%g", _actualTheta ) ) ;

    _ge_rate[ _ith ]      =  new TGraph( vectorNPointE.at( _ith ) );
    _ge_rate[ _ith ]      -> SetName( Form( "ge_rate_theta_%g", _actualTheta ) ) ;

    gStyle -> SetPalette( 55 ) ; /* kRainBow */

    Int_t _color = 1 + _ith ;

    Int_t _i_multiple = (int)( _ith/99 ) ;

    if( _color > 99 ) _color -= 99*_i_multiple ;

    /* 19 is a tiny gray*/
    if(      _color == 10 ) _color = 38 ;
    else if( _color == 19 ) _color = 1  ;

    // cout << " Color set to " << _color << "\t" << _ith << "\t" << _i_multiple << endl ;

    _ge_sigma_rad[ _ith ] -> SetLineColor(   _color );

    _ge_sigma_rad[ _ith ] -> SetMarkerColor( _color );

    _ge_sigma_rad[ _ith ] -> SetLineWidth( 2 );

    _ge_rate[ _ith ]      -> SetLineColor(   _color );

    _ge_rate[ _ith ]      -> SetMarkerColor( _color );

    _ge_rate[ _ith ]      -> SetLineWidth( 1.3 );

    Int_t _ip = 0 ;

    for( Int_t _ie = _istart ; _ie < _istart + vectorNPointTot.at( _ith ) ; _ie++ ){

      if( W.at( _ie ) > wCut && Q2.at( _ie ) > q2Cut ){

	_ge_sigma_rad[ _ith ] -> SetPoint( _ip, energyEle.at( _ie ), dsdThetadErad.at( _ie ) );

	hEnergyvsThetaKinCut      -> Fill( _actualTheta, energyEle.at( _ie )                         ) ;

	/* calculate rate */
	/* it is calculated as dsigma/dEdO integrated over phi, multiplied for L=1035 cm-2 s-1 and dE=0.01GeV*/
	/* factor 10000 comes from cross-section in microbarn times 10^35 luminosity */
	Double_t _lumi_factor  = Luminosity*pow( 10, -30 ) ;

	Double_t _rate = dsdThetadErad.at( _ie )*_lumi_factor*TMath::Sin( _actualTheta*TMath::DegToRad() )*thetaStep*TMath::DegToRad()*energyStep ;

	Double_t _phi_acc_weight = AcceptanceWeight( energyEle.at( _ie ), _actualTheta ) ;

	_ge_rate[      _ith ] -> SetPoint( _ip, energyEle.at( _ie ), _rate*_phi_acc_weight*2.0*TMath::Pi() );

	/* total rate & integrated cross-section */
	_total_rate = _total_rate + _rate*_phi_acc_weight*2.0*TMath::Pi() ;

	/* the integrated cross-section is calculated as dsigma/dEdOmega x dE x dOmega (2*pi*sin(Theta) dTheta ) */
	integratedCrossSection    +=  dsdThetadE.at(       _ie )*energyStep*(2*TMath::Pi())*TMath::Sin( _actualTheta*TMath::DegToRad() )*thetaStep*TMath::DegToRad();

	integratedCrossSectionRad +=  dsdThetadErad.at(    _ie )*energyStep*(2*TMath::Pi())*TMath::Sin( _actualTheta*TMath::DegToRad() )*thetaStep*TMath::DegToRad();

	if( TMath::Abs( _actualTheta - thetaEle.at( _ie )) > 0.00001 ) cout << " WARNING! Check theta! " << _actualTheta << " \t " << thetaEle.at( _ie ) << endl ;

	_ip++ ;

      }

    } /* it closes the loop on the energies */

    _leg          -> AddEntry(  _ge_sigma_rad[ _ith ], Form( "#theta = %g deg", thetaMin + thetaStep*_ith ), "l" );

    _mg_sigma_rad -> Add(       _ge_sigma_rad[ _ith ], "L" );

    _mg_rate      -> Add(       _ge_rate[      _ith ], "L" );

    output_file -> WriteTObject( _ge_sigma_rad[ _ith ] );

    output_file -> WriteTObject( _ge_rate[      _ith ] );

    } /* it closes the condition of having a non-zero number of energy points */

    _istart = _istart + vectorNPointTot.at( _ith ) ;

  }  /* it closes the loop on the angles */

  cout << " Total rate is " << _total_rate/1000.0 << " kHz \n" << endl ;

  cout << " Integrated rad cross-section is " << integratedCrossSectionRad << " mubarn \n" << endl ;

  TCanvas* _c = new TCanvas("c_mg", "", 500, 1000 );

  // _c -> SetLogy();

  _c -> Divide( 1, 2 );

  _mg_sigma_rad -> SetTitle( Form( "#theta = %g - %g (deg), E_{e'} = %g - %g (GeV), E_{beam} = %g (GeV) ;E (GeV);d#sigma/dEd#Omega #mu b/GeV/str", thetaMin, thetaMax, energyMin, energyMax, energyBeam ) );

  _mg_rate      -> SetTitle( Form( "#theta = %g - %g (deg), E_{e'} = %g - %g (GeV), E_{beam} = %g (GeV) ;E (GeV);Inclusive Electron Rate/%g GeV/%g deg/sec", thetaMin, thetaMax, energyMin, energyMax, energyBeam, energyStep, thetaStep ) );

  // _mg_sigma_rad -> SetTitleFont(   132     );

  _c -> cd( 1 ) ;

  gPad -> SetLogy();

  _mg_sigma_rad -> Draw("A");

  _mg_sigma_rad -> GetXaxis() -> SetLabelFont(   132     );

  _mg_sigma_rad -> GetYaxis() -> SetLabelFont(   132     );

  _mg_sigma_rad -> GetXaxis() -> SetTitleFont(   132     );

  _mg_sigma_rad -> GetYaxis() -> SetTitleFont(   132     );

  _mg_sigma_rad -> SetMinimum( 0.0000001 ) ;

  _mg_sigma_rad -> SetMaximum( 100    ) ;

  _mg_sigma_rad -> Draw("A");

  TPaveText *_text_temp = new TPaveText( 5, 0.1, 6, 1 );

  _text_temp -> SetTextColor( kRed );

  _text_temp -> SetTextFont(  132 );

  _text_temp -> SetFillStyle(  0 );

  _text_temp -> SetBorderSize( 0 );

  _text_temp -> AddText( Form("%g (deg)", thetaMin ) );

  // _text_temp -> Draw("same");

  TPaveText *_text_temp_max = new TPaveText( 2.7, 0.000005 , 4.0, 0.0001 );

  _text_temp_max -> SetTextColor( kGreen + 2 );

  _text_temp_max -> SetTextFont(  132 );

  _text_temp_max -> SetFillStyle(  0 );

  _text_temp_max -> SetBorderSize( 0 );

  _text_temp_max -> AddText( Form("%g (deg)", thetaMax ) );

  // _text_temp_max -> Draw("same");

  _leg           -> Draw();

  _c -> cd( 2 ) ;

  gPad -> SetLogy();

  _mg_rate -> Draw("A");

  _mg_rate -> GetXaxis() -> SetLabelFont(   132     );

  _mg_rate -> GetYaxis() -> SetLabelFont(   132     );

  _mg_rate -> GetXaxis() -> SetTitleFont(   132     );

  _mg_rate -> GetYaxis() -> SetTitleFont(   132     );

  _mg_rate -> SetMinimum( 0.000001 ) ;

  _mg_rate -> SetMaximum( 10 ) ;

  _mg_rate -> Draw("A");


  output_file -> WriteTObject( _mg_sigma_rad );

  output_file -> WriteTObject( _mg_rate );

  output_file -> WriteTObject( _c );

  return ;

}

/* -------------------------- */
/*    Set kinematical limits  */
/* -------------------------- */

Double_t InclusiveRateAnalysis::AcceptanceWeight( Double_t _energy, Double_t _theta )
{

  Double_t _d_phi = -999.0 ;

  Double_t _theta_acc_1   =    5.0    ;

  Double_t _theta_acc_2   =   17.0    ;

  Double_t _theta_acc_3   =   35.0    ;

  Double_t _theta_acc_4   =      0.25 ;

  Double_t _theta_acc_5   =      0.2  ;

  Double_t _torus_max     = 3375.0    ;

  Double_t _torus_current = 3375.0    ; /* max torus current is 3375.0 */

  Double_t _p_norm        = _energy*_torus_max/TMath::Abs( _torus_current ) ;

  Double_t _theta_cut = _theta_acc_1 + _theta_acc_2/(( _energy + 0.05 )*_torus_max/TMath::Abs( _torus_current ) ) ;

  // cout << "_theta - theta_cut = " << _theta << "\t" <<  _theta_cut << endl ;

  if( _theta > _theta_cut){

    Double_t _exponent = _theta_acc_4*TMath::Power( _p_norm, _theta_acc_5 ) ;

    _d_phi = _theta_acc_3*TMath::Power( TMath::Sin( ( _theta - _theta_cut)*TMath::DegToRad() ), _exponent )/30.0 ;

    /* the /30 is due to the fact that d_phi is in the sector degrees */

  } else _d_phi = 0.0 ;

  // cout << " dphi is " << _d_phi << endl ;

  return _d_phi ;

}

/* -------------------- */
/* Prepare output file  */
/* -------------------- */
void InclusiveRateAnalysis::PrepareOutput( const char* _file_suffix )
{

  TString _out_name = file_name ;

  _out_name.ReplaceAll( ".dat", Form( "%s.root", _file_suffix ) );
  // _out_name.ReplaceAll( ".dat", "_torus1500A.root");

  output_file = new TFile( _out_name.Data(), "RECREATE");  gDirectory -> cd();

  return;

}

/* -------------------- */
/*   Save output file   */
/* -------------------- */
void InclusiveRateAnalysis::SaveOutput()
{

  output_file -> Write();

  output_file -> Close();

  return;

}
#endif
