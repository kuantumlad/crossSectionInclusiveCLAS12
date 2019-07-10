/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#ifndef _INCLUSIVERATEANALYSIS_H_
#define _INCLUSIVERATEANALYSIS_H_

#include<TH1.h>
#include<TH2.h>
#include<TH3.h>
#include<TGraph.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TLine.h>
#include<TStyle.h>
#include<TDirectory.h>
#include<TFile.h>
#include<TLorentzVector.h>
#include<TRandom.h>
#include<TRandom3.h>
#include<TPaveText.h>
#include<TLegend.h>
#include<TMultiGraph.h>
#include<map>
#include <iostream>
#include <fstream>
#include <stdio.h>     
#include <stdlib.h>    
#include <math.h>      
#include <TMath.h>     
#include <cctype>
#include <stdexcept>

#include <GenFunctions.h>
#include <DisFunctions.h>

using namespace std;


const int n_hadron_types = 16 ;

class InclusiveRateAnalysis {

 private:

  const char* file_name ;

  TFile* output_file ;

  map<int, int>    mapLundIdToNumber;

  map<int, double> mapLundIdToMass;

  /* steps and limits in the explored kinematics */

  Double_t thetaMin,  thetaMax,  thetaStep ;

  Double_t energyMin, energyMax, energyStep ;

  Double_t sigmaMin, sigmaMax ;

  Double_t energyBeam ;

  Int_t thetaPointN  ;

  Int_t energyPointN ;

  Int_t fileNLines ;

  Double_t wCut ;

  Double_t q2Cut ;

  Double_t integratedCrossSection    ;

  Double_t integratedCrossSectionRad ;

  Double_t integratedDISCrossSection ;

  /* Luminosity for rate estimate */
  Double_t Luminosity ;

  /* array with map of the hadron types to their lund id */
  Int_t  hadron_types[ n_hadron_types ];

  /* file format is the following */

  /* 35.00 3.300   0.00000360   0.00000299 13.13  1.48   0.00000049   0.00000041 */
  /* Theta  E'  Sigma/dE/dOmega   Sigma_radiated/dE/dOmega   Q2   W Sigma/dW/dQ2  Sigma_radiated/dW/dQ2 */

  vector<Double_t> thetaEle       ;
  vector<Double_t> energyEle      ;
  vector<Double_t> dsdThetadE     ;
  vector<Double_t> dsdThetadErad  ;
  vector<Double_t> W              ;
  vector<Double_t> Q2             ;
  vector<Double_t> dsdWdQ2        ;
  vector<Double_t> dsdWdQ2rad     ;

  /* cross-section times area */
  vector<Double_t> dsdThetadEradVolume ;

  /* vector to store java output */
  vector<Double_t> phiEleRec      ;
  vector<Double_t> thetaEleRec    ;
  vector<Double_t> energyEleRec   ;

  vector<Double_t> phiEleGen      ;
  vector<Double_t> thetaEleGen    ;
  vector<Double_t> energyEleGen   ;

  /* array to store the weights */
  Double_t _xs_weights[ 61 ][ 2101 ];

  /* vector< vector<Double_t> > dataValues ; */

  TH2D* hThetavsPhiGen        ;
  TH2D* hEnergyvsThetaGen     ;
  TH2D* hEnergyvsTheta        ;
  TH2D* hThetavsPhi           ;
  TH2D* hEnergyvsThetaKinCut  ;
  TH2D* hEnergyvsThetaFastMC  ;
  TH2D* hThetavsPhiFastMC     ;

  TH2D* hGenHadronEnergyvsTheta[ n_hadron_types ] ;
  TH2D* hGenHadronThetavsPhi[    n_hadron_types ] ;

  TH2D* hRecHadronEnergyvsTheta[ n_hadron_types ] ;
  TH2D* hRecHadronThetavsPhi[    n_hadron_types ] ;

  TH2D* hQ2vsW    ;
  TH2D* hQ2vsWGen ;

  TH1D* hEnergy    ;
  TH1D* hEnergyGen ;
  TH1D* hTheta     ;
  TH1D* hThetaGen  ;
  TH1D* hW         ;
  TH1D* hWGen      ;

  Int_t histosNBins    ;


 public:

  InclusiveRateAnalysis();

  virtual ~InclusiveRateAnalysis() {};

  void SetFile( const char* _file_name );

  void SetKinematicLimits( Double_t _theta_min, Double_t _theta_max, Double_t _theta_step, Double_t _energy_min, Double_t _energy_max, Double_t _energy_step, Int_t _file_n_lines, Double_t _e_beam = 11.0, Int_t _histos_n_bins = 200  );

  void ReadDataFile();

  void SetKinCuts( Double_t _q2_cut = 1.0, Double_t _w_cut = 2.0 ){ q2Cut = _q2_cut; wCut = _w_cut ; cout << Form( "Applying kinematics cuts Q^{2} > %g GeV^{2}, W > %g GeV", q2Cut, wCut ) << endl ; };

  void PlotCrossSectionDependences();

  void ClearVectors();

  void PrepareOutput( const char* _file_suffix = "" );

  void SaveOutput();

  void BookHistos();

  void GeneratePseudoData( Int_t _n_events = -1, Int_t _n_events_per_lund_file = 20000, const char* _dir_for_lund_files = "lund_files" );

  void TestGeneratedLund( Int_t _n_lund_files = 1, const char* _lund_file_body = "lund_files/out_lund_files_n" );

  /* void CreateTreeFromJavaDataFiles(); */

  void ReadJavaDataFiles( const char* _java_file_name, Int_t _n_java_files = 1 );

  void CompareGenEventsToOriginal();

  void CompareJavaEventsToGen( const char* _gen_root_file, const char* _java_root_file );

  void CompareFastMCToJava(    const char* _fastmc_root_file, const char* _java_root_file );

  void CompareNewFastMCToJava( const char* _fastmc_root_file, const char* _java_root_file, const char* _new_fastmc_root_file );

  void ReadJavaRecAndFastMC(   const char* _java_file_name, Int_t _n_java_files ) ;

  void ReadHadronJavaFiles( const char* _java_file_name, Int_t _n_java_files, Double_t _int_dis_cross_section );

  Double_t AcceptanceWeight( Double_t _energy, Double_t _theta );

  void SetLuminosity( Double_t _lumi_for_operations ){ Luminosity = _lumi_for_operations ; cout << Form( "\n Luminosity set to %g cm^-2 s^-1 \n", Luminosity ) << endl ; };

  ClassDef(InclusiveRateAnalysis,1)
    };

#endif
