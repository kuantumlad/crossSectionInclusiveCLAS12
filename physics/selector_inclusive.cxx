#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include<vector>
#include <algorithm>
#include <functional>
#include <TCutG.h>
using namespace std;

#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;


// similiar code to selector elastic but for inclusive analysis ( no W cut )


/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

//double Ebeam = 6.42313;
double Ebeam = 7.546;
//double Ebeam = 10.594;

int process_Events = -1;            // process all events
//int process_Events = 500000;     // process given number of events


/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// common variables
///

Float_t Pival = 3.14159265359;

Float_t m_e = 0.000511;
Float_t m_p = 0.93827;
Float_t m_n = 0.9396;
Float_t m_pip = 0.1396;
Float_t m_pim = 0.1396;
Float_t m_Kp = 0.4937;
Float_t m_Km = 0.4937;
Float_t c = 299792458;

// vectors with components of the Lorentzvector to fill into the tree
vector<int> *p4_mc_pid = 0;
vector<double> *p4_mc_px  = 0;
vector<double> *p4_mc_py = 0;
vector<double> *p4_mc_pz = 0;
vector<double> *p4_mc_vx = 0;
vector<double> *p4_mc_vy = 0;
vector<double> *p4_mc_vz = 0;


int helicity;
float fcup;
vector<int > *sectorE = 0;
vector<double> *p4_ele_px = 0;
vector<double> *p4_ele_py = 0;
vector<double> *p4_ele_pz = 0;
vector<double> *p4_ele_E = 0;
vector<double> *p4_prot_px = 0;
vector<double> *p4_prot_py = 0;
vector<double> *p4_prot_pz = 0;
vector<double> *p4_prot_E = 0;
vector<double> *p4_neutr_px = 0;
vector<double> *p4_neutr_py = 0;
vector<double> *p4_neutr_pz = 0;
vector<double> *p4_neutr_E = 0;
vector<double> *p4_pip_px = 0;
vector<double> *p4_pip_py = 0;
vector<double> *p4_pip_pz = 0;
vector<double> *p4_pip_E = 0;
vector<double> *p4_pim_px = 0;
vector<double> *p4_pim_py = 0;
vector<double> *p4_pim_pz = 0;
vector<double> *p4_pim_E = 0;
vector<double> *p4_Kp_px = 0;
vector<double> *p4_Kp_py = 0;
vector<double> *p4_Kp_pz = 0;
vector<double> *p4_Kp_E = 0;
vector<double> *p4_Km_px = 0;
vector<double> *p4_Km_py = 0;
vector<double> *p4_Km_pz = 0;
vector<double> *p4_Km_E = 0;
vector<double> *p4_phot_px = 0;
vector<double> *p4_phot_py = 0;
vector<double> *p4_phot_pz = 0;
vector<double> *p4_phot_E = 0;
vector<double> *prot_beta_final = 0;
vector<double> *Kp_beta_final = 0;
vector<double> *Km_beta_final = 0;
vector<double> *el_sector = 0;
vector<int> *eventNumber = 0;
vector<double> *prot_sector = 0;
vector<double> *Kp_sector = 0;
vector<double> *Km_sector = 0;


/// varibales for particles:

const static int BUFFER = 20;

int ele_sector[BUFFER];
TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER]; 
TLorentzVector neutr[BUFFER]; 
TLorentzVector pip[BUFFER]; 
TLorentzVector pim[BUFFER]; 
TLorentzVector Kp[BUFFER]; 
TLorentzVector Km[BUFFER]; 
TLorentzVector phot[BUFFER];
TLorentzVector mc_ele[BUFFER];
TLorentzVector mc_prot[BUFFER];

//added for beta vs p
double prot_beta[BUFFER];
double kaonP_beta[BUFFER];
double kaonM_beta[BUFFER];

int sect_el[BUFFER];
int sect_pr[BUFFER];
int sect_kp[BUFFER];
int sect_km[BUFFER];
int event[BUFFER];

///  counting variables:

double ele_count;
double prot_count;
double neutr_count;
double pip_count;
double pim_count;
double Kp_count;
double Km_count;
double phot_count;

double mc_count;

// kinematic variables

double W, Q2, nu, x, y;
double t1, cmphi1, cmcostheta1, pt1, eta1, z1;
double M_e_p_X_miss, M_e_p_X_miss2;

double W_out, Q2_out, nu_out, x_out, y_out;
double t1_out, cmphi1_out, cmcostheta1_out, pt1_out, eta1_out, z1_out;
double perp_mntm, missing_e, missing_mm2, epX_mm2, epX_mm;
 

double W_2, Q2_2, nu_2, x_2, y_2;
double t1_2, cmphi1_2, cmcostheta1_2, pt1_2, eta1_2, z1_2;


///////////////////////////////////////////////////////////////////////////////////////////////////////
///  selected particles:

Int_t evcat; 
Double_t E_ele; Double_t px_ele; Double_t py_ele; Double_t pz_ele;
Double_t E_prot_1; Double_t px_prot_1; Double_t py_prot_1; Double_t pz_prot_1;

//ADDED 
Double_t E_kaonP_1; Double_t px_kaonP_1; Double_t py_kaonP_1; Double_t pz_kaonP_1;
Double_t E_kaonM_1; Double_t px_kaonM_1; Double_t py_kaonM_1; Double_t pz_kaonM_1;

Double_t prot_beta_out; Double_t kaonP_beta_out; Double_t kaonM_beta_out;
Int_t el_sect; Int_t pr_sect; Int_t kp_sect; Int_t km_sect;


// Monte Carlo Particles for elastic events 
Double_t E_ele_mc; Double_t px_ele_mc; Double_t py_ele_mc; Double_t pz_ele_mc;
Double_t E_prot_1_mc; Double_t px_prot_1_mc; Double_t py_prot_1_mc; Double_t pz_prot_1_mc;


/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  index of selected particle:

int select_ele;
int select_prot_1; 

//NEW
int select_kaonP_1;
int select_kaonM_1;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  input and output files:

TFile *out;

char name[200];
char title[200];

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// define functions

double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double kin_x(TLorentzVector ele);
double kin_y(TLorentzVector ele);
double kin_nu(TLorentzVector ele);

double kin_pT(TLorentzVector ele, TLorentzVector hadron);
double kin_eta(TLorentzVector ele, TLorentzVector hadron);
double kin_z(TLorentzVector ele, TLorentzVector hadron);
double kin_t(TLorentzVector ele, TLorentzVector hadron);

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron);
double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron);

double kin_mismass(TLorentzVector ele, TLorentzVector hadron);
double kin_mismass2(TLorentzVector ele, TLorentzVector hadron);

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron);

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron);
double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2);

double Phi_Mass(TLorentzVector kaonP, TLorentzVector kaonM);

//ADDED 
double kin_epkXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp);
double kin_epKpKmXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp, TLorentzVector km);



/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// main
///

Int_t selector_inclusive( Char_t *inFile, Char_t *outputfile, int run, std::string data_type = "DATA")
{
		
  Char_t *inTree="out_tree";
  cout << "Initalize the input tree ... " << endl;

  /// ///////////////////////////////////////////////////////////////////////////
  ///  get input file:
  /// ///////////////////////////////////////////////////////////////////////////

  TFile *f; 
  Char_t tmpstr[80];
  Double_t fraction;

  f = new TFile(inFile,"");   // Input File

  if(f->IsZombie()){   // Check if TFile exists!
    cout<<"Input file " << inFile << "doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFile << endl;


  Char_t *mcTreeName="mc_tree";
  TTree *mcTree=(TTree *) f->Get(mcTreeName);
  std::string analysis_sim = "SIM";
  if(mcTree==0){  // Check if TTree exists!
    cout << " MC Tree doesn't exist!!!" << endl;
    cout <<"Exit program" << endl;
    //return 0;
  }

  if( data_type == analysis_sim ){    
  std::cout << " SETTING MC VARIABLES " << std::endl;
  mcTree->SetBranchAddress("p4_mc_pid",&p4_mc_pid);
  mcTree->SetBranchAddress("p4_mc_px",&p4_mc_px);
  mcTree->SetBranchAddress("p4_mc_py",&p4_mc_py);
  mcTree->SetBranchAddress("p4_mc_pz",&p4_mc_pz);
  mcTree->SetBranchAddress("p4_mc_vx",&p4_mc_vx);
  mcTree->SetBranchAddress("p4_mc_vy",&p4_mc_vy);
  mcTree->SetBranchAddress("p4_mc_vz",&p4_mc_vz);
  }


  TTree *anaTree=(TTree *) f->Get(inTree);

  if(anaTree==0){  // Check if TTree exists!
    cout << "Tree " << inTree << " doesn't exist!!!" << endl;
    cout <<"Exit program" << endl;
    return 0;
  }
 
    

    
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  get branches from input file:
  /// ///////////////////////////////////////////////////////////////////////////// 
    
  anaTree->SetBranchAddress("helicity", &helicity);
  anaTree->SetBranchAddress("fcup", &fcup);
  anaTree->SetBranchAddress("sectorE",&sectorE);
  anaTree->SetBranchAddress("eventNumber",&eventNumber);
  anaTree->SetBranchAddress("p4_ele_px", &p4_ele_px);
  anaTree->SetBranchAddress("p4_ele_py", &p4_ele_py);
  anaTree->SetBranchAddress("p4_ele_pz", &p4_ele_pz);
  anaTree->SetBranchAddress("p4_ele_E", &p4_ele_E);
  anaTree->SetBranchAddress("p4_prot_px", &p4_prot_px);
  anaTree->SetBranchAddress("p4_prot_py", &p4_prot_py);
  anaTree->SetBranchAddress("p4_prot_pz", &p4_prot_pz);
  anaTree->SetBranchAddress("p4_prot_E", &p4_prot_E);
  anaTree->SetBranchAddress("p4_neutr_px", &p4_neutr_px);
  anaTree->SetBranchAddress("p4_neutr_py", &p4_neutr_py);
  anaTree->SetBranchAddress("p4_neutr_pz", &p4_neutr_pz);
  anaTree->SetBranchAddress("p4_neutr_E", &p4_neutr_E);
  anaTree->SetBranchAddress("p4_pip_px", &p4_pip_px);
  anaTree->SetBranchAddress("p4_pip_py", &p4_pip_py);
  anaTree->SetBranchAddress("p4_pip_pz", &p4_pip_pz);
  anaTree->SetBranchAddress("p4_pip_E", &p4_pip_E);
  anaTree->SetBranchAddress("p4_pim_px", &p4_pim_px);
  anaTree->SetBranchAddress("p4_pim_py", &p4_pim_py);
  anaTree->SetBranchAddress("p4_pim_pz", &p4_pim_pz);
  anaTree->SetBranchAddress("p4_pim_E", &p4_pim_E);
  anaTree->SetBranchAddress("p4_Kp_px", &p4_Kp_px);
  anaTree->SetBranchAddress("p4_Kp_py", &p4_Kp_py);
  anaTree->SetBranchAddress("p4_Kp_pz", &p4_Kp_pz);
  anaTree->SetBranchAddress("p4_Kp_E", &p4_Kp_E);
  anaTree->SetBranchAddress("p4_Km_px", &p4_Km_px);
  anaTree->SetBranchAddress("p4_Km_py", &p4_Km_py);
  anaTree->SetBranchAddress("p4_Km_pz", &p4_Km_pz);
  anaTree->SetBranchAddress("p4_Km_E", &p4_Km_E);
  anaTree->SetBranchAddress("p4_phot_px", &p4_phot_px);
  anaTree->SetBranchAddress("p4_phot_py", &p4_phot_py);
  anaTree->SetBranchAddress("p4_phot_pz", &p4_phot_pz);
  anaTree->SetBranchAddress("p4_phot_E", &p4_phot_E);  

  // }
  
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  SET OUTPUT FILE 
  /// ///////////////////////////////////////////////////////////////////////////// 
  out = new TFile(outputfile, "RECREATE");
  
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  SET histograms
  /// ///////////////////////////////////////////////////////////////////////////// 

  double el_theta_max = 30.0;
  double pr_theta_max = 50.0;

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;

  TH1F *hist_all_el_energy;

  hist_all_electron_p = new TH1F("hist_all_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_electron_p->GetYaxis()->SetTitle("counts");
  hist_all_electron_theta = new TH1F("hist_all_electron_theta", "electron #Theta", 140,0,el_theta_max);   
  hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_theta->GetYaxis()->SetTitle("counts");
  hist_all_electron_phi = new TH1F("hist_all_electron_phi", "electron #phi", 73,-180,180);   
  hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_phi->GetYaxis()->SetTitle("counts");
  hist_all_electron_p_vs_theta = new TH2F("hist_all_electron_p_vs_theta", "electron p vs #Theta", 140,0,el_theta_max,500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_theta_vs_phi = new TH2F("hist_all_electron_theta_vs_phi", "electron #Theta vs #phi", 73, -180,180, 30 , 0, el_theta_max);   
  hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");
  hist_all_el_energy = new TH1F("hist_all_el_energy","hist_all_el_energy",200, -0.5, 0.8 );
  hist_all_el_energy->GetXaxis()->SetTitle("#delta E [GeV]");

  hist_all_proton_p = new TH1F("hist_all_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_proton_p->GetYaxis()->SetTitle("counts");
  hist_all_proton_theta = new TH1F("hist_all_proton_theta", "proton #Theta", 140,0,el_theta_max);   
  hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_theta->GetYaxis()->SetTitle("counts");
  hist_all_proton_phi = new TH1F("hist_all_proton_phi", "proton #phi", 180,-180,180);   
  hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_phi->GetYaxis()->SetTitle("counts");
  hist_all_proton_p_vs_theta = new TH2F("hist_all_proton_p_vs_theta", "proton p vs #Theta", 140,0,pr_theta_max,500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_p_vs_phi = new TH2F("hist_all_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_theta_vs_phi = new TH2F("hist_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 140,0,pr_theta_max);   
  hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  out->mkdir("particle_histograms_selected");				
  out->cd ("particle_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi; 
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;
  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi; 
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;

  hist_electron_p = new TH1F("hist_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_electron_p->GetYaxis()->SetTitle("counts");
  hist_electron_theta = new TH1F("hist_electron_theta", "electron #Theta", 140,0,40);   
  hist_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_theta->GetYaxis()->SetTitle("counts");
  hist_electron_phi = new TH1F("hist_electron_phi", "electron #phi", 73,-180,180);   
  hist_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_phi->GetYaxis()->SetTitle("counts");
  hist_electron_p_vs_theta = new TH2F("hist_electron_p_vs_theta", "electron p vs #Theta", 140,0,40,500,0,Ebeam+0.3);   
  hist_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_electron_p_vs_phi = new TH2F("hist_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_electron_theta_vs_phi = new TH2F("hist_electron_theta_vs_phi", "electron #Theta vs phi", 73,-180,180, 30,0,el_theta_max);   
  hist_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_proton_p = new TH1F("hist_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_proton_p->GetYaxis()->SetTitle("counts");
  hist_proton_theta = new TH1F("hist_proton_theta", "proton #Theta", 140,0, el_theta_max);   
  hist_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_theta->GetYaxis()->SetTitle("counts");
  hist_proton_phi = new TH1F("hist_proton_phi", "proton #phi", 180,-180,180);   
  hist_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_phi->GetYaxis()->SetTitle("counts");
  hist_proton_p_vs_theta = new TH2F("hist_proton_p_vs_theta", "proton p vs #Theta", 140,0,pr_theta_max,500,0,Ebeam+0.3);   
  hist_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_proton_p_vs_phi = new TH2F("hist_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_proton_theta_vs_phi = new TH2F("hist_proton_theta_vs_phi", "proton #Theta vs phi", 180,-180,180, 140,0,pr_theta_max);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");
  
  out->mkdir("kinematics");				
  out->cd ("kinematics");

  // Look at total phase space here as supplied to the detector
  TH1F *hist_W;
  TH1F *hist_Q2;
  TH1F *hist_x;
  TH1F *hist_y;
  TH1F *hist_nu;
  TH1F *hist_t;
  TH1F *hist_cmphi;
  TH1F *hist_cmcostheta;

  //ADDED 
  TH2F *hist_Q2_x;
  TH2F *hist_x_t;
  TH2F *hist_t_phi;
  TH2F *hist_Q2_phi;
  TH2F *hist_Q2_vs_W;
  TH2F *hist_W_vs_phi;

  std::vector< TH1F* > h_el_w_sect;
  std::vector< TH1F* > h_el_theta_sect;
  std::vector< TH2F* > h_el_ptheta_sect;

  std::vector< TH1F* > h_el_p_sect_final;
  std::vector< TH1F* > h_el_theta_sect_final;
  std::vector< TH1F* > h_el_phi_sect_final;

  std::vector< TH2F* > h_el_ptheta_sect_final;
  std::vector< TH2F* > h_el_phitheta_sect_final;

  TH2F *h_el_phitheta_final = new TH2F(Form("h_el_phitheta_final_run%d",run), Form("h_el_phitheta_final_run%d",run), 73, -180.0, 180.0, 30, 0.0, el_theta_max );

  for( int s = 0; s < 7; s++ ){
    h_el_w_sect.push_back( new TH1F(Form("h_el_w_s%d",s),Form("h_el_w_s%d",s), 100, 0.9, 1.2) );
    h_el_theta_sect.push_back( new TH1F(Form("h_el_theta_s%d",s),Form("h_el_theta_s%d",s), 100, 0.0, 40.0) );
   
    h_el_ptheta_sect.push_back( new TH2F(Form("h_el_ptheta_s%d",s),Form("h_el_ptheta_s%d",s), 200, 0.0, 2.5, 200, 0.0, 40.0) );
    
    h_el_p_sect_final.push_back( new TH1F(Form("h_el_p_s%d_final",s),Form("h_el_p_s%d_final",s),100, 0.0, 2.5) );
    h_el_theta_sect_final.push_back( new TH1F(Form("h_el_theta_s%d_final",s),Form("h_el_theta_s%d_final",s),30, 0.0, el_theta_max) );
    h_el_phi_sect_final.push_back( new TH1F(Form("h_el_phi_s%d_final",s),Form("h_el_phi_s%d_final",s), 73, -180.0, 180.0 ) );

    h_el_ptheta_sect_final.push_back( new TH2F(Form("h_el_ptheta_s%d_final",s),Form("h_el_ptheta_s%d_final",s), 200, 0.0, 2.5, 200, 0.0, el_theta_max) );
    h_el_phitheta_sect_final.push_back( new TH2F(Form("h_el_phitheta_s%d_final",s),Form("h_el_phitheta_s%d_final",s), 73, -180.0, 180.0, 30, 0.0, el_theta_max) );
  }


  out->mkdir("kinematics_inclusive");				
  out->cd ("kinematics_inclusive");

  std::vector< TH2F* > h_el_q2w_sect;
  std::vector<TH1D*> h_el_w_q2bin;

  for( int s = 0; s < 7; s++ ){      
    h_el_q2w_sect.push_back( new TH2F(Form("h_el_q2w_s%d",s),Form("h_el_q2w_s%d",s), 200, 1.1, 2.1, 200, 0.0, 4.0) );
  }

  /// /////////////////////////////////////////////////////////////////////////////    
  ///  create output tree:
  /// /////////////////////////////////////////////////////////////////////////////


  // categroy 1:  e p -->  e p X

  TTree out_tree1("events_eX","events_eX");
  
  out_tree1.Branch("W", &W);
  out_tree1.Branch("Q2", &Q2);
  out_tree1.Branch("x", &x);
  out_tree1.Branch("y", &y);
  out_tree1.Branch("nu", &nu);
  out_tree1.Branch("minus_t", &t1);
  out_tree1.Branch("cmphi", &cmphi1);
  out_tree1.Branch("cmcostheta", &cmcostheta1);
  out_tree1.Branch("pt", &pt1);
  out_tree1.Branch("eta", &eta1);
  out_tree1.Branch("z", &z1);
  out_tree1.Branch("missmass", &M_e_p_X_miss);
  out_tree1.Branch("missmass2", &M_e_p_X_miss2);
  out_tree1.Branch("helicity", &helicity);
  out_tree1.Branch("E_ele", &E_ele);
  //out_tree1.Branch("eventNumber",&eventNumber);
  out_tree1.Branch("px_ele", &px_ele);
  out_tree1.Branch("py_ele", &py_ele);
  out_tree1.Branch("pz_ele", &pz_ele);
  out_tree1.Branch("E_prot", &E_prot_1);
  out_tree1.Branch("px_prot", &px_prot_1);
  out_tree1.Branch("py_prot", &py_prot_1);
  out_tree1.Branch("pz_prot", &pz_prot_1);
 
  out_tree1.Branch("E_ele_mc", &E_ele_mc);
  out_tree1.Branch("px_ele_mc", &px_ele_mc);
  out_tree1.Branch("py_ele_mc", &py_ele_mc);
  out_tree1.Branch("pz_ele_mc", &pz_ele_mc);
  out_tree1.Branch("E_prot_mc", &E_prot_1_mc);
  out_tree1.Branch("px_prot_mc", &px_prot_1_mc);
  out_tree1.Branch("py_prot_mc", &py_prot_1_mc);
  out_tree1.Branch("pz_prot_mc", &pz_prot_1_mc);
  
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  Set cut limits on W and other variables
  /// /////////////////////////////////////////////////////////////////////////////
    

  
  cout << "Analysing Tree: " << inTree << endl;
  cout << "Event Loop starting ... " << endl;
  cout << endl;
        
  //std::cout << " Number of MC events " << mcTree->GetEntriesFast() << std::endl;
  for(Int_t k=0; k < anaTree->GetEntriesFast(); k++){    
      
    anaTree->GetEntry(k);

    if(process_Events > 0 && k == process_Events) break;

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// progress:

    if(k % 100000 == 0){

      double events = anaTree->GetEntriesFast();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
    }

    //pr_sect=sect_pr[select_prot_1];//
    //kp_sect=sect_kp[select_kaonP_1];//
    //km_sect=sect_kp[select_kaonM_1];

    //std::cout << " pr sector " << pr_sect << std::endl;
	//	std::cout << " kp sector " << kp_sect << std::endl;
    //	std::cout << " km sector " << km_sect << std::endl;


    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// assign number of particles of each type PER EVENT:
    // (SIZE OF ARRAYS)
    ele_count = p4_ele_px->size();
    prot_count = p4_prot_px->size();
    neutr_count = p4_neutr_px->size();
    pip_count = p4_pip_px->size();
    pim_count = p4_pim_px->size();
    Kp_count = p4_Kp_px->size();
    Km_count = p4_Km_px->size();
    phot_count = p4_phot_px->size();

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// initalisatize variables:

    TLorentzVector beam(0,0,Ebeam,Ebeam);
    TLorentzVector target(0,0,0,0.93827);

    ///  index of selected particle for the different categories:
    // out tree variables
    select_ele = 0;  
    select_prot_1 = 0;
    evcat = 0; 
    E_ele = 0; px_ele = 0; py_ele = 0; pz_ele = 0;
    E_prot_1 = 0; px_prot_1 = 0; py_prot_1 = 0; pz_prot_1 = 0;
    //ADDED
    E_kaonP_1 = 0; px_kaonP_1 = 0; py_kaonP_1 = 0; pz_kaonP_1 = 0;
    E_kaonM_1 = 0; px_kaonM_1 = 0; py_kaonM_1 = 0; pz_kaonM_1 = 0;

    //Moncte carlo variables (GEN)
    E_ele_mc = 0; px_ele_mc = 0; py_ele_mc = 0; pz_ele_mc = 0;
    E_prot_1_mc = 0; px_prot_1_mc = 0; py_prot_1_mc = 0; pz_prot_1_mc = 0;

    W = 0; Q2 = 0; nu = 0; x = 0; y = 0;
    t1 = 0; cmphi1 = 0; cmcostheta1 = 0; pt1 = 0; eta1 = 0; z1 = 0;

    M_e_p_X_miss = 0; M_e_p_X_miss2 = 0;

    el_sect=0; pr_sect=0; kp_sect=0; km_sect=0;
    
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Assign tree components to Lorentz vectors:
    /// //////////////////////////////////////////////////////////////////

    for(Int_t i = 0; i < BUFFER; i++){ 
      if(ele_count > i){ele[i].SetPxPyPzE(p4_ele_px->at(i),p4_ele_py->at(i),p4_ele_pz->at(i),p4_ele_E->at(i)); ele_sector[i]=sectorE->at(i); }//event[i]=eventNumber->at(i);}
      else{ele[i].SetPxPyPzE(0, 0, 0, 0); ele_sector[i] = -1;}
      if(prot_count > i){prot[i].SetPxPyPzE(p4_prot_px->at(i),p4_prot_py->at(i),p4_prot_pz->at(i),p4_prot_E->at(i));}
      else{prot[i].SetPxPyPzE(0, 0, 0, 0);}
      if(neutr_count > i){neutr[i].SetPxPyPzE(p4_neutr_px->at(i),p4_neutr_py->at(i),p4_neutr_pz->at(i),p4_neutr_E->at(i));}
      else{neutr[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pip_count > i){pip[i].SetPxPyPzE(p4_pip_px->at(i),p4_pip_py->at(i),p4_pip_pz->at(i),p4_pip_E->at(i));}
      else{pip[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pim_count > i){pim[i].SetPxPyPzE(p4_pim_px->at(i),p4_pim_py->at(i),p4_pim_pz->at(i),p4_pim_E->at(i));}
      else{pim[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Kp_count > i){Kp[i].SetPxPyPzE(p4_Kp_px->at(i),p4_Kp_py->at(i),p4_Kp_pz->at(i),p4_Kp_E->at(i));}
      else{Kp[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Km_count > i){Km[i].SetPxPyPzE(p4_Km_px->at(i),p4_Km_py->at(i),p4_Km_pz->at(i),p4_Km_E->at(i));}
      else{Km[i].SetPxPyPzE(0, 0, 0, 0);}
      if(phot_count > i){phot[i].SetPxPyPzE(p4_phot_px->at(i),p4_phot_py->at(i),p4_phot_pz->at(i),p4_phot_E->at(i));}
      else{phot[i].SetPxPyPzE(0, 0, 0, 0);}      
    }

    double mass_el = 0.000511;

    /*  
    mcTree->GetEntry(event[0]);
    	  
    for( int i = 0; i < p4_mc_pid->size(); i++ ){
      if( p4_mc_pid->at(i) == 11 ){
	int ii = i;
	double mc_ele_e = sqrt( p4_mc_px->at(ii)*p4_mc_px->at(ii) + 
				p4_mc_py->at(ii)*p4_mc_py->at(ii) +
				p4_mc_pz->at(ii)*p4_mc_pz->at(ii) + 
				mass_el * mass_el );
	mc_ele[ii].SetPxPyPzE(p4_mc_px->at(ii), p4_mc_py->at(ii), p4_mc_pz->at(ii), mc_ele_e );
      }
    }
    */				 
    
    // Assign beta from tree to array
    for( Int_t i = 0; i < BUFFER; i++ ){
      //if( prot_count > i ){ prot_beta[i]=prot_beta_final->at(i); }
      //if( Kp_count > i ){ kaonP_beta[i]=Kp_beta_final->at(i); }
      //if( Km_count > i ){ kaonM_beta[i]=Km_beta_final->at(i); }

      // and assign DC sector for hadrons from tree to array
      //if( prot_count > i ){ sect_pr[i]=prot_sector->at(i); }
      //if( Kp_count > i ){ sect_kp[i]=Kp_sector->at(i); }
      //if( Km_count > i ){ sect_km[i]=Km_sector->at(i); }

    }




    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Build event:
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(ele_count == 1){    // one detected electron is the basic trigger condition for all topologies 

      select_ele = 0;   // if more than one electron is detected (< 0.2 percent of the events), take the one with the hihest momentum

      W  = kin_W(ele[select_ele]);
      Q2 = kin_Q2(ele[select_ele]);
      x  = kin_x(ele[select_ele]);
      y  = kin_y(ele[select_ele]);
      nu = kin_nu(ele[select_ele]);
      
      E_ele  = ele[select_ele].E();
      px_ele = ele[select_ele].Px();
      py_ele = ele[select_ele].Py();
      pz_ele = ele[select_ele].Pz();

      // Fill all electrons

      for(Int_t i = 0; i < BUFFER; i++){ 
	if(ele[i].P() > 0) hist_all_electron_p->Fill(ele[i].P());
	if(ele[i].Theta() > 0) hist_all_electron_theta->Fill(ele[i].Theta()*180/Pival);
	if(ele[i].Phi() != 0) hist_all_electron_phi->Fill(ele[i].Phi()*180/Pival);
	if(ele[i].P() > 0) hist_all_electron_p_vs_theta->Fill(ele[i].Theta()*180/Pival, ele[i].P());
	if(ele[i].P() > 0) hist_all_electron_p_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].P());
	if(ele[i].Theta() > 0 && ele[i].Phi() != 0) hist_all_electron_theta_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].Theta()*180/Pival);
      }


      int el_sect = ele_sector[select_ele];
      if( el_sect >= 1 ){
	h_el_ptheta_sect[el_sect]->Fill( ele[select_ele].P(), ele[select_ele].Theta()*180/Pival );
	//	h_el_w_sect[el_sect]->Fill( W );
	h_el_theta_sect[el_sect]->Fill(ele[select_ele].Theta()*180/Pival);
      	
	double calculated_energy = Ebeam/( 1 + (Ebeam/0.938)*(1 - cos(ele[select_ele].Theta())) );
	double measured_energy = ele[select_ele].E();
	double delta_energy = calculated_energy - measured_energy;
       
	if( true ){ //W > w_cut_min && W < w_cut_max ){
	  if( (ele[select_ele].Theta()*180/Pival) > 5.0 ){

	    // Fill best electron
	    W  = kin_W(ele[select_ele]);
	    Q2 = kin_Q2(ele[select_ele]);
	    x  = kin_x(ele[select_ele]);
	    y  = kin_y(ele[select_ele]);
	    nu = kin_nu(ele[select_ele]);

	    h_el_q2w_sect[el_sect]->Fill( W, Q2 );


	    E_ele  = ele[select_ele].E();
	    px_ele = ele[select_ele].Px();
	    py_ele = ele[select_ele].Py();
	    pz_ele = ele[select_ele].Pz();

	    E_ele_mc  = mc_ele[0].E();
	    px_ele_mc = mc_ele[0].Px();
	    py_ele_mc = mc_ele[0].Py();
	    pz_ele_mc = mc_ele[0].Pz();
	    out_tree1.Fill();
	   
	    if(ele[select_ele].P() > 0) hist_electron_p->Fill(ele[select_ele].P());
	    if(ele[select_ele].Theta() > 0) hist_electron_theta->Fill(ele[select_ele].Theta()*180/Pival);
	    if(ele[select_ele].Phi() != 0) hist_electron_phi->Fill(ele[select_ele].Phi()*180/Pival);
	    if(ele[select_ele].P() > 0) hist_electron_p_vs_theta->Fill(ele[select_ele].Theta()*180/Pival, ele[select_ele].P());
	    if(ele[select_ele].P() > 0) hist_electron_p_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].P());
	    if(ele[select_ele].Theta() > 0 && ele[select_ele].Phi() != 0) hist_electron_theta_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].Theta()*180/Pival);


	    hist_all_el_energy->Fill(delta_energy);
	    h_el_w_sect[el_sect]->Fill( W );
	    h_el_p_sect_final[el_sect]->Fill(ele[select_ele].P()); 
	    h_el_theta_sect_final[el_sect]->Fill(ele[select_ele].Theta()*180/Pival); 
	    h_el_phi_sect_final[el_sect]->Fill(ele[select_ele].Phi()*180/Pival);
	    h_el_ptheta_sect_final[el_sect]->Fill( ele[select_ele].P(), ele[select_ele].Theta()*180/Pival );   
	    h_el_phitheta_sect_final[el_sect]->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].Theta()*180/Pival);
	    h_el_phitheta_final->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].Theta()*180/Pival); 
	    //

	  }

	  //}
  
	  

	}
      }

 
    }

  }
  

  cout << endl;
  cout << "Tree successfully analysed!" << endl;
  cout << "Writing the output file ... " << endl;
  out->Write(); // Saving Histograms
  out_tree1.Write();
  cout << "Histograms saved in File: " << outputfile << endl;
  out->Close(); // Closing Output File
  f->Close();   // Closing Input File
  cout << "... Completed!" << endl;


  
  return 1;

}

///////////////////////////////////////////////////////////////
/// angles between particles

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;

  return ele.Vect().Angle(miss.Vect());
}

double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2){

  return particle1.Vect().Angle(particle2.Vect());
}


///////////////////////////////////////////////////////////////
///  kinematics

double kin_W(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}

double kin_Q2(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2();
}

double kin_x(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2()/(2*0.938272*fGamma.E());
}


double kin_y(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E()/Ebeam;
}

double kin_nu(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E();
}

double kin_pT(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Perp(); 
}

double kin_eta(TLorentzVector ele, TLorentzVector hadron){ 
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Rapidity(); 
}

double kin_z(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;  
  return hadron.T()/fGamma.T();
}

double kin_t(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector transfer = fGamma - hadron;
  return transfer.M2();
}

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();  // get boost vector to CM frame

  TLorentzVector gammaClone(fGamma);    // clone particles
  TLorentzVector eleClone(ele); 
  TLorentzVector hadronClone(hadron); 

  gammaClone.Boost(boost);     // boost particles to CM frame
  eleClone.Boost(boost);
  hadronClone.Boost(boost);

  TVector3 v0, v1;
  v0 = gammaClone.Vect().Cross(eleClone.Vect());
  v1 = gammaClone.Vect().Cross(hadronClone.Vect());    // hadronClone.Vect().Cross(gammaClone.Vect());

  Double_t c0, c1, c2, c3;
  c0 = v0.Dot(hadronClone.Vect());
  c1 = v0.Dot(v1);
  c2 = v0.Mag();
  c3 = v1.Mag();

  return (c0/TMath::Abs(c0)) * TMath::ACos(c1 /(c2*c3));
}

double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();

  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);

  return TMath::Cos(hadronClone.Theta()); 
}


/////////////////////////////////////////////////////////////////////
///  missing mass and energy

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return miss;
}

double kin_mismass(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return sqrt(miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

double kin_mismass2(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

