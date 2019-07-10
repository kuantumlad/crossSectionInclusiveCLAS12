/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#ifndef _GenFunctions_h_
#define _GenFunctions_h_

#include <map>
#include <vector>
#include <string>
#include <TMath.h>
#include <TVector3.h>
#include <TString.h>
#include <TH1.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <sstream>


void    MapLundIdToParticleNumber( map<int, int>    &mapLundIdToNumber  );

void    MapLundIdToParticleMass(   map<int, double> &mapLundIdToMass    );

void    MapLundIdToParticleName(   map<int, string> &mapLundIdToName    );

Bool_t  CheckParticle( vector<Int_t> _vec_lundId, Int_t _ParticleLundId, Int_t &_ParticleIndex);

Bool_t  CheckSIDIS( Double_t _q2, Double_t _w2, Double_t _y, Double_t _z_hadron );

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *_x, Double_t *_par) {

  return (0.5*_par[0]*_par[1]/TMath::Pi()) / 
    TMath::Max( 1.e-10,(_x[0]-_par[2])*(_x[0]-_par[2]) 
		+ .25*_par[1]*_par[1]);

}

// Gaussian Peak function
Double_t gaussianPeak(Double_t *_x, Double_t *_par) {


  return _par[0]*exp( -0.5*((_x[0] - _par[1])/_par[2])*((_x[0] - _par[1])/_par[2]) );

  /*  = p0*exp(-0.5*pow((x-p1)/p2),2) */

}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {

  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];

}

// pol3 background function
Double_t backgroundpol3(Double_t *x, Double_t *par) {

  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];

}

// Sum of a pol3 background and gaussian peak
Double_t fitGausBackPol3(Double_t *x, Double_t *par) {

  return backgroundpol3(x,par) + gaussianPeak(x,&par[4]);

}

// Sum of background and peak function
Double_t fitLorBackPol2(Double_t *x, Double_t *par) {

  return background(x,par) + lorentzianPeak(x,&par[3]);

}

#endif
