/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#ifndef _GenFunctions_cxx_
#define _GenFunctions_cxx_

#include <GenFunctions.h>


/* --------------------- */
/*     Check particles   */
/* --------------------- */
Bool_t  CheckParticle( vector<Int_t> _vec_lundId, Int_t _ParticleLundId, Int_t &_ParticleIndex)
{

  Bool_t _ThereIsTheParticle = kFALSE ;

  for( Int_t _iParticle = 0; _iParticle < (int)(_vec_lundId.size()); _iParticle++ ){

    if(      _vec_lundId.at( _iParticle ) == _ParticleLundId ){

      _ThereIsTheParticle  =  kTRUE ;

      _ParticleIndex       = _iParticle ;

      break ;

    }

  }

  return _ThereIsTheParticle;

}

/* ------------------------------ */
/*   Check if it is a SIDIS event */
/* ------------------------------ */
Bool_t  CheckSIDIS( Double_t _q2, Double_t _w2, Double_t _y, Double_t _z_hadron )
{

  Bool_t _IsSIDIS = kFALSE ;

  if( _q2 > 1 && _w2 > 4 && _y > 0.1 && _y < 0.85 && _z_hadron > 0.3 ) _IsSIDIS = kTRUE ;

  return _IsSIDIS;

}

/* --------------------- */
/*   Lund Id to mass map */
/* --------------------- */
void    MapLundIdToParticleNumber( map<int, int> &mapLundIdToNumber )
{

  mapLundIdToNumber[   11 ]  =  0 ; /* electron */
  mapLundIdToNumber[  211 ]  =  1 ; /* pi plus  */
  mapLundIdToNumber[ -211 ]  =  2 ; /* pi minus */
  mapLundIdToNumber[  321 ]  =  3 ; /* k plus   */
  mapLundIdToNumber[ -321 ]  =  4 ; /* k minus  */
  mapLundIdToNumber[ 2212 ]  =  5 ;
  mapLundIdToNumber[-2212 ]  =  6 ;
  mapLundIdToNumber[ 2112 ]  =  7 ;
  mapLundIdToNumber[-2112 ]  =  8 ;
  mapLundIdToNumber[  -11 ]  =  9 ;
  mapLundIdToNumber[   13 ]  = 10 ;
  mapLundIdToNumber[  -13 ]  = 11 ;
  mapLundIdToNumber[   12 ]  = 12 ;
  mapLundIdToNumber[   22 ]  = 13 ;
  mapLundIdToNumber[  111 ]  = 14 ;
  mapLundIdToNumber[  130 ]  = 15 ; /* 311 is k0, but it decays */
  // mapLundIdToNumber[  ]  =  ;

  return;

}
/* --------------------- */
/*   Lund Id to mass map */
/* --------------------- */
void    MapLundIdToParticleMass( map<int, double> &mapLundIdToMass )
{

  /* fill a map with the mass value associated to a specific particle, identified with its lund id */

  /* masses are given in GeV */


  /* -------- Gauge bosons --------- */


  /* photon */

  mapLundIdToMass[ 22 ]  = 0.0;


  /* --------- Leptons ----------- */


  /* electron */

  mapLundIdToMass[ 11 ]  = 0.0005 ;

  //   mapLundIdToMass[ 11 ]  = 0.000510998928 ;

  /* muon */

  mapLundIdToMass[ 13 ]  = 0.10565837 ;

  /* neutrino e- */

  mapLundIdToMass[ 12 ]  = 0.0 ;

  /* neutrino muonic */

  mapLundIdToMass[ 14 ]  = 0.0 ;


  /* --------- mesons ---------- */


  /* charged pion */

  mapLundIdToMass[ 211 ] = 0.13957 ;

  /* neutral pion */

  mapLundIdToMass[ 111 ] = 0.1349766 ;

  /* eta */

  mapLundIdToMass[ 221 ] = 0.547853 ;

  /* charged kaons */

  mapLundIdToMass[ 321 ] = 0.493677 ;

  /* neutral kaons */

  mapLundIdToMass[ 311 ] = 0.497614 ;

//   /* k long */

//   mapLundIdToMass[ 130 ] = 0.5 ; /* it is not the real value */


  /* ---------- baryons ---------- */


  /* protons */

  mapLundIdToMass[ 2212 ] = 0.938272046 ;

  /* neutrons */

  mapLundIdToMass[ 2112 ] = 0.939565379 ;


  return;

}

/* ----------------------- */
/*   Lund Id to name map   */
/* -------------------..-- */

void    MapLundIdToParticleName( map<int, string> &mapLundIdToName )
{

  /* fill a map with the mass value associated to a specific particle, identified with its lund id */

  /* masses are given in GeV */


  /* -------- Gauge bosons --------- */


  /* photon */

  mapLundIdToName[ 22 ]  = "photon";


  /* --------- Leptons ----------- */


  /* electron */

  mapLundIdToName[ 11 ]  = "electron" ;

  /* muon */

  mapLundIdToName[ 13 ]  = "muon" ;

  /* neutrino e- */

  mapLundIdToName[ 12 ]  = "neutrino (e^-)" ;

  /* neutrino muonic */

  mapLundIdToName[ 14 ]  =  "neutrino (#mu^-)" ;


  /* --------- mesons ---------- */


  /* charged pion */

  mapLundIdToName[ 211 ] = "charged pion" ;

  /* neutral pion */

  mapLundIdToName[ 111 ] = "neutral pion" ;

  /* eta */

  mapLundIdToName[ 221 ] = "eta" ;

  /* charged kaons */

  mapLundIdToName[ 321 ] = "charged kaon" ;

  /* neutral kaons */

  mapLundIdToName[ 311 ] = "neutral kaon" ;

  /* neutral kaons */

  mapLundIdToName[ 130 ] = "k long" ;

//   /* k long */

//   mapLundIdToName[ 130 ] = "k long" ; /* it is not the real value */


  /* ---------- baryons ---------- */


  /* protons */

  mapLundIdToName[ 2212 ] = "proton" ;

  /* neutrons */

  mapLundIdToName[ 2112 ] = "neutron" ;


  return;

}

#endif
