/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#include "DisFunctions.h"

/* Phih */
Double_t CalcPhih(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Phad)
{

  Double_t Phihad;

  TVector3 Boost = -(Qgam + Target).BoostVector();
  Double_t c0, c1, c2, c3;
  TVector3 v0, v1;

  /*Boost the vectors*/
  Qgam.Boost(Boost);
  InElec.Boost(Boost);
  Phad.Boost(Boost);
  Target.Boost(Boost);

  /* Calcolo phi_h */
  v0 = Qgam.Vect().Cross(InElec.Vect());
  v1 = Qgam.Vect().Cross(Phad.Vect());
  c0 = v0.Dot(Phad.Vect());
  c1 = v0.Dot(v1);
  c2 = v0.Mag();
  c3 = v1.Mag();

  Phihad = AngleConversion*(c0/TMath::Abs(c0)) * TMath::ACos(c1 / (c2*c3) );

  return Phihad;
}

/* PhiR */
Double_t CalcPhiR(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim)
{

  TVector3 Boost = -(Qgam + Target).BoostVector();

  /*Boost the vectors*/
  Ppip.Boost(Boost);
  Ppim.Boost(Boost);
  Qgam.Boost(Boost);
  InElec.Boost(Boost);
  Target.Boost(Boost);

  TLorentzVector _Phad = Ppip + Ppim;
  TLorentzVector _R    = 0.5*(Ppip - Ppim);

  /* Now define the 3-vectors */
  TVector3 _PElein = InElec.Vect() ;
  TVector3 _Pgvec  = Qgam.Vect()   ;
  TVector3 _Rvec   = _R.Vect()      ;
  TVector3 _Phvec  = _Phad.Vect()   ;

  /* versor n */
  TVector3 v1 = _Pgvec.Cross(_PElein);
  TVector3 A0 = v1.Unit();
  TVector3 v0 = _Phvec.Cross(A0);
  TVector3 B0 = (v0.Cross(A0)).Unit();
  Double_t c0 = v1.Dot(_Phvec);
  Double_t c1 = v1.Mag();
  Double_t c2 = _Phvec.Mag();

  Double_t delta = 0.5*TMath::Pi() - TMath::ACos(c0 / (c1 * c2));

  Double_t c3 = TMath::Cos(delta);
  Double_t c4 = TMath::Sin(delta);

  TVector3 N0 = A0 * c3 + B0 * c4;

  /* R_T */
  Double_t _long_comp = _Rvec.Dot(_Phvec) ;
  Double_t _phmag     =  TMath::Abs(1.0/(_Phvec.Mag()*_Phvec.Mag()));
  TVector3 RhadL      = _Phvec*_long_comp*_phmag;
  TVector3 RhadT      = _Rvec - RhadL;

  /* PhiR2 */

  Double_t _f1   = N0.Dot(RhadT);
  TVector3 _vec1 = _Phvec.Cross(RhadT);
  Double_t _f2   = N0.Dot( _vec1 );
  Double_t _f3   = _vec1.Mag();

  Double_t PhiR = AngleConversion*(_f1/TMath::Abs(_f1)) * TMath::ACos( _f2 /_f3 );

  return PhiR;

}

Double_t CalcPhiRPerp(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim)
{

  TVector3 Boost = -(Qgam + Target).BoostVector();

  /*Boost the vectors*/
  Ppip.Boost(Boost);
  Ppim.Boost(Boost);
  Qgam.Boost(Boost);
  InElec.Boost(Boost);
//   Target.Boost(Boost);

  TLorentzVector _Phad = Ppip + Ppim;
  TLorentzVector _R    = 0.5*(Ppip - Ppim);

  /* Now define the 3-vectors */

  TVector3 _PElein =  InElec.Vect() ;
  TVector3 _Pgvec  =  Qgam.Vect()   ;
  TVector3 _Rvec   = _R.Vect()      ;
  TVector3 _Phvec  = _Phad.Vect()   ;

    /* R_T */
  Double_t _long_comp = _Rvec.Dot(_Phvec) ;
  Double_t _phmag     =  TMath::Abs(1.0/(_Phvec.Mag()*_Phvec.Mag()));
  TVector3 RhadL      = _Phvec*_long_comp*_phmag;
  TVector3 RhadT      = _Rvec - RhadL;

  TVector3 _vec1 = _Pgvec.Cross(_PElein);

  TVector3 _vec2 = _Pgvec.Cross(RhadT);

  Double_t _c1   = _vec1.Dot(RhadT);

  Double_t _c2   = _vec1.Dot( _vec2 );

  Double_t _c3   = _vec1.Mag();

  Double_t _c4   = _vec2.Mag();

  // cout << " Pgvec " << _Pgvec.Mag() << endl ;

  // cout << " _c1 = " << _c1 << endl ;

  // cout << " _c2 = " << _c2 << endl ;

  // cout << " _c3 = " << _c3 << endl ;

  // cout << " _c4 = " << _c4 << endl ;




  Double_t PhiRPerp = AngleConversion * (_c1/TMath::Abs(_c1)) * TMath::ACos(_c2 / (_c3*_c4) );

  return PhiRPerp;

}

Double_t CalcPhiRPerpNew(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim)
{

  TVector3 Boost = -(Qgam + Target).BoostVector();

  /*Boost the vectors*/
  Ppip.Boost(Boost);
  Ppim.Boost(Boost);
  Qgam.Boost(Boost);
  InElec.Boost(Boost);
//   Target.Boost(Boost);

  TLorentzVector _Phad = Ppip + Ppim;
  TLorentzVector _R    = 0.5*(Ppip - Ppim);

  /* Now define the 3-vectors */

  TVector3 _PElein = InElec.Vect() ;
  TVector3 _Pgvec  = Qgam.Vect()   ;
//   TVector3 _Rvec   = _R.Vect()      ;
//   TVector3 _Phvec  = _Phad.Vect()   ;

    /* R_T */
//   Double_t _long_comp = _Rvec.Dot(_Phvec) ;
//   Double_t _phmag     =  TMath::Abs(1.0/(_Phvec.Mag()*_Phvec.Mag()));
//   TVector3 RhadL      = _Phvec*_long_comp*_phmag;
//   TVector3 RhadT      = _Rvec - RhadL;

  Double_t       _long_comp_4d = (Double_t)( _R.Dot(_Phad) );
  Double_t       _phmag_4d     = TMath::Abs(1.0/(_Phad.Mag()*_Phad.Mag()));
  TLorentzVector  RhadL_4d     = _Phad*_long_comp_4d*_phmag_4d;
  TLorentzVector  RhadT_4d     = _R - RhadL_4d;

  TVector3        RhadT        = RhadT_4d.Vect();


  TVector3 _vec1 = _Pgvec.Cross(_PElein);

  TVector3 _vec2 = _Pgvec.Cross(RhadT);

  Double_t _c1   = _vec1.Dot(RhadT);

  Double_t _c2   = _vec1.Dot( _vec2 );

  Double_t _c3   = _vec1.Mag();

  Double_t _c4   = _vec2.Mag();

  Double_t PhiRPerp = AngleConversion * (_c1/TMath::Abs(_c1)) * TMath::ACos(_c2 / (_c3*_c4) );

  return PhiRPerp;

}

void DefDISVar(TLorentzVector InElec4Vector, TLorentzVector Elec4Vector, Double_t &Q2, Double_t &xB, Double_t &W, Double_t &y)
{
  Double_t Eprime, ThetaEl, Ebeam;
  Eprime  = Elec4Vector.E();
  ThetaEl = Elec4Vector.Theta();
  Ebeam   = InElec4Vector.E();

  Q2 = 4*Ebeam*Eprime*TMath::Sin(ThetaEl/2)*TMath::Sin(ThetaEl/2);
  xB = Q2/(2*ProtonMass*(Ebeam-Eprime));
  W  = sqrt(ProtonMass*ProtonMass+2*(Ebeam-Eprime)*ProtonMass-Q2);
  y  = (InElec4Vector-Elec4Vector).E() / InElec4Vector.E();

}

Double_t GetPPerp(TLorentzVector Qgam, TLorentzVector Phad)
{

  Double_t _pt = -999.0;

  TVector3 pgamma3 = Qgam.Vect();
  TVector3 phad3   = Phad.Vect();

  _pt = phad3.Perp(pgamma3);

  return _pt;

}

Bool_t CheckDIS(Double_t q2, Double_t w)
{

  Bool_t iDIS = kFALSE;

  if(q2 > 1 && w > 2) iDIS = kTRUE;

  return iDIS;

}

Bool_t CheckClasDisComp(Double_t q2, Double_t w, Double_t xB, Double_t ElTheta)
{
  Bool_t iCompCutsOk=kFALSE;

  if((ElTheta>14)    
     &&(q2>0.85)                                                                                                                                                                                                                   
     &&(q2<20)                                                                                                                                                                                                                          
     &&(w>2)                                                                                                                                                                                                                          
     &&(w<TMath::Sqrt(50))                                                                                                                                                                                                              
     &&(xB > 0.05)                                                                                                                                                                                                                      
     &&(xB < 0.95)                                                                                                                                                                                                                      
     //     &&()                                                                                                                                                                                                                        
     ) iCompCutsOk=kTRUE;

  return iCompCutsOk;
}

Double_t GetFeynmanX(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad) {

  TVector3 boost = -(Pgam+Ptgt).BoostVector();

  Pgam.Boost(boost);
  Ptgt.Boost(boost);
  Phad.Boost(boost);

  /* Max hadron momentum */
  // Double_t mX2 = (Pgam+Ptgt-Phad).Mag2();
  // Double_t mH2 = Phad.Mag2();
  // Double_t Ein2 = (Pgam+Ptgt).E() * (Pgam+Ptgt).E();
  // Double_t Pmax = TMath::Sqrt( ( (Ein2-mX2-mH2)*(Ein2-mX2-mH2) - 4.*mX2*mH2 ) / (4.*Ein2) );

  /* xF definition according to Kotzinian */
  Double_t W = (Pgam+Ptgt).M();
  Double_t xF = 2. * Pgam.Vect().Dot( Phad.Vect() ) / (W*Pgam.Vect().Mag());

  /* xF definition according to Hermes */
  //Double_t xF = Pgam.Vect().Dot( Phad.Vect() ) / (Pmax*Pgam.Vect().Mag());                                                                                                                                                               

  /* Using cos(theta_lambda_cm) instead of xF */
  //Double_t xF = Phad.CosTheta();                                                                                                                                                                                                         
  return xF;
}

Double_t GetZHarut(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad) {

  TVector3 boost = -(Pgam+Ptgt).BoostVector();

  Pgam.Boost(boost);
  Ptgt.Boost(boost);
  Phad.Boost(boost);

  /* z definition according to Harut snapshot */
  Double_t _z_num = Phad.E() + (Phad.Vect()).Dot(Pgam.Vect())/((Pgam.Vect()).Mag());
  Double_t _z_den = Pgam.E() + Pgam.Pz();
  Double_t _zh    = _z_num/_z_den;

  return _zh;
}

Double_t GetRapidity(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad) {

  TVector3 boost = -(Pgam+Ptgt).BoostVector();

  // Pgam.Boost(boost);
  // Ptgt.Boost(boost);
  Phad.Boost(boost);

  Double_t eta = Phad.Rapidity();

  return eta;
}

Double_t GetTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta){

  /* the last argument, _iTheta, tells what Theta to use, if the standard one or the one symmetrized with respect to pi/2 */
  /* if _iTheta = 0 -> use the standard theta; if _iTheta = 1 use the symmetric one */

  Double_t _theta = -99999.0;

  TVector3 boost_gp = -(Pgam+Ptgt).BoostVector();

  TVector3 boost_ph = -(Ph).BoostVector();

  /* now boost the vectors */
  Ppip.Boost(boost_ph);

  Ph.Boost(boost_gp);

  if(      !_iTheta      ) _theta = (Ppip.Vect()).Angle(Ph.Vect());
  else if(  _iTheta == 1 ) _theta = TMath::Abs( TMath::Abs(((Ppip.Vect()).Angle(Ph.Vect())) - TMath::Pi()/2.) - TMath::Pi()/2. ); /* theta_prime is symmetric with respect to pi/2 */

  return _theta;

}

Double_t GetSinTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta){

  /* the last argument, _iTheta, tells what Theta to use, if the standard one or the one symmetrized with respect to pi/2 */
  /* if _iTheta = 0 -> use the standard theta; if _iTheta = 1 use the symmetric one */

  Double_t _sin_theta = -99999.0;

  TVector3 boost_gp = -(Pgam+Ptgt).BoostVector();

  TVector3 boost_ph = -(Ph).BoostVector();

  /* now boost the vectors */
  Ppip.Boost(boost_ph);

  Ph.Boost(boost_gp);

  Double_t _theta       = (Ppip.Vect()).Angle(Ph.Vect());

  Double_t _theta_prime = TMath::Abs( TMath::Abs(_theta - TMath::Pi()/2.) - TMath::Pi()/2. ); /* theta_prime is symmetric with respect to pi/2 */

  if(      !_iTheta      ) _sin_theta = TMath::Sin( _theta       );
  else if(  _iTheta == 1 ) _sin_theta = TMath::Sin( _theta_prime );
  else                      cout << " WARNING! Wrong theta request! " << endl;

  return _sin_theta;

}

Double_t GetCosTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta){

  /* the last argument, _iTheta, tells what Theta to use, if the standard one or the one symmetrized with respect to pi/2 */
  /* if _iTheta = 0 -> use the standard theta; if _iTheta = 1 use the symmetric one */

  Double_t _cos_theta = -99999.0;

  TVector3 boost_gp = -(Pgam+Ptgt).BoostVector();

  TVector3 boost_ph = -(Ph).BoostVector();

  /* now boost the vectors */
  Ppip.Boost(boost_ph);

  Ph.Boost(boost_gp);

  Double_t _theta       = (Ppip.Vect()).Angle(Ph.Vect());

  Double_t _theta_prime = TMath::Abs( TMath::Abs(_theta - TMath::Pi()/2.) - TMath::Pi()/2. ); /* theta_prime is symmetric with respect to pi/2 */

  if(      !_iTheta      ) _cos_theta = TMath::Cos( _theta       );
  else if(  _iTheta == 1 ) _cos_theta = TMath::Cos( _theta_prime );
  else                      cout << " WARNING! Wrong theta request! " << endl;

  return _cos_theta;

}
// /* PhiR by Van der Nat */
// Double_t CalcPhiRVdN(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Phad, TLorentzVector R)
// {

//   Double_t PhiR = -9999999.9 ;

//   TVector3 Boost = -(Qgam + Target).BoostVector();

//   /*Boost the vectors*/
//   Qgam.Boost(Boost);
//   InElec.Boost(Boost);
//   Phad.Boost(Boost);
//   Target.Boost(Boost);

//   /* Now define the 3-vectors */
//   TVector3 _PElein = InElec.Vect() ;
//   TVector3 _Pgvec  = Qgam.Vect()   ;
//   TVector3 _Rvec   = R.Vect()      ;
//   TVector3 _Phvec  = Phad.Vect()   ;

//   /* define n */
//   TVector3 _Aver   = _Pgvec.Cross( _PElein );

//   Double_t _avmag  = 1/_Pgvec.Mag()*_PElein.Mag();

//   _Aver            = _Aver*_avmag;

//   TVector3 _Bver   = ( _Phvec.Cross( _Aver ) ).Cross( _Aver );

//   Double_t _bvmag  = 1/_Phvec.Mag();

//   _Bver            = _Bver*_bvmag;

//   Double_t _delta  = 0.5*TMath::Pi() - TMath::ACos( (_Pgvec.Cross( _PElein ) ).Dot( _Phvec )*_avmag*_bvmag );

//   TVector3 _n      = _Aver*( TMath::Cos(_delta) ) + _Bver*(TMath::Sin(_delta) );

//   /* R_T */
//   Double_t _long_comp = _Rvec.Dot(_Phvec.Unit()) ;
//   TVector3 RhadL      = (_Phvec.Unit())*_long_comp;
//   TVector3 RhadT      = _Rvec - RhadL;

//   /* PhiR */

//   Double_t _f1 = _n.Dot( RhadT ) ;

//   Double_t _s1 = RhadT.Mag();

//   Double_t _f2 = _n.Dot( _Phvec.Cross( RhadT ) );

//   Double_t _s2 = ( _Phvec.Cross( RhadT ) ).Mag();

//   PhiR = AngleConversion*( _f1/TMath::Abs(_s1) ) * TMath::ACos( _f2/TMath::Abs(_s2) );

//   return PhiR;

// }

/* PhiRperp */


// /* PhiR */
// Double_t CalcPhiR(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim)
// {

//   TVector3 Boost = -(Qgam + Target).BoostVector();
//   Double_t c0, c1, c2, c3;
//   TVector3 v0, v1;

//   /*Boost the vectors*/
//   Ppip.Boost(Boost);
//   Ppim.Boost(Boost);
//   Qgam.Boost(Boost);
//   InElec.Boost(Boost);
// //   Target.Boost(Boost);

//   TLorentzVector _Phad = Ppip + Ppim;
//   TLorentzVector _R    = 0.5*(Ppip - Ppim);

//   /* Now define the 3-vectors */
//   TVector3 _PElein = InElec.Vect() ;
//   TVector3 _Pgvec  = Qgam.Vect()   ;
//   TVector3 _Rvec   = _R.Vect()      ;
//   TVector3 _Phvec  = _Phad.Vect()   ;

//   /* R_T */
//   Double_t _long_comp = _Rvec.Dot(_Phvec) ;
//   Double_t _phmag     =  TMath::Abs(1.0/(_Phvec.Mag()*_Phvec.Mag()));
//   TVector3 RhadL      = _Phvec*_long_comp*_phmag;
//   TVector3 RhadT      = _Rvec - RhadL;

//   /* PhiRPerp */

//   v0 = _Pgvec.Cross( _PElein );

//   v1 = _Pgvec.Cross(  RhadT  );

//   c0 =  v0.Dot(       RhadT  );

//   c1 =  v0.Dot( v1 );

//   c2 =  v0.Mag();

//   c3 =  v1.Mag();

//   Double_t PhiR = AngleConversion*(c0/TMath::Abs(c0)) * TMath::ACos(c1 / (c2*c3) );

//   return PhiR;

// }


// /* PhiRperp by Sergio */

// Double_t CalcPhiRperpS(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Phad, TLorentzVector R)
// {

//   Double_t PhiRPerp = -9999999.9 ;

//   TVector3 Boost = -(Qgam + Target).BoostVector();
//   TVector3 v0, v1;

//   /*Boost the vectors*/
//   Qgam.Boost(Boost);
//   InElec.Boost(Boost);
//   Phad.Boost(Boost);
//   Target.Boost(Boost);

//   /* Now define the 3-vectors */

//   TVector3 _PElein = InElec.Vect() ;
//   TVector3 _Pgvec  = Qgam.Vect()   ;
//   TVector3 _Rvec   = R.Vect()      ;
//   TVector3 _Phvec  = Phad.Vect()   ;

//   /* R_T */
//   Double_t _long_comp = _Rvec.Dot(_Phvec.Unit()) ;
//   TVector3 RhadL      = (_Phvec.Unit())*_long_comp;
//   TVector3 RhadT      = _Rvec - RhadL;

//   /* Define n */

//   Double_t _s0 = (_Rvec.Cross(_Phvec)).Dot(_PElein);

//   Double_t _f1 = ( _Pgvec.Cross( _PElein ) ).Dot( _Phvec.Cross(RhadT) );

//   Double_t _s1 =  (_Rvec.Cross(_Phvec)).Mag();

//   Double_t _s2 =  (_Phvec.Cross(RhadT)).Mag();

//   PhiRPerp     = _s0/(TMath::Abs( _s0 )) * TMath::ACos( _f1/( _s1*_s2) )*AngleConversion ;


//   return PhiRPerp;

// }

