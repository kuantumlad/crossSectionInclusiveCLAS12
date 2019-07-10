/* ---------------------------------- */
/*    Frascati, March 17th, 2016.     */
/*                                    */
/*       author: Silvia Pisano        */
/*        pisanos@jlab.org            */
/*    silvia.pisano@lnf.infn.it       */
/* ---------------------------------- */
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <DisConstants.h>

using namespace std;

Double_t CalcPhih(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Phad);

Double_t CalcPhiR(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim);

Double_t CalcPhiRPerp(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim);

Double_t CalcPhiRPerpNew(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim);

Bool_t   CheckDIS(Double_t q2, Double_t w);

Bool_t   CheckClasDisComp(Double_t q2, Double_t w, Double_t xB, Double_t ElTheta);

void     DefDISVar(TLorentzVector InElec4Vector, TLorentzVector Elec4Vector, Double_t &Q2, Double_t &xB, Double_t &W, Double_t &y);

Double_t GetFeynmanX(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad);

Double_t GetZHarut(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad);

Double_t GetRapidity(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Phad);

Double_t GetPPerp(TLorentzVector Qgam, TLorentzVector Phad);

Double_t GetTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta);

Double_t GetSinTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta);

Double_t GetCosTheta(TLorentzVector Pgam, TLorentzVector Ptgt, TLorentzVector Ph, TLorentzVector Ppip, Int_t _iTheta);

/* Double_t CalcPhiRperpS(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Phad, TLorentzVector R); */

/* Double_t CalcPhiR2(TLorentzVector Qgam, TLorentzVector Target, TLorentzVector InElec, TLorentzVector Ppip, TLorentzVector Ppim); */

