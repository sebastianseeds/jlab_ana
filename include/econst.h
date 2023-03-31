#ifndef ECONST_H
#define ECONST_H

#include <iostream>
#include <string>
#include <fstream>
#include "TMath.h"
#include "TString.h"

namespace econst {

  //// Detectors
  // HCAL 
  // Dimensions
  static const Int_t hcalchan = 288;
  static const Int_t hcalcol = 12;
  static const Int_t hcalrow = 24;
  static const Double_t hcalblk_w = 0.1524;       //m, width of a HCAL block
  static const Double_t hcalblk_h = 0.1524;       //m, height of a HCAL block
  static const Double_t hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
  static const Double_t hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
  static const Double_t hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const Double_t hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
  // Positions (mc)
  static const Double_t hcalposXi_mc = -2.3531;   //m, distance from beam center to top of HCal
  static const Double_t hcalposXf_mc = 1.45309;   //m, distance from beam center to bottom of HCal
  static const Double_t hcalposYi_mc = -0.931545; //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf_mc = 0.931545;  //m, distance from beam center to beam side of HCal
  // Positions (data)
  static const Double_t hcalposXi = -2.268095;    //m, distance from beam center to top of HCal
  static const Double_t hcalposXf = 1.538095;     //m, distance from beam center to bottom of HCal
  static const Double_t hcalposYi = -0.931545;    //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf = 0.931545;     //m, distance from beam center to beam side of HCal
  // Global
  static const Double_t hcalheight = -0.2897;      //m, height of the center of hcal above beam (m)

  // BBCAL
  static const Int_t shchan = 189;
  static const Int_t shcol = 7;
  static const Int_t shrow = 27;
  static const Int_t pschan = 52;
  static const Int_t pscol = 2;
  static const Int_t psrow = 26;  

  //// Experimental constants
  static const Int_t nkine = 6;
  // target
  static const Double_t l_tgt = 0.15;  //m 
  // LH2
  static const Double_t lh2tarrho = 0.0723;     //g/cc, target density
  static const Double_t lh2cthick = 0.02;       //cm, target cell thickness
  static const Double_t lh2uwallthick = 0.0145; //cm, upstream wall thickness
  static const Double_t lh2dwallthick = 0.015;  //cm, downstream wall thickness
  static const Double_t lh2dEdx = 0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
  // LD2
  static const Double_t ld2tarrho = 0.169;      //g/cc, target density
  // magnet
  static const Double_t sbsdipolegap = (48.0*2.54) / 100.;  // ~1.22 m
  static const Double_t sbsmaxfield = 3.1*atan(0.85 / (11.0-2.25-(sbsdipolegap/2.))) / (0.3*sbsdipolegap*0.7); // ~1.26 T (?)

  // shielding
  static const Double_t rho_al = 2.7; //g/cc
  static const Double_t aldEdx = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV

  // Following quantities vary with configuration
  Double_t     ebeam(Int_t config);        //GeV
  Double_t     bbtheta(Int_t config);      //deg
  Double_t     bbdist(Int_t config);       //m
  Double_t     sbstheta(Int_t config);     //deg
  Double_t     sbsdist(Int_t config);      //m
  Double_t     hcaltheta(Int_t config);    //deg
  Double_t     hcaldist(Int_t config);     //m
}

// a class for SBS config
class SBSconfig {
 public:

  Int_t        GetSBSconf()       const { return fSBSconf; }
  Int_t        GetSBSmag()        const { return fSBSmag; }
  Double_t     GetEbeam()         const { return fEbeam; }
  Double_t     GetBBtheta()       const { return fBBtheta; }
  Double_t     GetBBtheta_rad()   const { return fBBtheta_rad; }
  Double_t     GetBBdist()        const { return fBBdist; }
  Double_t     GetSBStheta()      const { return fSBStheta; }
  Double_t     GetSBStheta_rad()  const { return fSBStheta_rad; }
  Double_t     GetSBSdist()       const { return fSBSdist; }
  Double_t     GetHCALtheta()     const { return fHCALtheta; }
  Double_t     GetHCALtheta_rad() const { return fHCALtheta_rad; }
  Double_t     GetHCALdist()      const { return fHCALdist; }

  // constructor
  SBSconfig(Int_t conf, Int_t sbsmag) {
    fSBSconf       = conf;
    fSBSmag        = sbsmag;
    fEbeam         = econst::ebeam(conf);
    fBBtheta       = econst::bbtheta(conf);
    fBBtheta_rad   = econst::bbtheta(conf)*TMath::DegToRad();
    fBBdist        = econst::bbdist(conf);
    fSBStheta      = econst::sbstheta(conf);
    fSBStheta_rad  = econst::sbstheta(conf)*TMath::DegToRad();
    fSBSdist       = econst::sbsdist(conf);
    fHCALtheta     = econst::hcaltheta(conf);
    fHCALtheta_rad = econst::hcaltheta(conf)*TMath::DegToRad();
    fHCALdist      = econst::hcaldist(conf);
  }

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const SBSconfig& sbsconf) {
    out  << " -------------------------- "                                             << std::endl
	 << Form(" SBS Config: %d, "                   , sbsconf.fSBSconf)             << std::endl
	 << Form(" SBS Magnet Settings: %d (p), "      , sbsconf.fSBSmag)              << std::endl
    	 << Form(" Beam energy: %0.4f (GeV),"          , sbsconf.fEbeam)               << std::endl
    	 << Form(" BigBite angle: %0.1f (deg),"        , sbsconf.fBBtheta)             << std::endl
      	 << Form(" BigBite distance: %0.5f (m),"       , sbsconf.fBBdist)              << std::endl
    	 << Form(" Super BigBite angle: %0.1f (deg),"  , sbsconf.fSBStheta)            << std::endl
      	 << Form(" Super BigBite distance: %0.2f (m)," , sbsconf.fSBSdist)             << std::endl
	 << Form(" HCAL angle: %0.1f (deg),"           , sbsconf.fHCALtheta)           << std::endl
    	 << Form(" HCAL distance: %0.1f (m)"           , sbsconf.fHCALdist)            << std::endl
	 << " -------------------------- "                        << std::endl         << std::endl;
    return out;
  }

 private:
  Int_t        fSBSconf;             // SBS configuration number
  Int_t        fSBSmag;              // SBS magnet settings (%)
  Double_t     fEbeam;               // beam energy (better to get this from tree) (GeV)
  Double_t     fBBtheta;             // BigBite magnet angle (deg)
  Double_t     fBBtheta_rad;         // BigBite magnet angle (rad)
  Double_t     fBBdist;              // BigBite magnet distance from target (m)
  Double_t     fSBStheta;            // Super BigBite magnet angle (deg)
  Double_t     fSBStheta_rad;        // Super BigBite magnet angle (rad)
  Double_t     fSBSdist;             // Super BigBite magnet distance from target (m)
  Double_t     fHCALtheta;           // HCAL angle (deg)
  Double_t     fHCALtheta_rad;       // HCAL angle (rad)
  Double_t     fHCALdist;            // HCAL distance from target (m)
};

#endif
