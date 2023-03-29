#include "../include/vars.h"

namespace vars {
  
  //Central e' momentum
  Double_t pcentral(Double_t ebeam, Double_t etheta, std::string Ntype) {
    Double_t temp = 0.;
    if (Ntype.compare("p") == 0) 
      temp = ebeam/(1. + (ebeam/physconst::Mp)*(1.0 - cos(etheta)));
    else if (Ntype.compare("n") == 0) 
      temp = ebeam/(1. + (ebeam/physconst::Mn)*(1.0 - cos(etheta)));
    else if (Ntype.compare("np") == 0) {
      Double_t Nmass = 0.5*(physconst::Mn + physconst::Mp);
      temp = ebeam/(1. + (ebeam/Nmass)*(1.0 - cos(etheta)));
    }
    else
      std::cerr << "Error: [vars::pcentral] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return temp;
  }

  //e' theta wrt exit beamline
  Double_t etheta(TLorentzVector Peprime) {
    return acos(Peprime.Pz() / Peprime.E());
  }

  //e' phi wrt exit beamline
  Double_t ephi(TLorentzVector Peprime) {
    return atan2(Peprime.Py(), Peprime.Px());
  }

  //rest nucleon four-vector
  void SetPN(std::string Ntype, TLorentzVector &PN) {
    if (Ntype.compare("p") == 0) 
      PN.SetPxPyPzE(0., 0., 0., physconst::Mp);
    else if (Ntype.compare("n") == 0) 
      PN.SetPxPyPzE(0., 0., 0., physconst::Mn);
    else if (Ntype.compare("np") == 0) 
      PN.SetPxPyPzE(0., 0., 0., 0.5*(physconst::Mn+physconst::Mp));
    else
      std::cerr << "Error: [vars::setPN] cannot be calculated. Enter a valid nucleon type." << std::endl;
  } 

  //scattered nucleon expected momentum
  Double_t pN_expect(Double_t nu, std::string Ntype) {
    if (Ntype.compare("p") == 0)                      
      return sqrt(pow(nu, 2.) + 2. * physconst::Mp * nu);
    else if (Ntype.compare("n") == 0)      
      return sqrt(pow(nu, 2.) + 2. * physconst::Mn * nu);
    else if (Ntype.compare("np") == 0)      
      return sqrt(pow(nu, 2.) + 2. * 0.5*(physconst::Mn+physconst::Mp) * nu);
    else {
      std::cerr << "Error: [vars::pN_expect] cannot be calculated. Enter a valid nucleon type." << std::endl;
      return -1;
    }
  }

  //virtual photon q vector
  TVector3 qVect_unit(Double_t Ntheta, Double_t Nphi) {
    TVector3 pNhat(sin(Ntheta) * cos(Nphi), sin(Ntheta) * sin(Nphi), cos(Ntheta));
    return pNhat;
  }

  //sets hcal axes by kinematic
  void SetHCALaxes(Double_t sbstheta_rad,                           // SBS angle in radian 
		   vector<TVector3> &HCAL_axes) {
    TVector3 HCAL_zaxis(sin(-sbstheta_rad),0,cos(-sbstheta_rad)); // Clock-wise rotation about Y axis
    TVector3 HCAL_xaxis(0,-1,0);                                  // -Y axis of Hall CoS = X axis of HCAL CoS
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();
    HCAL_axes.push_back(HCAL_xaxis);
    HCAL_axes.push_back(HCAL_yaxis);
    HCAL_axes.push_back(HCAL_zaxis);
  }

  //gets expected location of scatterd nucleon assuming straight line projections from BB track
  void GetxyHCALexpect(TVector3 vertex, TVector3 pNhat, TVector3 HCAL_origin, 
		       vector<TVector3> HCAL_axes, vector<Double_t> &xyHCALexpect) {
    /* This function calculates the expected vertical (x) and horizontal (y) positions
     of the recoil nucleon at the face of HCAL. 
     input:
     1. vertex       : vertex vector [in Hall CoS?? It must be but haven't confirmed], 
     2. pNhat        : projected q vector, 
     3. HCAL_origin  : HCAL origin vector [in Hall CoS], 
     4. HCAL_axes    : HCAL CoS axes [in Hall CoS]
     output:
     1. xyHCALexpect : expected x and y positions
    */
    // Intersection of a ray with a plane
    Double_t sintersect = (HCAL_origin - vertex).Dot(HCAL_axes[2]) / (pNhat.Dot(HCAL_axes[2]));
    // ray from Hall origin onto the face of HCAL where the nucleon hit
    TVector3 HCAL_intersect = vertex + sintersect*pNhat; 

    Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[0]);
    Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot(HCAL_axes[1]);

    xyHCALexpect.push_back(xexpect_HCAL);
    xyHCALexpect.push_back(yexpect_HCAL);
  }

  //calculates inverse momentum transfer Q^2
  Double_t Q2(Double_t ebeam, Double_t eeprime, Double_t etheta) {
    return 2.0*ebeam*eeprime*(1.0-cos(etheta));
  }

  //calculates invariant mass squared W^2
  Double_t W2(Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype) {
    Double_t temp = 0.;
    if (Ntype.compare("p") == 0) 
      temp = pow(physconst::Mp,2.0) + 2.0*physconst::Mp*(ebeam-eeprime) - Q2;
    else if (Ntype.compare("n") == 0) 
      temp = pow(physconst::Mn,2.0) + 2.0*physconst::Mn*(ebeam-eeprime) - Q2;
    else if (Ntype.compare("np") == 0) 
      temp = pow(0.5*(physconst::Mn+physconst::Mp),2.0) + 2.0*0.5*(physconst::Mn+physconst::Mp)*(ebeam-eeprime) - Q2;
    else
      std::cerr << "Error: [vars::W2] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return temp;
  }

  //calculates invariant mass W
  Double_t W(Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype) {
    return max(0., sqrt(vars::W2(ebeam, eeprime, Q2, Ntype)));
  }

  //calculates luminosity
  Double_t Luminosity(Double_t ibeam, std::string targetType) {
    Double_t lumi = 0.;
    if (targetType.compare("LH2") == 0)
      lumi = ((ibeam/physconst::qe)*econst::l_tgt*econst::lh2tarrho*(physconst::N_A/physconst::H2_Amass));
    else if (targetType.compare("LD2") == 0)
      lumi = ((ibeam/physconst::qe)*econst::l_tgt*econst::ld2tarrho*(physconst::N_A/physconst::D2_Amass));
    else
      std::cerr << "Error: [vars::Luminosity] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return lumi;
  }

} //::vars
