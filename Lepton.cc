#include "Lepton.h"


Lepton& Lepton::operator=(const Lepton&) = default;


//Object Identification
bool Lepton::mu_isLoose_tthlep(){
    // check whether muon pass tth loose lepton selection  
    bool isloose = false;
    if( pt > 5 && fabs(eta) <2.4 && fabs(dxy_pv)<=0.05 && fabs(dz_pv)<0.1 && IP3Dsig_it <8 
        && loose && miniIsoRel<0.4) isloose = true;
    return isloose;
};


void Lepton::set_Wp_tthlep(){
    // set lepton working point 
    if(mu_isLoose_tthlep())cut=1;// (1, loose) (2, medium) (3, tight)
}
