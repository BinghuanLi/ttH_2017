#include "../interference/Tau.h"


Tau& Tau::operator=(const Tau&) = default;


//Object Identification
bool Tau::tau_isLoose_tthlep(){
    bool isloose = false;
    // check whether tau pass tth loose tau selection  
    if(pt >20 && fabs(eta)<2.3 && fabs(packedLeadTauCand_dz)<=0.2 
    && fabs(packedLeadTauCand_dxy)<=1000
    && byLooseIsolationMVArun2v1DBdR03oldDMwLT
    && decayModeFinding
    ) isloose = true;
    return isloose;
}


bool Tau::tau_isMedium_tthlep(){
    // check whether tau pass tth medium tau selection  
    bool ismedium = false;
    // check whether tau pass tth loose tau selection  
    if(pt >20 && fabs(eta)<2.3 && fabs(packedLeadTauCand_dz)<=0.2 
    && fabs(packedLeadTauCand_dxy)<=1000
    && byMediumIsolationMVArun2v1DBdR03oldDMwLT
    && decayModeFinding
    ) ismedium = true;
    return ismedium;
}
        
        
void Tau::set_Wp_tthlep(int& numLoose, int& numMedium){
    // set Tau working point
    if(tau_isLoose_tthlep() && tau_isMedium_tthlep()){
        cut = 2; // medium
        numMedium ++;
        numLoose ++;
    }else if(tau_isLoose_tthlep()){
        cut = 1;
        numLoose ++; // loose
    }
};
