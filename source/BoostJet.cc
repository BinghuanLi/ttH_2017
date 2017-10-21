#include "../interference/BoostJet.h"


BoostJet& BoostJet::operator=(const BoostJet&) = default;


//Object Identification
void BoostJet::set_Wp_Top(int& numSoftTop, int& numLooseTop, int& numMediumTop, int& numTightTop){
    if( pt > 200 && fabs(eta) <2.4 && tau32 < 0.46 && fabs(softdrop_mass-175)<25){
        topCut = 4;
        numTightTop++; 
        numMediumTop++; 
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 200 && fabs(eta) <2.4 && tau32 < 0.54 && fabs(softdrop_mass-175)<25){
        topCut = 3;
        numMediumTop++; 
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 200 && fabs(eta) <2.4 && tau32 < 0.65 && fabs(softdrop_mass-175)<25){
        topCut = 2;
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 200 && fabs(eta) <2.4 && tau32 < 0.8 && fabs(softdrop_mass-175)<25){
        topCut = 1;
        numSoftTop++; 
    }else{
        topCut = 0;
    }
}


void BoostJet::set_Wp_W(int& numLooseW, int& numTightW){
    if( pt > 200 && fabs(eta)<2.4 && tau21 < 0.45 && fabs(pruned_mass-85)<20){
        wCut = 2;
        numTightW++;
        numLooseW++;
    }else if( pt > 200 && fabs(eta)<2.4 && tau21 < 0.6 && fabs(pruned_mass-85)<20){
        wCut = 1;
        numLooseW++;
    }else{
        wCut = 0;
    }
};
