#include "../interference/Jet.h"


Jet& Jet::operator=(const Jet&) = default;


//Object Identification
bool Jet::jet_isLoose(){
    bool isloose = false;
    if(pt > 25 && fabs(eta) <2.4 && neutralHadEnergyFraction<0.99 
       && neutralEmEnergyFraction<0.99 && numberOfConstituents>1 && chargedHadronEnergyFraction>0 
       && chargedEmEnergyFraction<0.99 && chargedMultiplicity>0) isloose = true;
    return isloose;
};
        
        
void Jet::set_Wp_jets(int& numLoose){
    if(jet_isLoose()){
        cut = 1.; // isloose pt 25
        numLoose++;
    } 
};


void Jet::set_Wp_bdisc(int& numbLoose, int& numbMedium, int& numbTight){
    if(pfCombinedInclusiveSecondaryVertexV2BJetTags>0.970){
        numbTight++;
        numbMedium++;
        numbLoose++;
        bcut = 3; // isTight b-tag
    }else if(pfCombinedInclusiveSecondaryVertexV2BJetTags>0.8484){
        numbMedium++;
        numbLoose++;
        bcut = 2; // isTight b-tag
    }else if(pfCombinedInclusiveSecondaryVertexV2BJetTags>0.5426){
        numbLoose++;
        bcut = 1; // isloose b-tag
    }else{
        bcut = 0;
    }
};


// hadTop tagger
void Jet::jet_ishadtop(int jet_index, int bjet_index, int w1_index, int w2_index){
    if(jet_index == bjet_index || jet_index == w1_index || jet_index == w2_index) isToptag=1;
    else isToptag = 0;
}
