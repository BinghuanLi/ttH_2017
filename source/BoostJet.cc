#include "../interference/BoostJet.h"


BoostJet& BoostJet::operator=(const BoostJet&) = default;


//Object Identification
void BoostJet::set_Wp_Top(int& numSoftTop, int& numLooseTop, int& numMediumTop, int& numTightTop){
    if( pt > 400 && fabs(eta) <2.4 && tau32 < 0.46 && fabs(softdrop_mass-162.5)<57.5){
        topCut = 4;
        numTightTop++; 
        numMediumTop++; 
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 400 && fabs(eta) <2.4 && tau32 < 0.54 && fabs(softdrop_mass-162.5)<57.5){
        topCut = 3;
        numMediumTop++; 
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 400 && fabs(eta) <2.4 && tau32 < 0.65 && fabs(softdrop_mass-162.5)<57.5){
        topCut = 2;
        numLooseTop++; 
        numSoftTop++; 
    }else if( pt > 400 && fabs(eta) <2.4 && tau32 < 0.8 && fabs(softdrop_mass-162.5)<57.5){
        topCut = 1;
        numSoftTop++; 
    }else{
        topCut = 0;
    }
}


void BoostJet::set_Wp_W(int& numLooseW, int& numTightW){
    if( pt > 190 && fabs(eta)<2.4 && tau21 < 0.45 && fabs(pruned_mass-85)<20){
        wCut = 2;
        numTightW++;
        numLooseW++;
    }else if( pt > 190 && fabs(eta)<2.4 && tau21 < 0.6 && fabs(pruned_mass-85)<20){
        wCut = 1;
        numLooseW++;
    }else{
        wCut = 0;
    }
};


void BoostJet::set_Wp_W(int& numSoftW, int& numLooseW, int& numMediumW, int& numTightW){
    if( pt > 190 && fabs(eta)<2.4 && tau21 < 0.45 && fabs(pruned_mass-85)<20){
        wCut = 4;
        numTightW++; 
        numMediumW++; 
        numLooseW++; 
        numSoftW++; 
    }else if( pt > 190 && fabs(eta)<2.4 && tau21 < 0.6 && fabs(pruned_mass-85)<20){
        wCut = 3;
        numMediumW++; 
        numLooseW++; 
        numSoftW++; 
    }else if( pt > 190 && fabs(eta)<2.4 && tau21 < 0.6){
        wCut = 2;
        numLooseW++; 
        numSoftW++; 
    }else if( pt > 190 && fabs(eta)<2.4){
        wCut = 1;
        numSoftW++; 
    }else{
        wCut = 0;
    }
};


void BoostJet::match_genW(vector<double>* genW_pt, vector<double>* genW_eta, vector<double>* genW_phi,
        vector<double>* genW_energy, vector<double>* genW_mass, vector<double>* genW_motherId){
    double min_diff_pt = 999.;
    for(uint gen_en=0; gen_en < genW_pt->size(); gen_en++){
        if(deltaR(deltaPhi(genW_phi->at(gen_en),phi),deltaEta(genW_eta->at(gen_en),eta))<0.8){
            min_diff_pt = fabs(genW_pt->at(gen_en) - pt);
            matchW_pt = genW_pt->at(gen_en); 
            matchW_eta = genW_eta->at(gen_en); 
            matchW_phi = genW_phi->at(gen_en); 
            matchW_energy = genW_energy->at(gen_en); 
            matchW_mass = genW_mass->at(gen_en); 
            matchW_mother_pdgId = genW_motherId->at(gen_en); 
        }
    }
};

