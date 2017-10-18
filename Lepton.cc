#include "Lepton.h"


Lepton& Lepton::operator=(const Lepton&) = default;


//Object Identification
bool Lepton::mu_isLoose_tthlep(){
    // check whether muon pass tth loose lepton selection  
    bool isloose = false;
    if( pt > 5 && fabs(eta) <2.4 && fabs(dxy_pv)<=0.05 && fabs(dz_pv)<0.1 && IP3Dsig <8 
        && loose && miniIsoRel<0.4) isloose = true;
    return isloose;
};


bool Lepton::mu_isfake_tthlep(){
    // check whether muon pass tth fakable lepton selection  
    bool isfakeable = false;
    if (corrpt <= 10)return isfakeable;
    if (BDT > 0.90){
        if (jetcsv <0.8484 ) isfakeable = true;
    }
    else{
        if (jetcsv <0.3 && jetptratio > 0.5 && segmentCompatibility > 0.3)isfakeable = true;
    }
    return isfakeable;
};


bool Lepton::mu_isTight_tthlep(bool isMedium){
    // check whether muon pass tth tight lepton selection  
    bool istight = false;
    if (corrpt>10 && jetcsv <0.8484 && BDT > 0.9 && isMedium) istight = true;
    return istight;
};


bool Lepton::ele_isLoose_tthlep(){
    // check whether ele pass tth loose lepton selection  
    bool isloose = false;
    if(pt>7 && fabs(eta)<2.5 && fabs(dxy_pv)<=0.05 && fabs(dz_pv)<0.1 && IP3Dsig<8 && miniIsoRel<0.4 && fabs(expectedMissingInnerHits)<=1 ) isloose = true;
    return isloose;
}


bool Lepton::ele_isfake_tthlep(){
    // check whether ele pass tth fakeable lepton selection  
    bool isfakeable = false;
    if (corrpt <= 10)return isfakeable;
    bool eleMVAId = false;
    if((fabs(SCeta) <0.8 && mvaValue_HZZ > 0.0)||
        (0.8 <= fabs(SCeta) && fabs(SCeta) <1.479 && mvaValue_HZZ >0.0)||
        (1.479 <=fabs(SCeta) && fabs(SCeta)<500 && mvaValue_HZZ>0.7)
     ){
        eleMVAId = true;
    } 
    if (BDT > 0.9){
        if (ele_passCuts() && jetcsv <0.8484) isfakeable = true;
    }
    else{
        if (ele_passCuts() && jetcsv <0.3 && jetptratio > 0.5 && eleMVAId)isfakeable = true;
    }
    return isfakeable;
};


bool Lepton::ele_isTight_tthlep(){
    // check whether ele pass tth tight lepton selection  
    bool istight = false;
    if (ele_passCuts() && corrpt>10 && jetcsv <0.8484 && BDT > 0.9) istight = true;
    return istight;
};


bool Lepton::ele_passCuts(){
    bool passCuts = false;
    if (fabs(SCeta) < 0.8){
        passCuts = full5x5_sigmaIetaIeta < 0.011 &&
        hOverE < 0.10 &&
        fabs(dEtaIn) < 0.01 &&
        fabs(dPhiIn) < 0.04 &&
        ooEmooP > -0.05 &&
        ooEmooP < 0.010;
    }
    else if (fabs(SCeta) < 1.479) {
        passCuts = full5x5_sigmaIetaIeta < 0.011 &&
        hOverE < 0.10 &&
        fabs(dEtaIn) < 0.01 &&
        fabs(dPhiIn) < 0.04 &&
        ooEmooP > -0.05 &&
        ooEmooP < 0.010;
    }
    else if (fabs(SCeta) < 2.5) {
        passCuts = full5x5_sigmaIetaIeta < 0.030 &&
        hOverE < 0.07 &&
        fabs(dEtaIn) < 0.008 &&
        fabs(dPhiIn) < 0.07 &&
        ooEmooP > -0.05 &&
        ooEmooP < 0.005;
    }
    else passCuts = true;
    return passCuts;
};


void Lepton::set_Wp_tthlep(bool isMedium, int& numLoose, int& numFake, int& numtight){
    // set lepton working point 
    //Muon
    if(fabs(pdgId)==13){
        if(mu_isTight_tthlep(isMedium) && mu_isfake_tthlep() && mu_isLoose_tthlep()){ 
            cut=3;  // Tight
            numtight ++;
            numFake ++;
            numLoose ++;
        }else if(mu_isfake_tthlep() && mu_isLoose_tthlep()){ 
            cut=2;// fakeable
            numFake ++;
            numLoose ++;
        }else if(mu_isLoose_tthlep()){
            cut=1;
            numLoose ++;
            // loose
        }
    }
    //Electron
    if(fabs(pdgId)==11){
        if(ele_isTight_tthlep() && ele_isfake_tthlep() && ele_isLoose_tthlep()){ 
            cut=3;  // Tight
            numtight ++;
            numFake ++;
            numLoose ++;
        }else if(ele_isfake_tthlep() && ele_isLoose_tthlep()){ 
            cut=2;// fakeable
            numFake ++;
            numLoose ++;
        }else if(ele_isLoose_tthlep()){
            cut=1;
            numLoose ++;
            // loose
        }
    }
}


void Lepton::cal_conept(bool isMedium){
    // set conept: BDT, pt and jetpt needs to be set before calling this function
    //Muon
    if(fabs(pdgId)==13){
        corrpt = (isMedium && BDT > 0.9) ?  pt : 0.9 * jetpt;
    }
    //Electron
    if(fabs(pdgId)==11){
        corrpt = BDT > 0.9 ?  pt : 0.9 * jetpt;
    }
}


void Lepton::cal_tight_property(){
    // calculate tigit lepton variables
    //Muon
    if(fabs(pdgId)==13){
        // set electron related variables to 1
        passConversion = 1;
        passEleTightCharge = 1;
        passMissHit = 1;
        if(fabs(pTErrOVpT_it)<0.2) passMuTightCharge= 1;
        else passMuTightCharge = 0;
    }
    //Electron
    if(fabs(pdgId)==11){
        // set muon related variables to 1
        passMuTightCharge=1;
        if(expectedMissingInnerHits==0) passMissHit = 1;
        else passMissHit = 0;
        if((isGsfCtfScPixChargeConsistent + isGsfScPixChargeConsistent)>1){
             passEleTightCharge = 1;
        }else{
             passEleTightCharge = 0;
        }
    }
};
        
        
double Lepton::get_valX_valY_binContent(TH2F* h, double valX, double valY){
    // return TH2F bin content of bin value (valX, valY) 
    int binX  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(valX)));
    int binY = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(std::abs(valY))));
    return h->GetBinContent(binX, binY);
};


