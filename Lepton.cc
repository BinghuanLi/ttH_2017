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
    bool istight = false;
    if (corrpt>10 && jetcsv <0.8484 && BDT > 0.9 && isMedium) istight = true;
    return istight;
};


void Lepton::set_Wp_tthlep(bool isMedium){
    // set lepton working point 
    //Muon
    if(fabs(pdgId)==13){
        if(mu_isTight_tthlep(isMedium) && mu_isfake_tthlep() && mu_isLoose_tthlep()) cut=3;// Tight
        else if(mu_isfake_tthlep() && mu_isLoose_tthlep()) cut=2;// fakeable
        else if(mu_isLoose_tthlep())cut=1;// loose
    }
}


void Lepton::cal_conept(bool isMedium){
    // set conept: BDT, pt and jetpt needs to be set before calling this function
    //Muon
    if(fabs(pdgId)==13){
        corrpt = (isMedium && BDT > 0.9) ?  pt : 0.9 * jetpt;
    }
}


void Lepton::cal_tight_property(){
    //Muon
    if(fabs(pdgId)==13){
        // set electron related variables to 1
        passConversion = 1;
        passEleTightCharge = 1;
        passMissHit = 1;
        if(fabs(pTErrOVpT_it)<0.2)passMuTightCharge=1;
        else passMuTightCharge=0;
    }
};
        
        
double Lepton::get_valX_valY_binContent(TH2F* h, double valX, double valY){
    int binX  = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(valX)));
    int binY = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindBin(std::abs(valY))));
    return h->GetBinContent(binX, binY);
};


