#include "TLorentzVector.h"

class Lepton {
    // this class define leptons
    public:
        // constructor
        Lepton() = default;
        Lepton(const Lepton&) = default;
        Lepton& operator=(const Lepton&);
        ~Lepton() = default;

        // member
        //readable variables 
        // reco variables
        double pt;
        double eta;
        double phi;
        double energy;
        double dxy_pv;
        double dz_pv;
        double IP3Dsig_it;
        double miniIsoRel;
        double charge;
        double isGlobal;
        double chi2;
        double chi2LocalPosition;
        double trkKink;
        double validFraction;
        double segmentCompatibility;
        double jetpt;
        double jetptratio;
        double jetcsv;
        double lepjetchtrks;
        double miniIsoCh;
        double miniIsoPUsub;
        double ptrel;
        double pTErrOVpT_it;
        double px;
        double py;
        double pz;
        double jetdr;
        double mvaValue_HZZ;
        
        // lepton id
        double loose;
        double pdgId;
        
        // new variables
        double cut; 
        double BDT;
        double corrpt;
        double FR;
        double CF;
        double passConversion;
        double passMuTightCharge;
        double passEleTightCharge;
        double passMissHit;
        double isMatchRightCharge;
        
        // lepton gen
        double gen_pt;
        double gen_eta;
        double gen_phi;
        double gen_en;
        double gen_pdgId;
        double genMother_pt;
        double genMother_eta;
        double genMother_phi;
        double genMother_en;
        double genMother_pdgId;
        double genGrandMother_pt;
        double genGrandMother_eta;
        double genGrandMother_phi;
        double genGrandMother_en;
        double genGrandMother_pdgId;
        double gen_isPromptFinalState;
        double gen_isDirectPromptTauDecayProductFinalState;

        // Object Identification
        bool mu_isLoose_tthlep();
        bool mu_isfake_tthlep();
        bool mu_isTight_tthlep(bool isMedium);
        void set_Wp_tthlep(bool isMedium);
        void cal_conept(bool isMedium);
        void cal_tight_property();

};
