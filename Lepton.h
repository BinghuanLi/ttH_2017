#include "TH2F.h"

class Lepton {
    // this class define leptons
    public:
        // constructor
        Lepton() = default;
        Lepton(const Lepton&) = default;
        Lepton& operator=(const Lepton&);
        ~Lepton() = default;

        // member
        // readable variables 
        // common variables
        // reco variables
        double pt = -999.;
        double eta = -999.;
        double phi = -999.;
        double energy = -999.;
        double dxy_pv = -999.;
        double dz_pv = -999.;
        double IP3Dsig = -999.;
        double miniIsoRel = -999.;
        double charge = -999.;
        double jetpt = -999.;
        double jetptratio = -999.;
        double jetcsv = -999.;
        double lepjetchtrks = -999.;
        double miniIsoCh = -999.;
        double miniIsoPUsub = -999.;
        double ptrel = -999.;
        double px = -999.;
        double py = -999.;
        double pz = -999.;
        double jetdr = -999.;
        
        // lepton id
        double loose = -999.;
        double pdgId = -999.;
        
        // Muon variables
        double isGlobal = -999.;
        double chi2 = -999.;
        double chi2LocalPosition = -999.;
        double trkKink = -999.;
        double validFraction = -999.;
        double segmentCompatibility = -999.;
        double pTErrOVpT_it = -999.;
        
        // Electron variables
        double mvaValue_HZZ = -999.; 
        double SCeta = -999.;
        double expectedMissingInnerHits = -999.;
        double full5x5_sigmaIetaIeta = -999.;
        double hOverE = -999.;
        double dEtaIn = -999.;
        double dPhiIn = -999.;
        double ooEmooP = -999.;
        
        
        // New variables
        double cut = -999.; 
        double BDT = -999.;
        double corrpt = -999.;
        double FR = -999.;
        double CF = -999.;
        double passConversion = -999.;
        double passMuTightCharge = -999.;
        double passEleTightCharge = -999.;
        double passMissHit = -999.;
        double isMatchRightCharge = 1.; // lepton match gen right charge? default 1 (yes)
        
        // lepton gen
        double gen_pt = -999.;
        double gen_eta = -999.;
        double gen_phi = -999.;
        double gen_en = -999.;
        double gen_pdgId = -999.;
        double genMother_pt = -999.;
        double genMother_eta = -999.;
        double genMother_phi = -999.;
        double genMother_en = -999.;
        double genMother_pdgId = -999.;
        double genGrandMother_pt = -999.;
        double genGrandMother_eta = -999.;
        double genGrandMother_phi = -999.;
        double genGrandMother_en = -999.;
        double genGrandMother_pdgId = -999.;
        double gen_isPromptFinalState = -999.;
        double gen_isDirectPromptTauDecayProductFinalState = -999.;

        // Object Identification
        bool mu_isLoose_tthlep();
        bool mu_isfake_tthlep();
        bool mu_isTight_tthlep(bool isMedium);
        bool ele_isLoose_tthlep();
        bool ele_isfake_tthlep();
        bool ele_isTight_tthlep();
        void set_Wp_tthlep(bool isMedium);
        void set_Wp_tthlep(bool isMedium, int& numLoose, int& numfake, int& numtight);
        void cal_conept(bool isMedium);
        void cal_tight_property();
        double get_valX_valY_binContent(TH2F* h, double valX, double valY);

    private:
        bool ele_passCuts();
};
