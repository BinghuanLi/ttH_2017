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
        //lepton id
        double loose;
        double pdgId;

        // new variables
        double cut; 
        
        //Object Identification
        bool mu_isLoose_tthlep();
        void set_Wp_tthlep();

};
