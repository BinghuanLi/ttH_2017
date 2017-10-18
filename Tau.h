#include "TLorentzVector.h"

class Tau {
    // this class define leptons
    public:
        // constructor
        Tau() = default;
        Tau(const Tau&) = default;
        Tau& operator=(const Tau&);
        ~Tau() = default;

        // member
        //readable variables 
        // reco variables
        double pt = -999.;
        double eta = -999.;
        double phi = -999.;
        double energy = -999.;
        double charge = -999.;
        double packedLeadTauCand_dz = -999.;
        double packedLeadTauCand_dxy = -999.;
        
        //Tau id
        double byLooseIsolationMVArun2v1DBdR03oldDMwLT = -999.;
        double decayModeFinding = -999.;
        double byMediumIsolationMVArun2v1DBdR03oldDMwLT = -999.;

        // new variables
        double cut = -999;
        
        //Object Identification
        bool tau_isLoose_tthlep();
        bool tau_isMedium_tthlep();
        void set_Wp_tthlep(int& numLoose, int& numMedium);

};
