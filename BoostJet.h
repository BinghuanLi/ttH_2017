#include "TLorentzVector.h"

class BoostJet {
    // this class define leptons
    public:
        // constructor
        BoostJet() = default;
        BoostJet(const BoostJet&) = default;
        BoostJet& operator=(const BoostJet&);
        ~BoostJet() = default;

        // member
        //readable variables 
        // reco variables
        double pt = -999.;
        double eta = -999.;
        double phi = -999.;
        double energy = -999.;
        double Uncorr_pt = -999.;
        double JesSF = -999.;
        double JesSFup = -999.;
        double JesSFdown = -999.;
        double JerSF = -999.;
        double JerSFup = -999.;
        double JerSFdown = -999.;
        double tau1 = 999.;
        double tau2 = 999.;
        double tau3 = 999.;
        double softdrop_mass = -999.;
        double pruned_mass = -999.;

        //BoostJet id
        double pfCombinedInclusiveSecondaryVertexV2BJetTags = -999.;
        double pfCombinedMVAV2BJetTags = -999.;
    
        // new variables
        double topCut = -999.;
        double wCut = -999.;
        double tau32 = 999.;
        double tau21 = 999.;
       
        //Object Identification
        void set_Wp_Top(int& numSoftTop, int& numLooseTop, int& numMediumTop, int& numTightTop);
        void set_Wp_W(int& numLooseW, int& numTightW);
};
