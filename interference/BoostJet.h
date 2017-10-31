#include "TLorentzVector.h"

////
// utils
////
//deltaR calculation 
//defined in ../source/utils_functions
       
extern double deltaPhi(double phi1, double phi2);
extern double deltaEta(double eta1, double eta2);
extern double deltaR(double dphi, double deta);

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
        /////
        // gen match
        ////
        // gen match W
        double matchW_pt = -999.;
        double matchW_eta = -999.;
        double matchW_phi = -999.;
        double matchW_energy = -999.;
        double matchW_mass = -999.;
        double matchW_mother_pdgId = -999.;

        //Object Identification
        void set_Wp_Top(int& numSoftTop, int& numLooseTop, int& numMediumTop, int& numTightTop);
        void set_Wp_W(
            int& numSoftW, int& numLooseW, int& numMediumW, int& numTightW,
            int& numCleanSoftW, int& numCleanLooseW, int& numCleanMediumW, int& numCleanTightW,
            vector<double>* jet_pt, vector<double>* jet_phi
            );
        void set_Wp_W(int& numLooseW, int& numTightW);
        void match_genW(vector<double>* genW_pt, vector<double>* genW_eta, vector<double>* phi, 
            vector<double>* genW_energy, vector<double>* genW_mass, vector<double>* genW_motherId);
};
