#include "TLorentzVector.h"

class Jet {
    // this class define leptons
    public:
        // constructor
        Jet() = default;
        Jet(const Jet&) = default;
        Jet& operator=(const Jet&);
        ~Jet() = default;

        // member
        //readable variables 
        // reco variables
        double pt = -999.;
        double eta = -999.;
        double phi = -999.;
        double energy = -999.;
        double Uncorr_pt = -999.;
        double neutralHadEnergyFraction = -999.;
        double neutralEmEnergyFraction = -999.;
        double chargedMultiplicity = -999.;
        double numberOfConstituents = -999.;
        double chargedHadronEnergyFraction = -999.;
        double chargedEmEnergyFraction = -999.;
        double JesSF = -999.;
        double JesSFup = -999.;
        double JesSFdown = -999.;
        double JerSF = -999.;
        double JerSFup = -999.;
        double JerSFdown = -999.;
        double axis2 = -999.;
        double ptD = -999.;
        double mult = -999.;
        double btag_sf = -999.;
        double btag_jesup = -999.;
        double btag_jesdown = -999.;
        double btag_hfup = -999.;
        double btag_hfdown = -999.;
        double btag_hfstat1up = -999.;
        double btag_hfstat1down = -999.;
        double btag_hfstat2up = -999.;
        double btag_hfstat2down = -999.;
        double btag_lfup = -999.;
        double btag_lfdown = -999.;
        double btag_lfstat1up = -999.;
        double btag_lfstat1down = -999.;
        double btag_lfstat2up = -999.;
        double btag_lfstat2down = -999.;
        double btag_cerr1up = -999.;
        double btag_cerr1down = -999.;
        double btag_cerr2up = -999.;
        double btag_cerr2down = -999.;

        //Jet id
        double pfCombinedInclusiveSecondaryVertexV2BJetTags = -999.;
        double pfCombinedMVAV2BJetTags = -999.;
        double qg = -999.;

        //gen information
        double partonFlavour = -999.;
        double hadronFlavour = -999.;
        double genpt = -999.;
        double geneta = -999.;
        double genphi = -999.;
        double genenergy = -999.;
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
        
        // new variables
        double cut = -999;
        double bcut = -999;
       
        //Object Identification
        bool jet_isLoose();
        void set_Wp_jets(int& numLoose);
        void set_Wp_bdisc(int& numbLoose, int& numbMedium, int& numbTight);
};
