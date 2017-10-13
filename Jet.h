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
        double pt;
        double eta;
        double phi;
        double energy;
        double Uncorr_pt;
        double JesSF;
        //Jet id
        double pfCombinedInclusiveSecondaryVertexV2BJetTags;

        // new variables
       
        
        //Object Identification

};
