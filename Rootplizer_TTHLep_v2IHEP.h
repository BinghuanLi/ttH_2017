#include "TFile.h"
#include "TTree.h"
#include <iostream>

#include "Lepton.cc"

using namespace std;
////
//   Declare constants
/////
const int nentries = 10;//-1 is all entries  51 for first DiMuSR
const bool debug =false;
string synchro = "";
const double evt = 653077.;
const bool sync_top = false;
string sample = "mc";
bool isHIPSafe = true;

//Variable handleling
void rSetBranchAddress(TTree* readingtree, string sample);
void wSetBranchAddress(TTree* newtree, string sample);
void wClearInitialization(string sample);
void rGetEntry(Long64_t tentry, string sample);
void Muon_sel(string sample);
void Event_sel();
void Lep_sel();


//variables to be read
//Event
ULong64_t rEVENT_event; TBranch* b_rEVENT_event =0;
double rEVENT_genWeight; TBranch* b_rEVENT_genWeight =0;
int rHiggsDecay; TBranch* b_rHiggsDecay =0;
double rMet_type1PF_pt; TBranch* b_rMet_type1PF_pt =0;
double rMet_type1PF_px; TBranch* b_rMet_type1PF_px =0;
double rMet_type1PF_py; TBranch* b_rMet_type1PF_py =0;
double rMet_type1PF_pz; TBranch* b_rMet_type1PF_pz =0;
double rMet_type1PF_phi; TBranch* b_rMet_type1PF_phi =0;
double rMet_type1PF_sumEt; TBranch* b_rMet_type1PF_sumEt =0;
double rMet_type1PF_shiftedPtUp; TBranch* b_rMet_type1PF_shiftedPtUp =0;
double rMet_type1PF_shiftedPtDown; TBranch* b_rMet_type1PF_shiftedPtDown =0;
//Muon
vector<double>* rMuon_pt; TBranch* b_rMuon_pt =0;
vector<double>* rMuon_eta; TBranch* b_rMuon_eta =0;
vector<double>* rMuon_phi; TBranch* b_rMuon_phi =0;
vector<double>* rMuon_energy; TBranch* b_rMuon_energy =0;
vector<double>* rMuon_dxy_pv; TBranch* b_rMuon_dxy_pv =0;
vector<double>* rMuon_dz_pv; TBranch* b_rMuon_dz_pv =0;
vector<double>* rMuon_IP3Dsig_it; TBranch* b_rMuon_IP3Dsig_it =0;
vector<double>* rMuon_loose; TBranch* b_rMuon_loose =0;
vector<double>* rMuon_miniIsoRel; TBranch* b_rMuon_miniIsoRel =0;
vector<double>* rMuon_charge; TBranch* b_rMuon_charge =0;
vector<double>* rMuon_pdgId; TBranch* b_rMuon_pdgId =0;

//variables to be written
double EVENT_event;
double EVENT_genWeight;
double HiggsDecay;
double Met_type1PF_pt;
double Met_type1PF_px;
double Met_type1PF_py;
double Met_type1PF_pz;
double Met_type1PF_phi;
double Met_type1PF_sumEt;
double Met_type1PF_shiftedPtUp;
double Met_type1PF_shiftedPtDown;

//Lepton
vector<Lepton>* leptons = new std::vector<Lepton>;

vector<double>* Lepton_pt = new std::vector<double>;
vector<double>* Lepton_eta = new std::vector<double>;
vector<double>* Lepton_phi = new std::vector<double>;
vector<double>* Lepton_energy = new std::vector<double>;
vector<double>* Lepton_dxy_pv = new std::vector<double>;
vector<double>* Lepton_dz_pv = new std::vector<double>;
vector<double>* Lepton_IP3Dsig_it = new std::vector<double>;
vector<double>* Lepton_miniIsoRel = new std::vector<double>;
vector<double>* Lepton_loose = new std::vector<double>;
vector<double>* Lepton_charge = new std::vector<double>;
vector<double>* Lepton_pdgId = new std::vector<double>;

//new variables
vector<double>* Lepton_cut = new std::vector<double>;
