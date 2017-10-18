#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"

#include "Lepton.cc"

using namespace std;
////
//   Declare constants
/////
const int nentries = 1000;//-1 is all entries  51 for first DiMuSR
const bool debug =false;
string synchro = "";
const double evt = 653077.;
const bool sync_top = false;
string sample = "mc";
bool isHIPSafe = true;
string data_path = "/home/binghuan/Work/RootTestFiles/TTHLep_2017/data/weights/";


//Variable handleling
void rSetBranchAddress(TTree* readingtree, string sample);
void wSetBranchAddress(TTree* newtree, string sample);
void wClearInitialization(string sample);
void rGetEntry(Long64_t tentry, string sample);
void Muon_sel(string sample);
void patElectron_sel(string sample);
void Event_sel();
void Lep_sel();


/////
//utils
/////
double deltaPhi(double phi1, double phi2);
double deltaEta(double eta1, double eta2);
double deltaR(double dphi, double deta);
bool compare_pt(const Lepton& LeptonA, const Lepton& LeptonB);

//Read MVA's
void set_wgtMVA();
double get_LeptonMVA(Lepton lep); 

// lepton mva
TMVA::Reader *mu_reader_;
TMVA::Reader *ele_reader_;
Float_t varpt;
Float_t vareta;
Float_t varneuRelIso;
Float_t varchRelIso;
Float_t varjetPtRel_in;
Float_t varjetPtRatio_in;
Float_t varjetBTagCSV_in;
Float_t varjetNDauCharged_in;
Float_t varsip3d;
Float_t varmvaId;
Float_t vardxy;
Float_t vardz;
Float_t varSegCompat;

//Charge Flip
string ChargeFlipName = data_path + "QF_data_el.root";
TFile* CFfile = new TFile(ChargeFlipName.c_str(), "read");
TH2F* hist_cf = (TH2F*) CFfile->Get("chargeMisId");

//Fake Rate
string FakeRateName = data_path + "FR_data_ttH_mva.root";
TFile* FRfile = new TFile(FakeRateName.c_str(),"read");
TH2F* hist_mu_fr = (TH2F*) FRfile->Get("FR_mva090_mu_data_comb");
TH2F* hist_el_fr = (TH2F*) FRfile->Get("FR_mva090_el_data_comb_NC");


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
vector<double>* rMuon_isGlobal; TBranch* b_rMuon_isGlobal =0;
vector<double>* rMuon_chi2; TBranch* b_rMuon_chi2 =0;
vector<double>* rMuon_chi2LocalPosition; TBranch* b_rMuon_chi2LocalPosition =0;
vector<double>* rMuon_trkKink; TBranch* b_rMuon_trkKink =0;
vector<double>* rMuon_validFraction; TBranch* b_rMuon_validFraction =0;
vector<double>* rMuon_segmentCompatibility; TBranch* b_rMuon_segmentCompatibility =0;
vector<double>* rMuon_jetptratio; TBranch* b_rMuon_jetptratio =0;
vector<double>* rMuon_jetpt; TBranch* b_rMuon_jetpt =0;
vector<double>* rMuon_jetcsv; TBranch* b_rMuon_jetcsv =0;
vector<double>* rMuon_lepjetchtrks; TBranch* b_rMuon_lepjetchtrks =0;
vector<double>* rMuon_miniIsoCh; TBranch* b_rMuon_miniIsoCh =0;
vector<double>* rMuon_miniIsoPUsub; TBranch* b_rMuon_miniIsoPUsub =0;
vector<double>* rMuon_ptrel; TBranch* b_rMuon_ptrel =0;
vector<double>* rMuon_pTErrOVpT_it; TBranch* b_rMuon_pTErrOVpT_it =0;
vector<double>* rMuon_px; TBranch* b_rMuon_px =0;
vector<double>* rMuon_py; TBranch* b_rMuon_py =0;
vector<double>* rMuon_pz; TBranch* b_rMuon_pz =0;
vector<double>* rMuon_jetdr; TBranch* b_rMuon_jetdr =0;
vector<double>* rMuon_gen_pt; TBranch* b_rMuon_gen_pt =0;
vector<double>* rMuon_gen_eta; TBranch* b_rMuon_gen_eta =0;
vector<double>* rMuon_gen_phi; TBranch* b_rMuon_gen_phi =0;
vector<double>* rMuon_gen_en; TBranch* b_rMuon_gen_en =0;
vector<double>* rMuon_gen_pdgId; TBranch* b_rMuon_gen_pdgId =0;
vector<double>* rMuon_genMother_pt; TBranch* b_rMuon_genMother_pt =0;
vector<double>* rMuon_genMother_eta; TBranch* b_rMuon_genMother_eta =0;
vector<double>* rMuon_genMother_phi; TBranch* b_rMuon_genMother_phi =0;
vector<double>* rMuon_genMother_en; TBranch* b_rMuon_genMother_en =0;
vector<double>* rMuon_genMother_pdgId; TBranch* b_rMuon_genMother_pdgId =0;
vector<double>* rMuon_genGrandMother_pt; TBranch* b_rMuon_genGrandMother_pt =0;
vector<double>* rMuon_genGrandMother_eta; TBranch* b_rMuon_genGrandMother_eta =0;
vector<double>* rMuon_genGrandMother_phi; TBranch* b_rMuon_genGrandMother_phi =0;
vector<double>* rMuon_genGrandMother_en; TBranch* b_rMuon_genGrandMother_en =0;
vector<double>* rMuon_genGrandMother_pdgId; TBranch* b_rMuon_genGrandMother_pdgId =0;
vector<double>* rMuon_gen_isPromptFinalState; TBranch* b_rMuon_gen_isPromptFinalState =0;
vector<double>* rMuon_gen_isDirectPromptTauDecayProductFinalState; TBranch* b_rMuon_gen_isDirectPromptTauDecayProductFinalState =0;
//Electron
vector<double>* rpatElectron_pt; TBranch* b_rpatElectron_pt =0;
vector<double>* rpatElectron_eta; TBranch* b_rpatElectron_eta =0;
vector<double>* rpatElectron_phi; TBranch* b_rpatElectron_phi =0;
vector<double>* rpatElectron_energy; TBranch* b_rpatElectron_energy =0;
vector<double>* rpatElectron_IP3Dsig; TBranch* b_rpatElectron_IP3Dsig =0;
vector<double>* rpatElectron_miniIsoRel; TBranch* b_rpatElectron_miniIsoRel =0;
vector<double>* rpatElectron_charge; TBranch* b_rpatElectron_charge =0;
vector<double>* rpatElectron_pdgId; TBranch* b_rpatElectron_pdgId =0;
vector<double>* rpatElectron_gsfTrack_dxy_pv; TBranch* b_rpatElectron_gsfTrack_dxy_pv =0;
vector<double>* rpatElectron_gsfTrack_dz_pv; TBranch* b_rpatElectron_gsfTrack_dz_pv =0;
vector<double>* rpatElectron_jetptratio; TBranch* b_rpatElectron_jetptratio =0;
vector<double>* rpatElectron_jetcsv; TBranch* b_rpatElectron_jetcsv =0;
vector<double>* rpatElectron_passConversionVeto; TBranch* b_rpatElectron_passConversionVeto =0;
vector<double>* rpatElectron_jetpt; TBranch* b_rpatElectron_jetpt =0;
vector<double>* rpatElectron_lepjetchtrks; TBranch* b_rpatElectron_lepjetchtrks =0;
vector<double>* rpatElectron_miniIsoCh; TBranch* b_rpatElectron_miniIsoCh =0;
vector<double>* rpatElectron_miniIsoPUsub; TBranch* b_rpatElectron_miniIsoPUsub =0;
vector<double>* rpatElectron_ptrel; TBranch* b_rpatElectron_ptrel =0;
vector<double>* rpatElectron_px; TBranch* b_rpatElectron_px =0;
vector<double>* rpatElectron_py; TBranch* b_rpatElectron_py =0;
vector<double>* rpatElectron_pz; TBranch* b_rpatElectron_pz =0;
vector<double>* rpatElectron_jetdr; TBranch* b_rpatElectron_jetdr =0;
vector<double>* rpatElectron_gen_pt; TBranch* b_rpatElectron_gen_pt =0;
vector<double>* rpatElectron_gen_eta; TBranch* b_rpatElectron_gen_eta =0;
vector<double>* rpatElectron_gen_phi; TBranch* b_rpatElectron_gen_phi =0;
vector<double>* rpatElectron_gen_en; TBranch* b_rpatElectron_gen_en =0;
vector<double>* rpatElectron_gen_pdgId; TBranch* b_rpatElectron_gen_pdgId =0;
vector<double>* rpatElectron_genMother_pt; TBranch* b_rpatElectron_genMother_pt =0;
vector<double>* rpatElectron_genMother_eta; TBranch* b_rpatElectron_genMother_eta =0;
vector<double>* rpatElectron_genMother_phi; TBranch* b_rpatElectron_genMother_phi =0;
vector<double>* rpatElectron_genMother_en; TBranch* b_rpatElectron_genMother_en =0;
vector<double>* rpatElectron_genMother_pdgId; TBranch* b_rpatElectron_genMother_pdgId =0;
vector<double>* rpatElectron_genGrandMother_pt; TBranch* b_rpatElectron_genGrandMother_pt =0;
vector<double>* rpatElectron_genGrandMother_eta; TBranch* b_rpatElectron_genGrandMother_eta =0;
vector<double>* rpatElectron_genGrandMother_phi; TBranch* b_rpatElectron_genGrandMother_phi =0;
vector<double>* rpatElectron_genGrandMother_en; TBranch* b_rpatElectron_genGrandMother_en =0;
vector<double>* rpatElectron_genGrandMother_pdgId; TBranch* b_rpatElectron_genGrandMother_pdgId =0;
vector<double>* rpatElectron_gen_isPromptFinalState; TBranch* b_rpatElectron_gen_isPromptFinalState =0;
vector<double>* rpatElectron_gen_isDirectPromptTauDecayProductFinalState; TBranch* b_rpatElectron_gen_isDirectPromptTauDecayProductFinalState =0;
vector<double>* rpatElectron_SCeta; TBranch* b_rpatElectron_SCeta =0;
vector<double>* rpatElectron_mvaValue_HZZ; TBranch* b_rpatElectron_mvaValue_HZZ =0;
vector<double>* rpatElectron_expectedMissingInnerHits; TBranch* b_rpatElectron_expectedMissingInnerHits =0;
vector<double>* rpatElectron_full5x5_sigmaIetaIeta; TBranch* b_rpatElectron_full5x5_sigmaIetaIeta =0;
vector<double>* rpatElectron_hOverE; TBranch* b_rpatElectron_hOverE =0;
vector<double>* rpatElectron_dEtaIn; TBranch* b_rpatElectron_dEtaIn =0;
vector<double>* rpatElectron_dPhiIn; TBranch* b_rpatElectron_dPhiIn =0;
vector<double>* rpatElectron_ooEmooP; TBranch* b_rpatElectron_ooEmooP =0;
vector<double>* rpatElectron_isGsfCtfScPixChargeConsistent; TBranch* b_rpatElectron_isGsfCtfScPixChargeConsistent =0;
vector<double>* rpatElectron_isGsfScPixChargeConsistent; TBranch* b_rpatElectron_isGsfScPixChargeConsistent =0;

//variables to be written

// Event level variables
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
double Muon_numLoose;
double Muon_numFake;
double Muon_numTight;
double patElectron_numLoose;
double patElectron_numFake;
double patElectron_numTight;


//Lepton
vector<Lepton>* leptons = new std::vector<Lepton>;

vector<double>* Lepton_pt = new std::vector<double>;
vector<double>* Lepton_eta = new std::vector<double>;
vector<double>* Lepton_phi = new std::vector<double>;
vector<double>* Lepton_energy = new std::vector<double>;
vector<double>* Lepton_dxy_pv = new std::vector<double>;
vector<double>* Lepton_dz_pv = new std::vector<double>;
vector<double>* Lepton_IP3Dsig = new std::vector<double>;
vector<double>* Lepton_miniIsoRel = new std::vector<double>;
vector<double>* Lepton_loose = new std::vector<double>;
vector<double>* Lepton_charge = new std::vector<double>;
vector<double>* Lepton_pdgId = new std::vector<double>;
vector<double>* Lepton_isGlobal = new std::vector<double>;
vector<double>* Lepton_chi2 = new std::vector<double>;
vector<double>* Lepton_chi2LocalPosition = new std::vector<double>;
vector<double>* Lepton_trkKink = new std::vector<double>;
vector<double>* Lepton_validFraction = new std::vector<double>;
vector<double>* Lepton_segmentCompatibility = new std::vector<double>;
vector<double>* Lepton_jetptratio = new std::vector<double>;
vector<double>* Lepton_jetpt = new std::vector<double>;
vector<double>* Lepton_jetcsv = new std::vector<double>;
vector<double>* Lepton_lepjetchtrks = new std::vector<double>;
vector<double>* Lepton_miniIsoCh = new std::vector<double>;
vector<double>* Lepton_miniIsoPUsub = new std::vector<double>;
vector<double>* Lepton_ptrel = new std::vector<double>;
vector<double>* Lepton_pTErrOVpT_it = new std::vector<double>;
vector<double>* Lepton_px = new std::vector<double>;
vector<double>* Lepton_py = new std::vector<double>;
vector<double>* Lepton_pz = new std::vector<double>;
vector<double>* Lepton_jetdr = new std::vector<double>;
vector<double>* Lepton_gen_pt = new std::vector<double>;
vector<double>* Lepton_gen_eta = new std::vector<double>;
vector<double>* Lepton_gen_phi = new std::vector<double>;
vector<double>* Lepton_gen_en = new std::vector<double>;
vector<double>* Lepton_gen_pdgId = new std::vector<double>;
vector<double>* Lepton_genMother_pt = new std::vector<double>;
vector<double>* Lepton_genMother_eta = new std::vector<double>;
vector<double>* Lepton_genMother_phi = new std::vector<double>;
vector<double>* Lepton_genMother_en = new std::vector<double>;
vector<double>* Lepton_genMother_pdgId = new std::vector<double>;
vector<double>* Lepton_genGrandMother_pt = new std::vector<double>;
vector<double>* Lepton_genGrandMother_eta = new std::vector<double>;
vector<double>* Lepton_genGrandMother_phi = new std::vector<double>;
vector<double>* Lepton_genGrandMother_en = new std::vector<double>;
vector<double>* Lepton_genGrandMother_pdgId = new std::vector<double>;
vector<double>* Lepton_gen_isPromptFinalState = new std::vector<double>;
vector<double>* Lepton_gen_isDirectPromptTauDecayProductFinalState = new std::vector<double>;
vector<double>* Lepton_SCeta = new std::vector<double>;
vector<double>* Lepton_mvaValue_HZZ = new std::vector<double>;
vector<double>* Lepton_expectedMissingInnerHits = new std::vector<double>;
vector<double>* Lepton_full5x5_sigmaIetaIeta = new std::vector<double>;
vector<double>* Lepton_hOverE = new std::vector<double>;
vector<double>* Lepton_dEtaIn = new std::vector<double>;
vector<double>* Lepton_dPhiIn = new std::vector<double>;
vector<double>* Lepton_ooEmooP = new std::vector<double>;

//new variables
vector<double>* Lepton_cut = new std::vector<double>;
vector<double>* Lepton_BDT = new std::vector<double>;
vector<double>* Lepton_corrpt = new std::vector<double>;
vector<double>* Lepton_FR = new std::vector<double>;
vector<double>* Lepton_CF = new std::vector<double>;
vector<double>* Lepton_passConversion = new std::vector<double>;
vector<double>* Lepton_passMuTightCharge = new std::vector<double>;
vector<double>* Lepton_passEleTightCharge = new std::vector<double>;
vector<double>* Lepton_passMissHit = new std::vector<double>;
vector<double>* Lepton_isMatchRightCharge = new std::vector<double>;
