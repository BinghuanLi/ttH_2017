#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"

#include "Lepton.cc"
#include "Tau.cc"
#include "Jet.cc"
#include "BoostJet.cc"

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
void Tau_sel();
void Jet_sel(string sample); 
void BoostedJet_sel();
void Event_sel( string OutputName);
void Lep_sel();
void GenParticle_sel();


/////
//utils
/////
double deltaPhi(double phi1, double phi2);
double deltaEta(double eta1, double eta2);
double deltaR(double dphi, double deta);
bool byPt(const Lepton& LeptonA, const Lepton& LeptonB);
bool byConept(const Lepton& LeptonA, const Lepton& LeptonB);

////
// event weights
///
double get_wgtlumi(string FileName);

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
double rPUWeight; TBranch* b_rPUWeight =0;
double rGen_type1PF_Met; TBranch* b_rGen_type1PF_Met =0;
double rGen_type1PF_Metpx; TBranch* b_rGen_type1PF_Metpx =0;
double rGen_type1PF_Metpy; TBranch* b_rGen_type1PF_Metpy =0;
double rGen_type1PF_Metpz; TBranch* b_rGen_type1PF_Metpz =0;
double rGen_type1PF_Meteta; TBranch* b_rGen_type1PF_Meteta =0;
double rGen_type1PF_Metphi; TBranch* b_rGen_type1PF_Metphi =0;
double rGen_type1PF_Meten; TBranch* b_rGen_type1PF_Meten =0;
double rMet_type1PF_pt; TBranch* b_rMet_type1PF_pt =0;
double rMet_type1PF_px; TBranch* b_rMet_type1PF_px =0;
double rMet_type1PF_py; TBranch* b_rMet_type1PF_py =0;
double rMet_type1PF_pz; TBranch* b_rMet_type1PF_pz =0;
double rMet_type1PF_phi; TBranch* b_rMet_type1PF_phi =0;
double rMet_type1PF_sumEt; TBranch* b_rMet_type1PF_sumEt =0;
double rMet_type1PF_shiftedPtUp; TBranch* b_rMet_type1PF_shiftedPtUp =0;
double rMet_type1PF_shiftedPtDown; TBranch* b_rMet_type1PF_shiftedPtDown =0;
int rHLT_DiMu9_Ele9_CaloIdL_TrackIdL; TBranch* b_rHLT_DiMu9_Ele9_CaloIdL_TrackIdL =0;
int rHLT_Mu8_DiEle12_CaloIdL_TrackIdL; TBranch* b_rHLT_Mu8_DiEle12_CaloIdL_TrackIdL =0;
int rHLT_TripleMu_12_10_5; TBranch* b_rHLT_TripleMu_12_10_5 =0;
int rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL; TBranch* b_rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL =0;
int rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL; TBranch* b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL =0;
int rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ; TBranch* b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ =0;
int rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL; TBranch* b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL =0;
int rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ; TBranch* b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ =0;
int rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ; TBranch* b_rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ =0;
int rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ; TBranch* b_rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ =0;
int rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ; TBranch* b_rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ =0;
int rHLT_IsoMu22; TBranch* b_rHLT_IsoMu22 =0;
int rHLT_IsoTkMu22; TBranch* b_rHLT_IsoTkMu22 =0;
int rHLT_IsoMu22_eta2p1; TBranch* b_rHLT_IsoMu22_eta2p1 =0;
int rHLT_IsoTkMu22_eta2p1; TBranch* b_rHLT_IsoTkMu22_eta2p1 =0;
int rHLT_IsoMu24; TBranch* b_rHLT_IsoMu24 =0;
int rHLT_IsoTkMu24; TBranch* b_rHLT_IsoTkMu24 =0;
int rHLT_Ele27_WPTight_Gsf; TBranch* b_rHLT_Ele27_WPTight_Gsf =0;
int rHLT_Ele25_eta2p1_WPTight_Gsf; TBranch* b_rHLT_Ele25_eta2p1_WPTight_Gsf =0;
int rHLT_Ele27_eta2p1_WPLoose_Gsf; TBranch* b_rHLT_Ele27_eta2p1_WPLoose_Gsf =0;

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

//Tau
vector<double>* rTau_pt; TBranch* b_rTau_pt =0;
vector<double>* rTau_eta; TBranch* b_rTau_eta =0;
vector<double>* rTau_phi; TBranch* b_rTau_phi =0;
vector<double>* rTau_energy; TBranch* b_rTau_energy =0;
vector<double>* rTau_charge; TBranch* b_rTau_charge =0;
vector<double>* rTau_packedLeadTauCand_dz; TBranch* b_rTau_packedLeadTauCand_dz =0;
vector<double>* rTau_packedLeadTauCand_dxy; TBranch* b_rTau_packedLeadTauCand_dxy =0;
vector<double>* rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT; TBranch* b_rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT =0;
vector<double>* rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT; TBranch* b_rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT =0;
vector<double>* rTau_decayModeFinding; TBranch* b_rTau_decayModeFinding =0;

//Jet
vector<double>* rJet_pt; TBranch* b_rJet_pt =0;
vector<double>* rJet_eta; TBranch* b_rJet_eta =0;
vector<double>* rJet_phi; TBranch* b_rJet_phi =0;
vector<double>* rJet_energy; TBranch* b_rJet_energy =0;
vector<double>* rJet_genMother_pt; TBranch* b_rJet_genMother_pt =0;
vector<double>* rJet_genMother_eta; TBranch* b_rJet_genMother_eta =0;
vector<double>* rJet_genMother_phi; TBranch* b_rJet_genMother_phi =0;
vector<double>* rJet_genMother_en; TBranch* b_rJet_genMother_en =0;
vector<double>* rJet_genMother_pdgId; TBranch* b_rJet_genMother_pdgId =0;
vector<double>* rJet_genGrandMother_pt; TBranch* b_rJet_genGrandMother_pt =0;
vector<double>* rJet_genGrandMother_eta; TBranch* b_rJet_genGrandMother_eta =0;
vector<double>* rJet_genGrandMother_phi; TBranch* b_rJet_genGrandMother_phi =0;
vector<double>* rJet_genGrandMother_en; TBranch* b_rJet_genGrandMother_en =0;
vector<double>* rJet_genGrandMother_pdgId; TBranch* b_rJet_genGrandMother_pdgId =0;
vector<double>* rJet_Uncorr_pt; TBranch* b_rJet_Uncorr_pt =0;
vector<double>* rJet_neutralHadEnergyFraction; TBranch* b_rJet_neutralHadEnergyFraction =0;
vector<double>* rJet_neutralEmEnergyFraction; TBranch* b_rJet_neutralEmEnergyFraction =0;
vector<double>* rJet_chargedMultiplicity; TBranch* b_rJet_chargedMultiplicity =0;
vector<double>* rJet_numberOfConstituents; TBranch* b_rJet_numberOfConstituents =0;
vector<double>* rJet_chargedHadronEnergyFraction; TBranch* b_rJet_chargedHadronEnergyFraction =0;
vector<double>* rJet_chargedEmEnergyFraction; TBranch* b_rJet_chargedEmEnergyFraction =0;
vector<double>* rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags; TBranch* b_rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags =0;
vector<double>* rJet_pfCombinedMVAV2BJetTags; TBranch* b_rJet_pfCombinedMVAV2BJetTags =0;
vector<double>* rJet_JesSF; TBranch* b_rJet_JesSF =0;
vector<double>* rJet_JesSFup; TBranch* b_rJet_JesSFup =0;
vector<double>* rJet_JesSFdown; TBranch* b_rJet_JesSFdown =0;
vector<double>* rJet_JerSF; TBranch* b_rJet_JerSF =0;
vector<double>* rJet_JerSFup; TBranch* b_rJet_JerSFup =0;
vector<double>* rJet_JerSFdown; TBranch* b_rJet_JerSFdown =0;
vector<double>* rJet_qg; TBranch* b_rJet_qg =0;
vector<double>* rJet_axis2; TBranch* b_rJet_axis2 =0;
vector<double>* rJet_ptD; TBranch* b_rJet_ptD =0;
vector<double>* rJet_mult; TBranch* b_rJet_mult =0;
vector<double>* rJet_partonFlavour; TBranch* b_rJet_partonFlavour =0;
vector<double>* rJet_hadronFlavour; TBranch* b_rJet_hadronFlavour =0;
vector<double>* rJet_genpt; TBranch* b_rJet_genpt =0;
vector<double>* rJet_geneta; TBranch* b_rJet_geneta =0;
vector<double>* rJet_genphi; TBranch* b_rJet_genphi =0;
vector<double>* rJet_genenergy; TBranch* b_rJet_genenergy =0;
vector<double>* rJet_btag_sf; TBranch* b_rJet_btag_sf =0;
vector<double>* rJet_btag_jesup; TBranch* b_rJet_btag_jesup =0;
vector<double>* rJet_btag_jesdown; TBranch* b_rJet_btag_jesdown =0;
vector<double>* rJet_btag_hfup; TBranch* b_rJet_btag_hfup =0;
vector<double>* rJet_btag_hfdown; TBranch* b_rJet_btag_hfdown =0;
vector<double>* rJet_btag_hfstat1up; TBranch* b_rJet_btag_hfstat1up =0;
vector<double>* rJet_btag_hfstat1down; TBranch* b_rJet_btag_hfstat1down =0;
vector<double>* rJet_btag_hfstat2up; TBranch* b_rJet_btag_hfstat2up =0;
vector<double>* rJet_btag_hfstat2down; TBranch* b_rJet_btag_hfstat2down =0;
vector<double>* rJet_btag_lfup; TBranch* b_rJet_btag_lfup =0;
vector<double>* rJet_btag_lfdown; TBranch* b_rJet_btag_lfdown =0;
vector<double>* rJet_btag_lfstat1up; TBranch* b_rJet_btag_lfstat1up =0;
vector<double>* rJet_btag_lfstat1down; TBranch* b_rJet_btag_lfstat1down =0;
vector<double>* rJet_btag_lfstat2up; TBranch* b_rJet_btag_lfstat2up =0;
vector<double>* rJet_btag_lfstat2down; TBranch* b_rJet_btag_lfstat2down =0;
vector<double>* rJet_btag_cerr1up; TBranch* b_rJet_btag_cerr1up =0;
vector<double>* rJet_btag_cerr1down; TBranch* b_rJet_btag_cerr1down =0;
vector<double>* rJet_btag_cerr2up; TBranch* b_rJet_btag_cerr2up =0;
vector<double>* rJet_btag_cerr2down; TBranch* b_rJet_btag_cerr2down =0;

//BoostJet
vector<double>* rBoostedJet_pt; TBranch* b_rBoostedJet_pt =0;
vector<double>* rBoostedJet_eta; TBranch* b_rBoostedJet_eta =0;
vector<double>* rBoostedJet_phi; TBranch* b_rBoostedJet_phi =0;
vector<double>* rBoostedJet_energy; TBranch* b_rBoostedJet_energy =0;
vector<double>* rBoostedJet_Uncorr_pt; TBranch* b_rBoostedJet_Uncorr_pt =0;
vector<double>* rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags; TBranch* b_rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags =0;
vector<double>* rBoostedJet_pfCombinedMVAV2BJetTags; TBranch* b_rBoostedJet_pfCombinedMVAV2BJetTags =0;
vector<double>* rBoostedJet_JesSF; TBranch* b_rBoostedJet_JesSF =0;
vector<double>* rBoostedJet_JesSFup; TBranch* b_rBoostedJet_JesSFup =0;
vector<double>* rBoostedJet_JesSFdown; TBranch* b_rBoostedJet_JesSFdown =0;
vector<double>* rBoostedJet_JerSF; TBranch* b_rBoostedJet_JerSF =0;
vector<double>* rBoostedJet_JerSFup; TBranch* b_rBoostedJet_JerSFup =0;
vector<double>* rBoostedJet_JerSFdown; TBranch* b_rBoostedJet_JerSFdown =0;
vector<double>* rBoostedJet_tau1; TBranch* b_rBoostedJet_tau1 =0;
vector<double>* rBoostedJet_tau2; TBranch* b_rBoostedJet_tau2 =0;
vector<double>* rBoostedJet_tau3; TBranch* b_rBoostedJet_tau3 =0;
vector<double>* rBoostedJet_softdrop_mass; TBranch* b_rBoostedJet_softdrop_mass =0;
vector<double>* rBoostedJet_pruned_mass; TBranch* b_rBoostedJet_pruned_mass =0;

//Gen
vector<double>* rGen_pdg_id; TBranch* b_rGen_pdg_id =0;
vector<double>* rGen_pt; TBranch* b_rGen_pt =0;
vector<double>* rGen_eta; TBranch* b_rGen_eta =0;
vector<double>* rGen_phi; TBranch* b_rGen_phi =0;
vector<double>* rGen_energy; TBranch* b_rGen_energy =0;
vector<double>* rGen_motherpdg_id; TBranch* b_rGen_motherpdg_id =0;
vector<double>* rGen_BmotherIndex; TBranch* b_rGen_BmotherIndex =0;
vector<double>* rGen_numMother; TBranch* b_rGen_numMother =0;
vector<double>* rGen_status; TBranch* b_rGen_status =0;

//variables to be written

// Event level variables
double lumi_wgt;
double EVENT_event;
double EVENT_genWeight;
double HiggsDecay;
double PUWeight;
double Gen_type1PF_Met;
double Gen_type1PF_Metpx;
double Gen_type1PF_Metpy;
double Gen_type1PF_Metpz;
double Gen_type1PF_Meteta;
double Gen_type1PF_Metphi;
double Gen_type1PF_Meten;
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
double Tau_numLoose;
double Tau_numMedium;
double Jet_numLoose;
double Jet_numbLoose;
double Jet_numbMedium;
double Jet_numbTight;
double BWeight;
double BWeightLFup;
double BWeightLFdown;
double BWeightHFup;
double BWeightHFdown;
double BWeightJESup;
double BWeightJESdown;
double BWeightLFStats1up;
double BWeightLFStats1down;
double BWeightLFStats2up;
double BWeightLFStats2down;
double BWeightHFStats1up;
double BWeightHFStats1down;
double BWeightHFStats2up;
double BWeightHFStats2down;
double BWeightCErr1up;
double BWeightCErr1down;
double BWeightCErr2up;
double BWeightCErr2down;
double Top_numSoft;
double Top_numLoose;
double Top_numMedium;
double Top_numTight;
double W_numLoose;
double W_numTight;
double TTHLep_2Mu;
double TTHLep_2Ele;
double TTHLep_MuEle;
double TTHLep_3L4L;
double metLD;
double mhtT_met;
double mht_met;
double mhtT;
double mht;


//Lepton
vector<Lepton>* leptons = new std::vector<Lepton>;

vector<double>* Lep_pt = new std::vector<double>;
vector<double>* Lep_eta = new std::vector<double>;
vector<double>* Lep_phi = new std::vector<double>;
vector<double>* Lep_energy = new std::vector<double>;
vector<double>* Lep_dxy_pv = new std::vector<double>;
vector<double>* Lep_dz_pv = new std::vector<double>;
vector<double>* Lep_IP3Dsig = new std::vector<double>;
vector<double>* Lep_miniIsoRel = new std::vector<double>;
vector<double>* Lep_loose = new std::vector<double>;
vector<double>* Lep_charge = new std::vector<double>;
vector<double>* Lep_pdgId = new std::vector<double>;
vector<double>* Lep_isGlobal = new std::vector<double>;
vector<double>* Lep_chi2 = new std::vector<double>;
vector<double>* Lep_chi2LocalPosition = new std::vector<double>;
vector<double>* Lep_trkKink = new std::vector<double>;
vector<double>* Lep_validFraction = new std::vector<double>;
vector<double>* Lep_segmentCompatibility = new std::vector<double>;
vector<double>* Lep_jetptratio = new std::vector<double>;
vector<double>* Lep_jetpt = new std::vector<double>;
vector<double>* Lep_jetcsv = new std::vector<double>;
vector<double>* Lep_lepjetchtrks = new std::vector<double>;
vector<double>* Lep_miniIsoCh = new std::vector<double>;
vector<double>* Lep_miniIsoPUsub = new std::vector<double>;
vector<double>* Lep_ptrel = new std::vector<double>;
vector<double>* Lep_pTErrOVpT_it = new std::vector<double>;
vector<double>* Lep_px = new std::vector<double>;
vector<double>* Lep_py = new std::vector<double>;
vector<double>* Lep_pz = new std::vector<double>;
vector<double>* Lep_jetdr = new std::vector<double>;
vector<double>* Lep_gen_pt = new std::vector<double>;
vector<double>* Lep_gen_eta = new std::vector<double>;
vector<double>* Lep_gen_phi = new std::vector<double>;
vector<double>* Lep_gen_en = new std::vector<double>;
vector<double>* Lep_gen_pdgId = new std::vector<double>;
vector<double>* Lep_genMother_pt = new std::vector<double>;
vector<double>* Lep_genMother_eta = new std::vector<double>;
vector<double>* Lep_genMother_phi = new std::vector<double>;
vector<double>* Lep_genMother_en = new std::vector<double>;
vector<double>* Lep_genMother_pdgId = new std::vector<double>;
vector<double>* Lep_genGrandMother_pt = new std::vector<double>;
vector<double>* Lep_genGrandMother_eta = new std::vector<double>;
vector<double>* Lep_genGrandMother_phi = new std::vector<double>;
vector<double>* Lep_genGrandMother_en = new std::vector<double>;
vector<double>* Lep_genGrandMother_pdgId = new std::vector<double>;
vector<double>* Lep_gen_isPromptFinalState = new std::vector<double>;
vector<double>* Lep_gen_isDirectPromptTauDecayProductFinalState = new std::vector<double>;
vector<double>* Lep_SCeta = new std::vector<double>;
vector<double>* Lep_mvaValue_HZZ = new std::vector<double>;
vector<double>* Lep_expectedMissingInnerHits = new std::vector<double>;
vector<double>* Lep_full5x5_sigmaIetaIeta = new std::vector<double>;
vector<double>* Lep_hOverE = new std::vector<double>;
vector<double>* Lep_dEtaIn = new std::vector<double>;
vector<double>* Lep_dPhiIn = new std::vector<double>;
vector<double>* Lep_ooEmooP = new std::vector<double>;

//new variables
vector<double>* Lep_cut = new std::vector<double>;
vector<double>* Lep_BDT = new std::vector<double>;
vector<double>* Lep_corrpt = new std::vector<double>;
vector<double>* Lep_FR = new std::vector<double>;
vector<double>* Lep_CF = new std::vector<double>;
vector<double>* Lep_passConversion = new std::vector<double>;
vector<double>* Lep_passMuTightCharge = new std::vector<double>;
vector<double>* Lep_passEleTightCharge = new std::vector<double>;
vector<double>* Lep_passMissHit = new std::vector<double>;
vector<double>* Lep_isMatchRightCharge = new std::vector<double>;

//Fakeable Leptons
vector<double>* FakeLep_pt = new std::vector<double>;
vector<double>* FakeLep_eta = new std::vector<double>;
vector<double>* FakeLep_phi = new std::vector<double>;
vector<double>* FakeLep_energy = new std::vector<double>;
vector<double>* FakeLep_dxy_pv = new std::vector<double>;
vector<double>* FakeLep_dz_pv = new std::vector<double>;
vector<double>* FakeLep_IP3Dsig = new std::vector<double>;
vector<double>* FakeLep_miniIsoRel = new std::vector<double>;
vector<double>* FakeLep_loose = new std::vector<double>;
vector<double>* FakeLep_charge = new std::vector<double>;
vector<double>* FakeLep_pdgId = new std::vector<double>;
vector<double>* FakeLep_isGlobal = new std::vector<double>;
vector<double>* FakeLep_chi2 = new std::vector<double>;
vector<double>* FakeLep_chi2LocalPosition = new std::vector<double>;
vector<double>* FakeLep_trkKink = new std::vector<double>;
vector<double>* FakeLep_validFraction = new std::vector<double>;
vector<double>* FakeLep_segmentCompatibility = new std::vector<double>;
vector<double>* FakeLep_jetptratio = new std::vector<double>;
vector<double>* FakeLep_jetpt = new std::vector<double>;
vector<double>* FakeLep_jetcsv = new std::vector<double>;
vector<double>* FakeLep_lepjetchtrks = new std::vector<double>;
vector<double>* FakeLep_miniIsoCh = new std::vector<double>;
vector<double>* FakeLep_miniIsoPUsub = new std::vector<double>;
vector<double>* FakeLep_ptrel = new std::vector<double>;
vector<double>* FakeLep_pTErrOVpT_it = new std::vector<double>;
vector<double>* FakeLep_px = new std::vector<double>;
vector<double>* FakeLep_py = new std::vector<double>;
vector<double>* FakeLep_pz = new std::vector<double>;
vector<double>* FakeLep_jetdr = new std::vector<double>;
vector<double>* FakeLep_gen_pt = new std::vector<double>;
vector<double>* FakeLep_gen_eta = new std::vector<double>;
vector<double>* FakeLep_gen_phi = new std::vector<double>;
vector<double>* FakeLep_gen_en = new std::vector<double>;
vector<double>* FakeLep_gen_pdgId = new std::vector<double>;
vector<double>* FakeLep_genMother_pt = new std::vector<double>;
vector<double>* FakeLep_genMother_eta = new std::vector<double>;
vector<double>* FakeLep_genMother_phi = new std::vector<double>;
vector<double>* FakeLep_genMother_en = new std::vector<double>;
vector<double>* FakeLep_genMother_pdgId = new std::vector<double>;
vector<double>* FakeLep_genGrandMother_pt = new std::vector<double>;
vector<double>* FakeLep_genGrandMother_eta = new std::vector<double>;
vector<double>* FakeLep_genGrandMother_phi = new std::vector<double>;
vector<double>* FakeLep_genGrandMother_en = new std::vector<double>;
vector<double>* FakeLep_genGrandMother_pdgId = new std::vector<double>;
vector<double>* FakeLep_gen_isPromptFinalState = new std::vector<double>;
vector<double>* FakeLep_gen_isDirectPromptTauDecayProductFinalState = new std::vector<double>;
vector<double>* FakeLep_SCeta = new std::vector<double>;
vector<double>* FakeLep_mvaValue_HZZ = new std::vector<double>;
vector<double>* FakeLep_expectedMissingInnerHits = new std::vector<double>;
vector<double>* FakeLep_full5x5_sigmaIetaIeta = new std::vector<double>;
vector<double>* FakeLep_hOverE = new std::vector<double>;
vector<double>* FakeLep_dEtaIn = new std::vector<double>;
vector<double>* FakeLep_dPhiIn = new std::vector<double>;
vector<double>* FakeLep_ooEmooP = new std::vector<double>;

//new variables
vector<double>* FakeLep_cut = new std::vector<double>;
vector<double>* FakeLep_BDT = new std::vector<double>;
vector<double>* FakeLep_corrpt = new std::vector<double>;
vector<double>* FakeLep_FR = new std::vector<double>;
vector<double>* FakeLep_CF = new std::vector<double>;
vector<double>* FakeLep_passConversion = new std::vector<double>;
vector<double>* FakeLep_passMuTightCharge = new std::vector<double>;
vector<double>* FakeLep_passEleTightCharge = new std::vector<double>;
vector<double>* FakeLep_passMissHit = new std::vector<double>;
vector<double>* FakeLep_isMatchRightCharge = new std::vector<double>;

// Tau
vector<Tau>* taus = new std::vector<Tau>;

vector<double>* Tau_pt = new std::vector<double>;
vector<double>* Tau_eta = new std::vector<double>;
vector<double>* Tau_phi = new std::vector<double>;
vector<double>* Tau_energy = new std::vector<double>;
vector<double>* Tau_charge = new std::vector<double>;
vector<double>* Tau_packedLeadTauCand_dz = new std::vector<double>;
vector<double>* Tau_packedLeadTauCand_dxy = new std::vector<double>;
vector<double>* Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT = new std::vector<double>;
vector<double>* Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT = new std::vector<double>;
vector<double>* Tau_decayModeFinding = new std::vector<double>;
vector<double>* Tau_cut = new std::vector<double>;

// Jet
vector<Jet>* jets = new std::vector<Jet>;

vector<double>* Jet_cut = new std::vector<double>;
vector<double>* Jet_bcut = new std::vector<double>;
vector<double>* Jet_pt = new std::vector<double>;
vector<double>* Jet_eta = new std::vector<double>;
vector<double>* Jet_phi = new std::vector<double>;
vector<double>* Jet_energy = new std::vector<double>;
vector<double>* Jet_genMother_pt = new std::vector<double>;
vector<double>* Jet_genMother_eta = new std::vector<double>;
vector<double>* Jet_genMother_phi = new std::vector<double>;
vector<double>* Jet_genMother_en = new std::vector<double>;
vector<double>* Jet_genMother_pdgId = new std::vector<double>;
vector<double>* Jet_genGrandMother_pt = new std::vector<double>;
vector<double>* Jet_genGrandMother_eta = new std::vector<double>;
vector<double>* Jet_genGrandMother_phi = new std::vector<double>;
vector<double>* Jet_genGrandMother_en = new std::vector<double>;
vector<double>* Jet_genGrandMother_pdgId = new std::vector<double>;
vector<double>* Jet_Uncorr_pt = new std::vector<double>;
vector<double>* Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags = new std::vector<double>;
vector<double>* Jet_pfCombinedMVAV2BJetTags = new std::vector<double>;
vector<double>* Jet_qg = new std::vector<double>;
vector<double>* Jet_axis2 = new std::vector<double>;
vector<double>* Jet_ptD = new std::vector<double>;
vector<double>* Jet_mult = new std::vector<double>;
vector<double>* Jet_partonFlavour = new std::vector<double>;
vector<double>* Jet_hadronFlavour = new std::vector<double>;
vector<double>* Jet_genpt = new std::vector<double>;
vector<double>* Jet_geneta = new std::vector<double>;
vector<double>* Jet_genphi = new std::vector<double>;
vector<double>* Jet_genenergy = new std::vector<double>;

//BoostJet
vector<BoostJet>* boostjets = new std::vector<BoostJet>;

vector<double>* BoostedJet_pt = new std::vector<double>;
vector<double>* BoostedJet_eta = new std::vector<double>;
vector<double>* BoostedJet_phi = new std::vector<double>;
vector<double>* BoostedJet_energy = new std::vector<double>;
vector<double>* BoostedJet_Uncorr_pt = new std::vector<double>;
vector<double>* BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags = new std::vector<double>;
vector<double>* BoostedJet_pfCombinedMVAV2BJetTags = new std::vector<double>;
vector<double>* BoostedJet_JesSF = new std::vector<double>;
vector<double>* BoostedJet_JesSFup = new std::vector<double>;
vector<double>* BoostedJet_JesSFdown = new std::vector<double>;
vector<double>* BoostedJet_JerSF = new std::vector<double>;
vector<double>* BoostedJet_JerSFup = new std::vector<double>;
vector<double>* BoostedJet_JerSFdown = new std::vector<double>;
vector<double>* BoostedJet_tau1 = new std::vector<double>;
vector<double>* BoostedJet_tau2 = new std::vector<double>;
vector<double>* BoostedJet_tau3 = new std::vector<double>;
vector<double>* BoostedJet_softdrop_mass = new std::vector<double>;
vector<double>* BoostedJet_pruned_mass = new std::vector<double>;
vector<double>* BoostedJet_tau21 = new std::vector<double>;
vector<double>* BoostedJet_tau32 = new std::vector<double>;
vector<double>* BoostedJet_wCut = new std::vector<double>;
vector<double>* BoostedJet_topCut = new std::vector<double>;

//Gen
vector<double>* Gen_pdg_id = new std::vector<double>;
vector<double>* Gen_pt = new std::vector<double>;
vector<double>* Gen_eta = new std::vector<double>;
vector<double>* Gen_phi = new std::vector<double>;
vector<double>* Gen_energy = new std::vector<double>;
vector<double>* Gen_motherpdg_id = new std::vector<double>;
vector<double>* Gen_BmotherIndex = new std::vector<double>;
vector<double>* Gen_numMother = new std::vector<double>;
vector<double>* Gen_status = new std::vector<double>;
//new hadronic Top and W
vector<double>* hadTop_Gen_pt = new std::vector<double>;
vector<double>* hadTop_Gen_eta = new std::vector<double>;
vector<double>* hadTop_Gen_phi = new std::vector<double>;
vector<double>* hadTop_Gen_energy = new std::vector<double>;
vector<double>* hadTop_Gen_pdgId = new std::vector<double>;
vector<double>* hadTop_Gen_Index = new std::vector<double>;
vector<double>* hadW_Gen_pt = new std::vector<double>;
vector<double>* hadW_Gen_eta = new std::vector<double>;
vector<double>* hadW_Gen_phi = new std::vector<double>;
vector<double>* hadW_Gen_energy = new std::vector<double>;
vector<double>* hadW_Gen_pdgId = new std::vector<double>;
vector<double>* hadW_Gen_Index = new std::vector<double>;
