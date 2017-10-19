#include "Rootplizer_TTHLep_v2IHEP.h"
/////
//   Main function
/////
void Rootplizer_TTHLep_v2IHEP(const char * Input = "", const char * Output ="", const char* sync_type =""){
    synchro = std::string(sync_type); 
 
    //Call input file
    TFile *inputfile = TFile::Open(Input);
    TTree *readingtree = new TTree("readingtree","readingtree"); readingtree = (TTree*) inputfile->Get("TNT/BOOM");
    if((string(Output).find("SMu") != std::string::npos) || (string(Output).find("SEle") != std::string::npos)
    || (string(Output).find("DMu") != std::string::npos) || (string(Output).find("DEleGm") != std::string::npos) || (string(Output).find("MuEleGm") != std::string::npos)
      ) sample = "data";
    if(
       (string(Output).find("SMuG") != std::string::npos) || (string(Output).find("SEleG") != std::string::npos) || (string(Output).find("SMuH") != std::string::npos) || (string(Output).find("SEleH") != std::string::npos)
    || (string(Output).find("DMuH") != std::string::npos) || (string(Output).find("DEleGmH") != std::string::npos) || (string(Output).find("MuEleGmH") != std::string::npos)
    || (string(Output).find("DMuG") != std::string::npos) || (string(Output).find("DEleGmG") != std::string::npos) || (string(Output).find("MuEleGmG") != std::string::npos)
     ) isHIPSafe = false;
    rSetBranchAddress(readingtree,sample);
 
    //Define output file
    TFile *newfile = new TFile(Output,"recreate");
    TTree* newtree = new TTree("BOOM","BOOM");
    newtree->SetMaxTreeSize(99000000000);
    wSetBranchAddress(newtree,sample);
 
    // Set MVA Weights
    set_wgtMVA();
     
    //Fill new branches
    int nen = nentries; if(nentries==-1) nen = readingtree->GetEntries();
    for(Int_t en=0; en<nen; en++){
        //Initialization
        wClearInitialization(sample);
        Long64_t tentry = readingtree->LoadTree(en); 
        rGetEntry(tentry, sample);
        //Muon
        Muon_sel(sample);
        //Electron
        patElectron_sel(sample);
        //Tau
        Tau_sel();
        //Jet
        Jet_sel(sample);
        //BoostJet
        BoostedJet_sel();
        Lep_sel();
        //Event
        Event_sel();
        newtree->Fill();
    }
 
    //Save new file
    newfile->cd();
    newfile->Write();
    newfile->Close();
}
////
// Anaysis functions
////
void set_wgtMVA(){
    string ele_wgt = data_path + "el_BDTG.weights.xml";
    string mu_wgt = data_path + "mu_BDTG.weights.xml";
    mu_reader_ = new TMVA::Reader("!Color:!Silent");
    ele_reader_ = new TMVA::Reader("!Color:!Silent");
    std::vector<TMVA::Reader *> mvas = { ele_reader_, mu_reader_ };
    for (auto &m : mvas) {
        m->AddVariable("LepGood_pt", &varpt);
        m->AddVariable("LepGood_eta", &vareta);
        m->AddVariable("LepGood_jetNDauChargedMVASel", &varjetNDauCharged_in);
        m->AddVariable("LepGood_miniRelIsoCharged", &varchRelIso);
        m->AddVariable("LepGood_miniRelIsoNeutral", &varneuRelIso);
        m->AddVariable("LepGood_jetPtRelv2", &varjetPtRel_in);
        m->AddVariable("LepGood_jetPtRatio := min(LepGood_jetPtRatiov2,1.5)", &varjetPtRatio_in);
        m->AddVariable("LepGood_jetBTagCSV := max(LepGood_jetBTagCSV,0)", &varjetBTagCSV_in);
        m->AddVariable("LepGood_sip3d", &varsip3d);
        m->AddVariable("LepGood_dxy := log(abs(LepGood_dxy))", &vardxy);
        m->AddVariable("LepGood_dz := log(abs(LepGood_dz))", &vardz);
    }
    ele_reader_->AddVariable("LepGood_mvaIdSpring16HZZ", &varmvaId);
    mu_reader_->AddVariable("LepGood_segmentCompatibility", &varSegCompat);
    ele_reader_->BookMVA("BDTG method",ele_wgt); 
    mu_reader_->BookMVA("BDTG method",mu_wgt); 
}; 


double get_LeptonMVA(Lepton lep){ 
    varpt = lep.pt;
    vareta = lep.eta;
    varchRelIso = lep.miniIsoCh/lep.pt;
    varneuRelIso = lep.miniIsoPUsub/lep.pt;
    varjetPtRel_in = lep.ptrel;
    varjetPtRatio_in = min(lep.jetptratio,1.5);
    varjetBTagCSV_in = max(lep.jetcsv,0.); 
    varjetNDauCharged_in =lep.lepjetchtrks; 
    varsip3d = lep.IP3Dsig;
    vardxy = log(abs(lep.dxy_pv)); 
    vardz =  log(abs(lep.dz_pv)); 
    varSegCompat = lep.segmentCompatibility; 
    varmvaId = lep.mvaValue_HZZ;
    if( fabs(lep.pdgId) == 13) return mu_reader_->EvaluateMVA("BDTG method");
    else if( fabs(lep.pdgId) == 11) return ele_reader_->EvaluateMVA("BDTG method");
    else return -1.;
}

// Analysis functions

bool mu_isMedium( double isGlobal, double chi_square, double chi2_localposition, double trkKink, double mu_isloose, double validFrac, double segment ,bool ishipsave){
    bool GoodGlobal = (isGlobal==1 && chi_square <3 && chi2_localposition < 12 && trkKink < 20);
    bool isMedium = mu_isloose ==1 && validFrac > (ishipsave? 0.49 : 0.80) && segment > ( GoodGlobal? 0.303 : 0.451);
    return isMedium;
};


//utils
double deltaPhi(double phi1, double phi2){
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}
double deltaEta(double eta1, double eta2){
    return (eta1-eta2);
};
double deltaR(double dphi, double deta){
    return sqrt(pow(dphi,2)+pow(deta,2));
};
bool byPt(const Lepton& LeptonA, const Lepton& LeptonB){
    return LeptonA.pt > LeptonB.pt;
};


////
//Object Selection
////
//Muon
void Muon_sel(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    int mu_numLoose = 0;
    int mu_numFake = 0;
    int mu_numTight = 0;
    for(uint mu_en = 0; mu_en<rMuon_pt->size(); mu_en++){
        Lepton Muon;
        
        // initialize needed variables for mu_isLoose_tthlep()
        Muon.pt= rMuon_pt->at(mu_en);
        Muon.eta= rMuon_eta->at(mu_en);
        Muon.phi= rMuon_phi->at(mu_en);
        Muon.energy= rMuon_energy->at(mu_en);
        Muon.dxy_pv= rMuon_dxy_pv->at(mu_en);
        Muon.dz_pv= rMuon_dz_pv->at(mu_en);
        Muon.IP3Dsig= rMuon_IP3Dsig_it->at(mu_en);
        Muon.loose= rMuon_loose->at(mu_en);
        Muon.miniIsoRel= rMuon_miniIsoRel->at(mu_en);
        
        // check whether Muon pass tth loose lepton selection
        if(!(Muon.mu_isLoose_tthlep())) continue;
        
        Muon.isGlobal= rMuon_isGlobal->at(mu_en);
        Muon.chi2= rMuon_chi2->at(mu_en);
        Muon.chi2LocalPosition= rMuon_chi2LocalPosition->at(mu_en);
        Muon.trkKink= rMuon_trkKink->at(mu_en);
        Muon.validFraction= rMuon_validFraction->at(mu_en);
        Muon.segmentCompatibility= rMuon_segmentCompatibility->at(mu_en);
        Muon.jetptratio= rMuon_jetptratio->at(mu_en);
        Muon.jetcsv= rMuon_jetcsv->at(mu_en);
        Muon.jetpt= rMuon_jetpt->at(mu_en);
        Muon.lepjetchtrks= rMuon_lepjetchtrks->at(mu_en);
        Muon.miniIsoCh= rMuon_miniIsoCh->at(mu_en);
        Muon.miniIsoPUsub= rMuon_miniIsoPUsub->at(mu_en);
        Muon.ptrel= rMuon_ptrel->at(mu_en);
        Muon.pTErrOVpT_it= rMuon_pTErrOVpT_it->at(mu_en);
        Muon.px= rMuon_px->at(mu_en);
        Muon.py= rMuon_py->at(mu_en);
        Muon.pz= rMuon_pz->at(mu_en);
        Muon.jetdr= rMuon_jetdr->at(mu_en);
        Muon.charge= rMuon_charge->at(mu_en);
        Muon.pdgId= rMuon_pdgId->at(mu_en);
        
        // few variables needs to be seted before set BDT
        Muon.BDT = get_LeptonMVA(Muon);
        bool isMedium_ST = mu_isMedium( rMuon_isGlobal->at(mu_en), rMuon_chi2->at(mu_en), rMuon_chi2LocalPosition->at(mu_en) , rMuon_trkKink->at(mu_en), rMuon_loose->at(mu_en), rMuon_validFraction->at(mu_en) , rMuon_segmentCompatibility->at(mu_en) , isHIPSafe);
        
        // pt, jetpt and BDT of Muon has to be seted before calling conept
        Muon.cal_conept(isMedium_ST);
       
        Muon.gen_pt= rMuon_gen_pt->at(mu_en);
        Muon.gen_eta= rMuon_gen_eta->at(mu_en);
        Muon.gen_phi= rMuon_gen_phi->at(mu_en);
        Muon.gen_en= rMuon_gen_en->at(mu_en);
        Muon.gen_pdgId= rMuon_gen_pdgId->at(mu_en);
        Muon.genMother_pt= rMuon_genMother_pt->at(mu_en);
        Muon.genMother_eta= rMuon_genMother_eta->at(mu_en);
        Muon.genMother_phi= rMuon_genMother_phi->at(mu_en);
        Muon.genMother_en= rMuon_genMother_en->at(mu_en);
        Muon.genMother_pdgId= rMuon_genMother_pdgId->at(mu_en);
        Muon.genGrandMother_pt= rMuon_genGrandMother_pt->at(mu_en);
        Muon.genGrandMother_eta= rMuon_genGrandMother_eta->at(mu_en);
        Muon.genGrandMother_phi= rMuon_genGrandMother_phi->at(mu_en);
        Muon.genGrandMother_en= rMuon_genGrandMother_en->at(mu_en);
        Muon.genGrandMother_pdgId= rMuon_genGrandMother_pdgId->at(mu_en);
        Muon.gen_isPromptFinalState= rMuon_gen_isPromptFinalState->at(mu_en);
        Muon.gen_isDirectPromptTauDecayProductFinalState= rMuon_gen_isDirectPromptTauDecayProductFinalState->at(mu_en);
       
        // calculate new variables 
        Muon.set_Wp_tthlep( isMedium_ST, mu_numLoose, mu_numFake, mu_numTight );
        Muon.cal_tight_property();
        Muon.CF = 0.;
        Muon.FR = Muon.get_valX_valY_binContent(hist_mu_fr, Muon.corrpt, Muon.eta);
        if(!(Muon.gen_pdgId * Muon.charge<0 && Muon.gen_pdgId!=-999))Muon.isMatchRightCharge=0.;
        leptons->push_back(Muon);    
    }
    Muon_numLoose = mu_numLoose;
    Muon_numFake = mu_numFake;
    Muon_numTight = mu_numTight;
}


//Electron
void patElectron_sel(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    int ele_numLoose = 0;
    int ele_numFake = 0;
    int ele_numTight = 0;
    for(uint ele_en = 0; ele_en<rpatElectron_pt->size(); ele_en++){
        Lepton patElectron;
        
        // initialize needed variables for ele_isLoose_tthlep()
        patElectron.pt= rpatElectron_pt->at(ele_en);
        patElectron.eta= rpatElectron_eta->at(ele_en);
        patElectron.phi= rpatElectron_phi->at(ele_en);
        patElectron.energy= rpatElectron_energy->at(ele_en);
        patElectron.IP3Dsig= rpatElectron_IP3Dsig->at(ele_en);
        patElectron.miniIsoRel= rpatElectron_miniIsoRel->at(ele_en);
        patElectron.dxy_pv= rpatElectron_gsfTrack_dxy_pv->at(ele_en);
        patElectron.dz_pv= rpatElectron_gsfTrack_dz_pv->at(ele_en);
        patElectron.SCeta= rpatElectron_SCeta->at(ele_en);
        patElectron.expectedMissingInnerHits= rpatElectron_expectedMissingInnerHits->at(ele_en);
        
        // check whether ele overlaps with tth loose muon
        bool ismatched = false; 
        for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
            Lepton Muon = leptons->at(lep_en);
            if(!(fabs(Muon.pdgId)==13)) continue; // skip if it's not a muon
            if(deltaR(deltaPhi(Muon.phi,patElectron.phi),deltaEta(Muon.eta,patElectron.eta))<0.05){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
        
        // check whether ele pass tth loose lepton selection
        if(!(patElectron.ele_isLoose_tthlep())) continue;
        
        // few variables needs to be seted before set BDT
        patElectron.charge= rpatElectron_charge->at(ele_en);
        patElectron.pdgId= rpatElectron_pdgId->at(ele_en);
        patElectron.jetptratio= rpatElectron_jetptratio->at(ele_en);
        patElectron.jetcsv= rpatElectron_jetcsv->at(ele_en);
        patElectron.jetpt= rpatElectron_jetpt->at(ele_en);
        patElectron.lepjetchtrks= rpatElectron_lepjetchtrks->at(ele_en);
        patElectron.miniIsoCh= rpatElectron_miniIsoCh->at(ele_en);
        patElectron.miniIsoPUsub= rpatElectron_miniIsoPUsub->at(ele_en);
        patElectron.ptrel= rpatElectron_ptrel->at(ele_en);
        patElectron.px= rpatElectron_px->at(ele_en);
        patElectron.py= rpatElectron_py->at(ele_en);
        patElectron.pz= rpatElectron_pz->at(ele_en);
        patElectron.jetdr= rpatElectron_jetdr->at(ele_en);
        patElectron.mvaValue_HZZ= rpatElectron_mvaValue_HZZ->at(ele_en);
        patElectron.BDT = get_LeptonMVA(patElectron);
    
        // pt, jetpt and BDT of electron has to be seted before calling conept
        // pass ismedium boolean, always true for electron
        patElectron.cal_conept(true); 
       
        // few variables needs to be seted before check fabkeable selection
        patElectron.full5x5_sigmaIetaIeta= rpatElectron_full5x5_sigmaIetaIeta->at(ele_en);
        patElectron.hOverE= rpatElectron_hOverE->at(ele_en);
        patElectron.dEtaIn= rpatElectron_dEtaIn->at(ele_en);
        patElectron.dPhiIn= rpatElectron_dPhiIn->at(ele_en);
        patElectron.ooEmooP= rpatElectron_ooEmooP->at(ele_en);
       
        patElectron.isGsfCtfScPixChargeConsistent= rpatElectron_isGsfCtfScPixChargeConsistent->at(ele_en);
        patElectron.isGsfScPixChargeConsistent= rpatElectron_isGsfScPixChargeConsistent->at(ele_en);
        patElectron.passConversion= rpatElectron_passConversionVeto->at(ele_en);
        patElectron.gen_pt= rpatElectron_gen_pt->at(ele_en);
        patElectron.gen_eta= rpatElectron_gen_eta->at(ele_en);
        patElectron.gen_phi= rpatElectron_gen_phi->at(ele_en);
        patElectron.gen_en= rpatElectron_gen_en->at(ele_en);
        patElectron.gen_pdgId= rpatElectron_gen_pdgId->at(ele_en);
        patElectron.genMother_pt= rpatElectron_genMother_pt->at(ele_en);
        patElectron.genMother_eta= rpatElectron_genMother_eta->at(ele_en);
        patElectron.genMother_phi= rpatElectron_genMother_phi->at(ele_en);
        patElectron.genMother_en= rpatElectron_genMother_en->at(ele_en);
        patElectron.genMother_pdgId= rpatElectron_genMother_pdgId->at(ele_en);
        patElectron.genGrandMother_pt= rpatElectron_genGrandMother_pt->at(ele_en);
        patElectron.genGrandMother_eta= rpatElectron_genGrandMother_eta->at(ele_en);
        patElectron.genGrandMother_phi= rpatElectron_genGrandMother_phi->at(ele_en);
        patElectron.genGrandMother_en= rpatElectron_genGrandMother_en->at(ele_en);
        patElectron.genGrandMother_pdgId= rpatElectron_genGrandMother_pdgId->at(ele_en);
        patElectron.gen_isPromptFinalState= rpatElectron_gen_isPromptFinalState->at(ele_en);
        patElectron.gen_isDirectPromptTauDecayProductFinalState= rpatElectron_gen_isDirectPromptTauDecayProductFinalState->at(ele_en);
        
        // calculate new variables 
        // pass ismedium boolean, always true for electron
        patElectron.set_Wp_tthlep(true, ele_numLoose, ele_numFake, ele_numTight ); 
        patElectron.cal_tight_property();
        patElectron.CF = patElectron.get_valX_valY_binContent(
            hist_cf, patElectron.corrpt,patElectron.eta);
        patElectron.FR = patElectron.get_valX_valY_binContent(
            hist_el_fr, patElectron.corrpt,patElectron.eta);
        if(!(patElectron.gen_pdgId * patElectron.charge<0 && patElectron.gen_pdgId!=-999))
            patElectron.isMatchRightCharge=0.;
        leptons->push_back(patElectron);    
    }
    patElectron_numLoose = ele_numLoose;
    patElectron_numFake =  ele_numFake;
    patElectron_numTight = ele_numTight;
}


//Tau
void Tau_sel(){
    int tau_numLoose = 0;
    int tau_numMedium = 0;
    for(uint tau_en = 0; tau_en<rTau_pt->size(); tau_en++){
        Tau tau;
        tau.pt= rTau_pt->at(tau_en);
        tau.eta= rTau_eta->at(tau_en);
        tau.phi= rTau_phi->at(tau_en);
        tau.energy= rTau_energy->at(tau_en);
        tau.charge= rTau_charge->at(tau_en);
        tau.packedLeadTauCand_dz= rTau_packedLeadTauCand_dz->at(tau_en);
        tau.packedLeadTauCand_dxy= rTau_packedLeadTauCand_dxy->at(tau_en);
        tau.byLooseIsolationMVArun2v1DBdR03oldDMwLT= rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(tau_en);
        tau.byMediumIsolationMVArun2v1DBdR03oldDMwLT= rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(tau_en);
        tau.decayModeFinding= rTau_decayModeFinding->at(tau_en);
        
        // check whether tau overlaps with tth loose leptons
        bool ismatched = false; 
        for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
            Lepton lep = leptons->at(lep_en);
            if(deltaR(deltaPhi(lep.phi,tau.phi),deltaEta(lep.eta,tau.eta))<0.4){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
        
        // check whether ele pass tth loose lepton selection
        if(!(tau.tau_isLoose_tthlep())) continue;

        // calculate new variables 
        tau.set_Wp_tthlep(tau_numLoose, tau_numMedium);
        Tau_pt->push_back(rTau_pt->at(tau_en));
        Tau_eta->push_back(rTau_eta->at(tau_en));
        Tau_phi->push_back(rTau_phi->at(tau_en));
        Tau_energy->push_back(rTau_energy->at(tau_en));
        Tau_charge->push_back(rTau_charge->at(tau_en));
        Tau_packedLeadTauCand_dxy->push_back(rTau_packedLeadTauCand_dxy->at(tau_en));
        Tau_packedLeadTauCand_dz->push_back(rTau_packedLeadTauCand_dz->at(tau_en));
        Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->push_back(rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(tau_en));
        Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->push_back(rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(tau_en));
        Tau_decayModeFinding->push_back(rTau_decayModeFinding->at(tau_en));
        Tau_cut->push_back(tau.cut);
        taus->push_back(tau);    
    }
    Tau_numMedium = tau_numMedium; 
    Tau_numLoose = tau_numLoose; 
}


void Jet_sel(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    int jet_numLoose   = 0;
    int jet_numbLoose = 0;
    int jet_numbMedium = 0;
    int jet_numbTight = 0;
    double sf_BWeight= 1.;
    double sf_BWeightLFup= 1.;
    double sf_BWeightLFdown= 1.;
    double sf_BWeightHFup= 1.;
    double sf_BWeightHFdown= 1.;
    double sf_BWeightHFStats1up= 1.;
    double sf_BWeightHFStats1down= 1.;
    double sf_BWeightLFStats1up= 1.;
    double sf_BWeightLFStats1down= 1.;
    double sf_BWeightHFStats2up= 1.;
    double sf_BWeightHFStats2down= 1.;
    double sf_BWeightLFStats2up= 1.;
    double sf_BWeightLFStats2down= 1.;
    double sf_BWeightCErr1up= 1.;
    double sf_BWeightCErr1down= 1.;
    double sf_BWeightCErr2up= 1.;
    double sf_BWeightCErr2down= 1.;
    double sf_BWeightJESup= 1.;
    double sf_BWeightJESdown= 1.;
    for(uint jet_en = 0; jet_en<rJet_pt->size(); jet_en++){
        Jet jet;
        jet.pt = rJet_Uncorr_pt->at(jet_en)*rJet_JesSF->at(jet_en);
        jet.energy = rJet_energy->at(jet_en)*rJet_Uncorr_pt->at(jet_en)/rJet_pt->at(jet_en)*rJet_JesSF->at(jet_en);
        jet.eta = rJet_eta->at(jet_en);
        jet.phi = rJet_phi->at(jet_en);
        
        // check whether jet overlaps with tth loose leptons
        bool ismatched = false; 
        for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
            Lepton lep = leptons->at(lep_en);
            if(lep.cut <=1) continue; // skip loose lepton
            if(deltaR(deltaPhi(lep.phi,jet.phi),deltaEta(lep.eta,jet.eta))<0.4){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
        
        // check whether jet overlaps with tau
        for(uint tau_en=0; tau_en < taus->size(); tau_en++){
            Tau tau = taus->at(tau_en);
            if(deltaR(deltaPhi(tau.phi,jet.phi),deltaEta(tau.eta,jet.eta))<0.4){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
       
        // set variables neeeded for jet_isLoose()
        jet.neutralHadEnergyFraction= rJet_neutralHadEnergyFraction->at(jet_en);
        jet.neutralEmEnergyFraction= rJet_neutralEmEnergyFraction->at(jet_en);
        jet.chargedMultiplicity= rJet_chargedMultiplicity->at(jet_en);
        jet.numberOfConstituents= rJet_numberOfConstituents->at(jet_en);
        jet.chargedHadronEnergyFraction= rJet_chargedHadronEnergyFraction->at(jet_en);
        jet.chargedEmEnergyFraction= rJet_chargedEmEnergyFraction->at(jet_en);
        
        // check whether jet pass loose jet selection 
        if(!jet.jet_isLoose()) continue;

        jet.JesSF= rJet_JesSF->at(jet_en);
        jet.JesSFup= rJet_JesSFup->at(jet_en);
        jet.JesSFdown= rJet_JesSFdown->at(jet_en);
        jet.JerSF= rJet_JerSF->at(jet_en);
        jet.JerSFup= rJet_JerSFup->at(jet_en);
        jet.JerSFdown= rJet_JerSFdown->at(jet_en);
        jet.btag_sf= rJet_btag_sf->at(jet_en);
        jet.btag_jesup= rJet_btag_jesup->at(jet_en);
        jet.btag_jesdown= rJet_btag_jesdown->at(jet_en);
        jet.btag_hfup= rJet_btag_hfup->at(jet_en);
        jet.btag_hfdown= rJet_btag_hfdown->at(jet_en);
        jet.btag_hfstat1up= rJet_btag_hfstat1up->at(jet_en);
        jet.btag_hfstat1down= rJet_btag_hfstat1down->at(jet_en);
        jet.btag_hfstat2up= rJet_btag_hfstat2up->at(jet_en);
        jet.btag_hfstat2down= rJet_btag_hfstat2down->at(jet_en);
        jet.btag_lfup= rJet_btag_lfup->at(jet_en);
        jet.btag_lfdown= rJet_btag_lfdown->at(jet_en);
        jet.btag_lfstat1up= rJet_btag_lfstat1up->at(jet_en);
        jet.btag_lfstat1down= rJet_btag_lfstat1down->at(jet_en);
        jet.btag_lfstat2up= rJet_btag_lfstat2up->at(jet_en);
        jet.btag_lfstat2down= rJet_btag_lfstat2down->at(jet_en);
        jet.btag_cerr1up= rJet_btag_cerr1up->at(jet_en);
        jet.btag_cerr1down= rJet_btag_cerr1down->at(jet_en);
        jet.btag_cerr2up= rJet_btag_cerr2up->at(jet_en);
        jet.btag_cerr2down= rJet_btag_cerr2down->at(jet_en);
        jet.genMother_pt= rJet_genMother_pt->at(jet_en);
        jet.genMother_eta= rJet_genMother_eta->at(jet_en);
        jet.genMother_phi= rJet_genMother_phi->at(jet_en);
        jet.genMother_en= rJet_genMother_en->at(jet_en);
        jet.genMother_pdgId= rJet_genMother_pdgId->at(jet_en);
        jet.genGrandMother_pt= rJet_genGrandMother_pt->at(jet_en);
        jet.genGrandMother_eta= rJet_genGrandMother_eta->at(jet_en);
        jet.genGrandMother_phi= rJet_genGrandMother_phi->at(jet_en);
        jet.genGrandMother_en= rJet_genGrandMother_en->at(jet_en);
        jet.genGrandMother_pdgId= rJet_genGrandMother_pdgId->at(jet_en);
        jet.Uncorr_pt= rJet_Uncorr_pt->at(jet_en);
        jet.pfCombinedInclusiveSecondaryVertexV2BJetTags= rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->at(jet_en);
        jet.pfCombinedMVAV2BJetTags= rJet_pfCombinedMVAV2BJetTags->at(jet_en);
        jet.qg= rJet_qg->at(jet_en);
        jet.axis2= rJet_axis2->at(jet_en);
        jet.ptD= rJet_ptD->at(jet_en);
        jet.mult= rJet_mult->at(jet_en);
        jet.partonFlavour= rJet_partonFlavour->at(jet_en);
        jet.hadronFlavour= rJet_hadronFlavour->at(jet_en);
        jet.genpt= rJet_genpt->at(jet_en);
        jet.geneta= rJet_geneta->at(jet_en);
        jet.genphi= rJet_genphi->at(jet_en);
        jet.genenergy= rJet_genenergy->at(jet_en);
        
        // calculate new variables 
        jet.set_Wp_jets(jet_numLoose);
        jet.set_Wp_bdisc(jet_numbLoose, jet_numbMedium, jet_numbTight);
         
        Jet_cut->push_back(jet.cut);
        Jet_bcut->push_back(jet.bcut);
        Jet_qg->push_back(rJet_qg->at(jet_en));
        Jet_axis2->push_back(rJet_axis2->at(jet_en));
        Jet_ptD->push_back(rJet_ptD->at(jet_en));
        Jet_mult->push_back(rJet_mult->at(jet_en));
        Jet_partonFlavour->push_back(rJet_partonFlavour->at(jet_en));
        Jet_hadronFlavour->push_back(rJet_hadronFlavour->at(jet_en));
        Jet_genpt->push_back(rJet_genpt->at(jet_en));
        Jet_geneta->push_back(rJet_geneta->at(jet_en));
        Jet_genphi->push_back(rJet_genphi->at(jet_en));
        Jet_genenergy->push_back(rJet_genenergy->at(jet_en));
        Jet_Uncorr_pt->push_back(rJet_Uncorr_pt->at(jet_en));
        Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags->push_back(rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->at(jet_en));
        Jet_pfCombinedMVAV2BJetTags->push_back(rJet_pfCombinedMVAV2BJetTags->at(jet_en));
        Jet_genMother_pt->push_back(rJet_genMother_pt->at(jet_en));
        Jet_genMother_eta->push_back(rJet_genMother_eta->at(jet_en));
        Jet_genMother_phi->push_back(rJet_genMother_phi->at(jet_en));
        Jet_genMother_en->push_back(rJet_genMother_en->at(jet_en));
        Jet_genMother_pdgId->push_back(rJet_genMother_pdgId->at(jet_en));
        Jet_genGrandMother_pt->push_back(rJet_genGrandMother_pt->at(jet_en));
        Jet_genGrandMother_eta->push_back(rJet_genGrandMother_eta->at(jet_en));
        Jet_genGrandMother_phi->push_back(rJet_genGrandMother_phi->at(jet_en));
        Jet_genGrandMother_en->push_back(rJet_genGrandMother_en->at(jet_en));
        Jet_genGrandMother_pdgId->push_back(rJet_genGrandMother_pdgId->at(jet_en));
        Jet_pt->push_back(jet.pt);
        Jet_eta->push_back(rJet_eta->at(jet_en));
        Jet_phi->push_back(rJet_phi->at(jet_en));
        Jet_energy->push_back(jet.energy);
        
        jets->push_back(jet);
  
        sf_BWeight *=rJet_btag_sf->at(jet_en);
        sf_BWeightLFup *= rJet_btag_lfup->at(jet_en);
        sf_BWeightLFdown *= rJet_btag_lfdown->at(jet_en);
        sf_BWeightHFup *= rJet_btag_hfup->at(jet_en);
        sf_BWeightHFdown *= rJet_btag_hfdown->at(jet_en);
        sf_BWeightHFStats1up *= rJet_btag_hfstat1up->at(jet_en);
        sf_BWeightHFStats1down *= rJet_btag_hfstat1down->at(jet_en);
        sf_BWeightLFStats1up *= rJet_btag_lfstat1up->at(jet_en);
        sf_BWeightLFStats1down *= rJet_btag_lfstat1down->at(jet_en);
        sf_BWeightHFStats2up *= rJet_btag_hfstat2up->at(jet_en);
        sf_BWeightHFStats2down *= rJet_btag_hfstat2down->at(jet_en);
        sf_BWeightLFStats2up *= rJet_btag_lfstat2up->at(jet_en);
        sf_BWeightLFStats2down *= rJet_btag_lfstat2down->at(jet_en);
        sf_BWeightCErr1up *= rJet_btag_cerr1up->at(jet_en);
        sf_BWeightCErr1down *= rJet_btag_cerr1down->at(jet_en);
        sf_BWeightCErr2up *= rJet_btag_cerr2up->at(jet_en);
        sf_BWeightCErr2down *= rJet_btag_cerr2down->at(jet_en);
        sf_BWeightJESup *= rJet_btag_jesup->at(jet_en);
        sf_BWeightJESdown *= rJet_btag_jesdown->at(jet_en);
    }
    BWeight = sf_BWeight;
    BWeightLFup = sf_BWeightLFup;
    BWeightLFdown = sf_BWeightLFdown;
    BWeightHFup = sf_BWeightHFup;
    BWeightHFdown = sf_BWeightHFdown;
    BWeightHFStats1up = sf_BWeightHFStats1up;
    BWeightHFStats1down = sf_BWeightHFStats1down;
    BWeightLFStats1up = sf_BWeightLFStats1up;
    BWeightLFStats1down = sf_BWeightLFStats1down;
    BWeightHFStats2up = sf_BWeightHFStats2up;
    BWeightHFStats2down = sf_BWeightHFStats2down;
    BWeightLFStats2up = sf_BWeightLFStats2up;
    BWeightLFStats2down = sf_BWeightLFStats2down;
    BWeightCErr1up = sf_BWeightCErr1up;
    BWeightCErr1down = sf_BWeightCErr1down;
    BWeightCErr2up = sf_BWeightCErr2up;
    BWeightCErr2down = sf_BWeightCErr2down;
    BWeightJESup = sf_BWeightJESup;
    BWeightJESdown = sf_BWeightJESdown;
    Jet_numLoose = jet_numLoose;
    Jet_numbLoose = jet_numbLoose;
    Jet_numbMedium = jet_numbMedium;
    Jet_numbTight = jet_numbTight;
}; 


void BoostedJet_sel(){
    int top_numSoft   = 0;
    int top_numLoose = 0;
    int top_numMedium = 0;
    int top_numTight = 0;
    int w_numLoose = 0;
    int w_numTight = 0;
    for(uint boostedjet_en = 0; boostedjet_en<rBoostedJet_pt->size(); boostedjet_en++){
        BoostJet boostjet;
        boostjet.pt = rBoostedJet_Uncorr_pt->at(boostedjet_en)*rBoostedJet_JesSF->at(boostedjet_en);
        boostjet.energy = rBoostedJet_energy->at(boostedjet_en)*rBoostedJet_Uncorr_pt->at(boostedjet_en)/rBoostedJet_pt->at(boostedjet_en)*rBoostedJet_JesSF->at(boostedjet_en);
        boostjet.eta = rBoostedJet_eta->at(boostedjet_en);
        boostjet.phi = rBoostedJet_phi->at(boostedjet_en);
        // check whether boostedjet overlaps with tth loose leptons
        bool ismatched = false; 
        for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
            Lepton lep = leptons->at(lep_en);
            if(lep.cut <=1) continue; // skip loose lepton
            if(deltaR(deltaPhi(lep.phi,boostjet.phi),deltaEta(lep.eta,boostjet.eta))<0.4){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
        
        // check whether boostedjet overlaps with tau
        for(uint tau_en=0; tau_en < taus->size(); tau_en++){
            Tau tau = taus->at(tau_en);
            if(deltaR(deltaPhi(tau.phi,boostjet.phi),deltaEta(tau.eta,boostjet.eta))<0.4){
                ismatched = true; 
                break; 
            }
        }
        if(ismatched) continue;
       
        boostjet.Uncorr_pt= rBoostedJet_Uncorr_pt->at(boostedjet_en);
        boostjet.pfCombinedInclusiveSecondaryVertexV2BJetTags= rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->at(boostedjet_en);
        boostjet.pfCombinedMVAV2BJetTags= rBoostedJet_pfCombinedMVAV2BJetTags->at(boostedjet_en);
        boostjet.JesSF= rBoostedJet_JesSF->at(boostedjet_en);
        boostjet.JesSFup= rBoostedJet_JesSFup->at(boostedjet_en);
        boostjet.JesSFdown= rBoostedJet_JesSFdown->at(boostedjet_en);
        boostjet.JerSF= rBoostedJet_JerSF->at(boostedjet_en);
        boostjet.JerSFup= rBoostedJet_JerSFup->at(boostedjet_en);
        boostjet.JerSFdown= rBoostedJet_JerSFdown->at(boostedjet_en);
        boostjet.tau1= rBoostedJet_tau1->at(boostedjet_en);
        boostjet.tau2= rBoostedJet_tau2->at(boostedjet_en);
        boostjet.tau3= rBoostedJet_tau3->at(boostedjet_en);
        boostjet.softdrop_mass= rBoostedJet_softdrop_mass->at(boostedjet_en);
        boostjet.pruned_mass= rBoostedJet_pruned_mass->at(boostedjet_en);
        if (boostjet.tau1!=0) boostjet.tau21 = boostjet.tau2/boostjet.tau1;
        if (boostjet.tau2!=0) boostjet.tau32 = boostjet.tau3/boostjet.tau2;

        // select boosted Top and W
        boostjet.set_Wp_Top(top_numSoft, top_numLoose, top_numMedium, top_numTight);
        boostjet.set_Wp_W(w_numLoose, w_numTight);
        
        BoostedJet_pt->push_back(boostjet.pt);
        BoostedJet_eta->push_back(rBoostedJet_eta->at(boostedjet_en));
        BoostedJet_phi->push_back(rBoostedJet_phi->at(boostedjet_en));
        BoostedJet_energy->push_back(boostjet.energy);
        BoostedJet_Uncorr_pt->push_back(rBoostedJet_Uncorr_pt->at(boostedjet_en));
        BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->push_back(rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->at(boostedjet_en));
        BoostedJet_pfCombinedMVAV2BJetTags->push_back(rBoostedJet_pfCombinedMVAV2BJetTags->at(boostedjet_en));
        BoostedJet_JesSF->push_back(rBoostedJet_JesSF->at(boostedjet_en));
        BoostedJet_JesSFup->push_back(rBoostedJet_JesSFup->at(boostedjet_en));
        BoostedJet_JesSFdown->push_back(rBoostedJet_JesSFdown->at(boostedjet_en));
        BoostedJet_JerSF->push_back(rBoostedJet_JerSF->at(boostedjet_en));
        BoostedJet_JerSFup->push_back(rBoostedJet_JerSFup->at(boostedjet_en));
        BoostedJet_JerSFdown->push_back(rBoostedJet_JerSFdown->at(boostedjet_en));
        BoostedJet_tau1->push_back(rBoostedJet_tau1->at(boostedjet_en));
        BoostedJet_tau2->push_back(rBoostedJet_tau2->at(boostedjet_en));
        BoostedJet_tau3->push_back(rBoostedJet_tau3->at(boostedjet_en));
        BoostedJet_softdrop_mass->push_back(rBoostedJet_softdrop_mass->at(boostedjet_en));
        BoostedJet_pruned_mass->push_back(rBoostedJet_pruned_mass->at(boostedjet_en));
        //new
        BoostedJet_tau21->push_back(boostjet.tau21);
        BoostedJet_tau32->push_back(boostjet.tau32);
        BoostedJet_wCut->push_back(boostjet.wCut);
        BoostedJet_topCut->push_back(boostjet.topCut);
       
        boostjets->push_back(boostjet);
    }
    Top_numSoft = top_numSoft;
    Top_numLoose = top_numLoose;
    Top_numMedium = top_numMedium;
    Top_numTight = top_numTight;
    W_numLoose = w_numLoose;
    W_numTight = w_numTight;
};

void Lep_sel(){
    sort(leptons->begin(), leptons->end(), byPt); 
    for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
        Lepton_pt->push_back(leptons->at(lep_en).pt);
        Lepton_eta->push_back(leptons->at(lep_en).eta);
        Lepton_phi->push_back(leptons->at(lep_en).phi);
        Lepton_energy->push_back(leptons->at(lep_en).energy);
        Lepton_dxy_pv->push_back(leptons->at(lep_en).dxy_pv);
        Lepton_dz_pv->push_back(leptons->at(lep_en).dz_pv);
        Lepton_IP3Dsig->push_back(leptons->at(lep_en).IP3Dsig);
        Lepton_loose->push_back(leptons->at(lep_en).loose);
        Lepton_miniIsoRel->push_back(leptons->at(lep_en).loose);
        Lepton_charge->push_back(leptons->at(lep_en).charge);
        Lepton_pdgId->push_back(leptons->at(lep_en).pdgId);
        Lepton_cut->push_back(leptons->at(lep_en).cut);
        Lepton_BDT->push_back(leptons->at(lep_en).BDT);
        Lepton_corrpt->push_back(leptons->at(lep_en).corrpt);
        Lepton_FR->push_back(leptons->at(lep_en).FR);
        Lepton_CF->push_back(leptons->at(lep_en).CF);
        Lepton_passConversion->push_back(leptons->at(lep_en).passConversion);
        Lepton_passMuTightCharge->push_back(leptons->at(lep_en).passMuTightCharge);
        Lepton_passEleTightCharge->push_back(leptons->at(lep_en).passEleTightCharge);
        Lepton_passMissHit->push_back(leptons->at(lep_en).passMissHit);
        Lepton_isMatchRightCharge->push_back(leptons->at(lep_en).isMatchRightCharge);
        Lepton_isGlobal->push_back(leptons->at(lep_en).isGlobal);
        Lepton_chi2->push_back(leptons->at(lep_en).chi2);
        Lepton_chi2LocalPosition->push_back(leptons->at(lep_en).chi2LocalPosition);
        Lepton_trkKink->push_back(leptons->at(lep_en).trkKink);
        Lepton_validFraction->push_back(leptons->at(lep_en).validFraction);
        Lepton_segmentCompatibility->push_back(leptons->at(lep_en).segmentCompatibility);
        Lepton_jetptratio->push_back(leptons->at(lep_en).jetptratio);
        Lepton_jetcsv->push_back(leptons->at(lep_en).jetcsv);
        Lepton_jetpt->push_back(leptons->at(lep_en).jetpt);
        Lepton_lepjetchtrks->push_back(leptons->at(lep_en).lepjetchtrks);
        Lepton_miniIsoCh->push_back(leptons->at(lep_en).miniIsoCh);
        Lepton_miniIsoPUsub->push_back(leptons->at(lep_en).miniIsoPUsub);
        Lepton_ptrel->push_back(leptons->at(lep_en).ptrel);
        Lepton_pTErrOVpT_it->push_back(leptons->at(lep_en).pTErrOVpT_it);
        Lepton_px->push_back(leptons->at(lep_en).px);
        Lepton_py->push_back(leptons->at(lep_en).py);
        Lepton_pz->push_back(leptons->at(lep_en).pz);
        Lepton_jetdr->push_back(leptons->at(lep_en).jetdr);
        Lepton_gen_pt->push_back(leptons->at(lep_en).gen_pt);
        Lepton_gen_eta->push_back(leptons->at(lep_en).gen_eta);
        Lepton_gen_phi->push_back(leptons->at(lep_en).gen_phi);
        Lepton_gen_en->push_back(leptons->at(lep_en).gen_en);
        Lepton_gen_pdgId->push_back(leptons->at(lep_en).gen_pdgId);
        Lepton_genMother_pt->push_back(leptons->at(lep_en).genMother_pt);
        Lepton_genMother_eta->push_back(leptons->at(lep_en).genMother_eta);
        Lepton_genMother_phi->push_back(leptons->at(lep_en).genMother_phi);
        Lepton_genMother_en->push_back(leptons->at(lep_en).genMother_en);
        Lepton_genMother_pdgId->push_back(leptons->at(lep_en).genMother_pdgId);
        Lepton_genGrandMother_pt->push_back(leptons->at(lep_en).genGrandMother_pt);
        Lepton_genGrandMother_eta->push_back(leptons->at(lep_en).genGrandMother_eta);
        Lepton_genGrandMother_phi->push_back(leptons->at(lep_en).genGrandMother_phi);
        Lepton_genGrandMother_en->push_back(leptons->at(lep_en).genGrandMother_en);
        Lepton_genGrandMother_pdgId->push_back(leptons->at(lep_en).genGrandMother_pdgId);
        Lepton_gen_isPromptFinalState->push_back(leptons->at(lep_en).gen_isPromptFinalState);
        Lepton_gen_isDirectPromptTauDecayProductFinalState->push_back(leptons->at(lep_en).gen_isDirectPromptTauDecayProductFinalState);
        Lepton_SCeta->push_back(leptons->at(lep_en).SCeta);
        Lepton_mvaValue_HZZ->push_back(leptons->at(lep_en).mvaValue_HZZ);
        Lepton_expectedMissingInnerHits->push_back(leptons->at(lep_en).expectedMissingInnerHits);
        Lepton_full5x5_sigmaIetaIeta->push_back(leptons->at(lep_en).full5x5_sigmaIetaIeta);
        Lepton_hOverE->push_back(leptons->at(lep_en).hOverE);
        Lepton_dEtaIn->push_back(leptons->at(lep_en).dEtaIn);
        Lepton_dPhiIn->push_back(leptons->at(lep_en).dPhiIn);
        Lepton_ooEmooP->push_back(leptons->at(lep_en).ooEmooP);
    }
};


void Event_sel(){
    //Write_variables
    EVENT_event = rEVENT_event;
    EVENT_genWeight = rEVENT_genWeight;
    HiggsDecay = rHiggsDecay;
    PUWeight = rPUWeight;
    Met_type1PF_py = rMet_type1PF_py;
    Met_type1PF_pz = rMet_type1PF_pz;
    Met_type1PF_phi = rMet_type1PF_phi;
    Met_type1PF_sumEt = rMet_type1PF_sumEt;
    Met_type1PF_shiftedPtUp = rMet_type1PF_shiftedPtUp;
    Met_type1PF_shiftedPtDown = rMet_type1PF_shiftedPtDown;
    for(uint gen_en=0; gen_en<rGen_pdg_id->size(); gen_en++){
        Gen_pt->push_back(rGen_pt->at(gen_en));
        Gen_eta->push_back(rGen_eta->at(gen_en));
        Gen_phi->push_back(rGen_phi->at(gen_en));
        Gen_energy->push_back(rGen_energy->at(gen_en));
        Gen_pdg_id->push_back(rGen_pdg_id->at(gen_en));
        Gen_motherpdg_id->push_back(rGen_motherpdg_id->at(gen_en));
        Gen_BmotherIndex->push_back(rGen_BmotherIndex->at(gen_en));
        Gen_numMother->push_back(rGen_numMother->at(gen_en));
        Gen_status->push_back(rGen_status->at(gen_en));
    }
    Gen_type1PF_Met = rGen_type1PF_Met;
    Gen_type1PF_Metpx = rGen_type1PF_Metpx;
    Gen_type1PF_Metpy = rGen_type1PF_Metpy;
    Gen_type1PF_Metpz = rGen_type1PF_Metpz;
    Gen_type1PF_Meteta = rGen_type1PF_Meteta;
    Gen_type1PF_Metphi = rGen_type1PF_Metphi;
    Gen_type1PF_Meten = rGen_type1PF_Meten;
    Met_type1PF_pt = rMet_type1PF_pt;
    Met_type1PF_px = rMet_type1PF_px;
    // Event new info
    //Triggers
    if( rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ==1 || rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ==1 ||rHLT_IsoMu22==1 || rHLT_IsoTkMu22==1 || rHLT_IsoMu22_eta2p1==1 || rHLT_IsoTkMu22_eta2p1 ==1||rHLT_IsoMu24 ==1 || rHLT_IsoTkMu24==1)TTHLep_2Mu=1;
    else TTHLep_2Mu=0;
    if( rHLT_Ele27_WPTight_Gsf==1 || rHLT_Ele25_eta2p1_WPTight_Gsf==1 || rHLT_Ele27_eta2p1_WPLoose_Gsf==1 || rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ==1 )TTHLep_2Ele=1;
    else TTHLep_2Ele=0;
    if( rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL ==1 || rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ==1 || rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ==1 || rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ==1 
    || rHLT_IsoMu22==1 || rHLT_IsoTkMu22==1 || rHLT_IsoMu22_eta2p1==1 || rHLT_IsoTkMu22_eta2p1 ==1||rHLT_IsoMu24 ==1 || rHLT_IsoTkMu24==1 
    || rHLT_Ele27_WPTight_Gsf==1 || rHLT_Ele25_eta2p1_WPTight_Gsf==1 || rHLT_Ele27_eta2p1_WPLoose_Gsf==1
    ) TTHLep_MuEle=1;
    else TTHLep_MuEle=0;
    if(rHLT_DiMu9_Ele9_CaloIdL_TrackIdL==1 || rHLT_Mu8_DiEle12_CaloIdL_TrackIdL ==1 || rHLT_TripleMu_12_10_5==1 || rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL==1 
    || rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ==1 || rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ==1 || rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ==1 || rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL ==1 || rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ==1 || rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ==1 || rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ==1
    || rHLT_IsoMu22==1 || rHLT_IsoTkMu22==1 || rHLT_IsoMu22_eta2p1==1 || rHLT_IsoTkMu22_eta2p1 ==1||rHLT_IsoMu24 ==1 || rHLT_IsoTkMu24==1 || rHLT_Ele27_WPTight_Gsf==1 || rHLT_Ele25_eta2p1_WPTight_Gsf==1 || rHLT_Ele27_eta2p1_WPLoose_Gsf==1
    ) TTHLep_3L4L=1;
    else TTHLep_3L4L=0;
    // ht and St
    TLorentzVector mht_lv(0,0,0,0);
    TLorentzVector mhtT_lv(0,0,0,0);
    for(uint lv_en = 0; lv_en<leptons->size(); lv_en++){
        TLorentzVector lv(0.,0.,0.,0.);
        lv.SetPtEtaPhiE(leptons->at(lv_en).pt,leptons->at(lv_en).eta,leptons->at(lv_en).phi,leptons->at(lv_en).energy);
        mht_lv = mht_lv + lv;
        if(leptons->at(lv_en).cut==3){
            mhtT_lv = mhtT_lv + lv;
        }
    }
    for(uint lv_en = 0; lv_en<taus->size(); lv_en++){
        TLorentzVector lv(0.,0.,0.,0.);
        lv.SetPtEtaPhiE(taus->at(lv_en).pt,taus->at(lv_en).eta,taus->at(lv_en).phi,taus->at(lv_en).energy);
        mht_lv = mht_lv + lv;
        if(taus->at(lv_en).cut==2){
            mhtT_lv = mhtT_lv + lv;
        }
    }
    for(uint lv_en = 0; lv_en<jets->size(); lv_en++){
        TLorentzVector lv(0.,0.,0.,0.);
        lv.SetPtEtaPhiE(jets->at(lv_en).pt,jets->at(lv_en).eta,jets->at(lv_en).phi,jets->at(lv_en).energy);
        mht_lv = mht_lv + lv;
        mhtT_lv = mhtT_lv + lv;
    }
    mht = mht_lv.Pt();
    mhtT = mhtT_lv.Pt();
    mht_met = mht + Met_type1PF_pt;
    mhtT_met = mhtT + Met_type1PF_pt;
    metLD =  0.00397*Met_type1PF_pt+0.00265*mht;  
}


//////
// Variables handling
//////
void rSetBranchAddress(TTree* readingtree, string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //read setbranchaddress
    //Event
    readingtree->SetBranchAddress("EVENT_event",&rEVENT_event,&b_rEVENT_event);
    readingtree->SetBranchAddress("EVENT_genWeight",&rEVENT_genWeight,&b_rEVENT_genWeight);
    readingtree->SetBranchAddress("HiggsDecay",&rHiggsDecay,&b_rHiggsDecay);
    readingtree->SetBranchAddress("PUWeight",&rPUWeight,&b_rPUWeight);
    readingtree->SetBranchAddress("Gen_type1PF_Met",&rGen_type1PF_Met,&b_rGen_type1PF_Met);
    readingtree->SetBranchAddress("Gen_type1PF_Metpx",&rGen_type1PF_Metpx,&b_rGen_type1PF_Metpx);
    readingtree->SetBranchAddress("Gen_type1PF_Metpy",&rGen_type1PF_Metpy,&b_rGen_type1PF_Metpy);
    readingtree->SetBranchAddress("Gen_type1PF_Metpz",&rGen_type1PF_Metpz,&b_rGen_type1PF_Metpz);
    readingtree->SetBranchAddress("Gen_type1PF_Meteta",&rGen_type1PF_Meteta,&b_rGen_type1PF_Meteta);
    readingtree->SetBranchAddress("Gen_type1PF_Metphi",&rGen_type1PF_Metphi,&b_rGen_type1PF_Metphi);
    readingtree->SetBranchAddress("Gen_type1PF_Meten",&rGen_type1PF_Meten,&b_rGen_type1PF_Meten);
    readingtree->SetBranchAddress("Met_type1PF_pt",&rMet_type1PF_pt,&b_rMet_type1PF_pt);
    readingtree->SetBranchAddress("Met_type1PF_px",&rMet_type1PF_px,&b_rMet_type1PF_px);
    readingtree->SetBranchAddress("Met_type1PF_py",&rMet_type1PF_py,&b_rMet_type1PF_py);
    readingtree->SetBranchAddress("Met_type1PF_pz",&rMet_type1PF_pz,&b_rMet_type1PF_pz);
    readingtree->SetBranchAddress("Met_type1PF_phi",&rMet_type1PF_phi,&b_rMet_type1PF_phi);
    readingtree->SetBranchAddress("Met_type1PF_sumEt",&rMet_type1PF_sumEt,&b_rMet_type1PF_sumEt);
    readingtree->SetBranchAddress("Met_type1PF_shiftedPtUp",&rMet_type1PF_shiftedPtUp,&b_rMet_type1PF_shiftedPtUp);
    readingtree->SetBranchAddress("Met_type1PF_shiftedPtDown",&rMet_type1PF_shiftedPtDown,&b_rMet_type1PF_shiftedPtDown);
    readingtree->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL",&rHLT_DiMu9_Ele9_CaloIdL_TrackIdL,&b_rHLT_DiMu9_Ele9_CaloIdL_TrackIdL);
    readingtree->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL",&rHLT_Mu8_DiEle12_CaloIdL_TrackIdL,&b_rHLT_Mu8_DiEle12_CaloIdL_TrackIdL);
    readingtree->SetBranchAddress("HLT_TripleMu_12_10_5",&rHLT_TripleMu_12_10_5,&b_rHLT_TripleMu_12_10_5);
    readingtree->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",&rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL,&b_rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
    readingtree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",&rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL,&b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
    readingtree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",&rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ,&b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
    readingtree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",&rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,&b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    readingtree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",&rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,&b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    readingtree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,&b_rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    readingtree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,&b_rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    readingtree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,&b_rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    readingtree->SetBranchAddress("HLT_IsoMu22",&rHLT_IsoMu22,&b_rHLT_IsoMu22);
    readingtree->SetBranchAddress("HLT_IsoTkMu22",&rHLT_IsoTkMu22,&b_rHLT_IsoTkMu22);
    readingtree->SetBranchAddress("HLT_IsoMu22_eta2p1",&rHLT_IsoMu22_eta2p1,&b_rHLT_IsoMu22_eta2p1);
    readingtree->SetBranchAddress("HLT_IsoTkMu22_eta2p1",&rHLT_IsoTkMu22_eta2p1,&b_rHLT_IsoTkMu22_eta2p1);
    readingtree->SetBranchAddress("HLT_IsoMu24",&rHLT_IsoMu24,&b_rHLT_IsoMu24);
    readingtree->SetBranchAddress("HLT_IsoTkMu24",&rHLT_IsoTkMu24,&b_rHLT_IsoTkMu24);
    readingtree->SetBranchAddress("HLT_Ele27_WPTight_Gsf",&rHLT_Ele27_WPTight_Gsf,&b_rHLT_Ele27_WPTight_Gsf);
    readingtree->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf",&rHLT_Ele25_eta2p1_WPTight_Gsf,&b_rHLT_Ele25_eta2p1_WPTight_Gsf);
    readingtree->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf",&rHLT_Ele27_eta2p1_WPLoose_Gsf,&b_rHLT_Ele27_eta2p1_WPLoose_Gsf);
    //Muon
    readingtree->SetBranchAddress("Muon_pt",&rMuon_pt,&b_rMuon_pt);
    readingtree->SetBranchAddress("Muon_eta",&rMuon_eta,&b_rMuon_eta);
    readingtree->SetBranchAddress("Muon_phi",&rMuon_phi,&b_rMuon_phi);
    readingtree->SetBranchAddress("Muon_energy",&rMuon_energy,&b_rMuon_energy);
    readingtree->SetBranchAddress("Muon_dxy_pv",&rMuon_dxy_pv,&b_rMuon_dxy_pv);
    readingtree->SetBranchAddress("Muon_dz_pv",&rMuon_dz_pv,&b_rMuon_dz_pv);
    readingtree->SetBranchAddress("Muon_IP3Dsig_it",&rMuon_IP3Dsig_it,&b_rMuon_IP3Dsig_it);
    readingtree->SetBranchAddress("Muon_loose",&rMuon_loose,&b_rMuon_loose);
    readingtree->SetBranchAddress("Muon_miniIsoRel",&rMuon_miniIsoRel,&b_rMuon_miniIsoRel);
    readingtree->SetBranchAddress("Muon_charge",&rMuon_charge,&b_rMuon_charge);
    readingtree->SetBranchAddress("Muon_pdgId",&rMuon_pdgId,&b_rMuon_pdgId);
    readingtree->SetBranchAddress("Muon_isGlobal",&rMuon_isGlobal,&b_rMuon_isGlobal);
    readingtree->SetBranchAddress("Muon_chi2",&rMuon_chi2,&b_rMuon_chi2);
    readingtree->SetBranchAddress("Muon_chi2LocalPosition",&rMuon_chi2LocalPosition,&b_rMuon_chi2LocalPosition);
    readingtree->SetBranchAddress("Muon_trkKink",&rMuon_trkKink,&b_rMuon_trkKink);
    readingtree->SetBranchAddress("Muon_validFraction",&rMuon_validFraction,&b_rMuon_validFraction);
    readingtree->SetBranchAddress("Muon_segmentCompatibility",&rMuon_segmentCompatibility,&b_rMuon_segmentCompatibility);
    readingtree->SetBranchAddress("Muon_jetptratio",&rMuon_jetptratio,&b_rMuon_jetptratio);
    readingtree->SetBranchAddress("Muon_jetcsv",&rMuon_jetcsv,&b_rMuon_jetcsv);
    readingtree->SetBranchAddress("Muon_jetpt",&rMuon_jetpt,&b_rMuon_jetpt);
    readingtree->SetBranchAddress("Muon_lepjetchtrks",&rMuon_lepjetchtrks,&b_rMuon_lepjetchtrks);
    readingtree->SetBranchAddress("Muon_miniIsoCh",&rMuon_miniIsoCh,&b_rMuon_miniIsoCh);
    readingtree->SetBranchAddress("Muon_miniIsoPUsub",&rMuon_miniIsoPUsub,&b_rMuon_miniIsoPUsub);
    readingtree->SetBranchAddress("Muon_ptrel",&rMuon_ptrel,&b_rMuon_ptrel);
    readingtree->SetBranchAddress("Muon_pTErrOVpT_it",&rMuon_pTErrOVpT_it,&b_rMuon_pTErrOVpT_it);
    readingtree->SetBranchAddress("Muon_px",&rMuon_px,&b_rMuon_px);
    readingtree->SetBranchAddress("Muon_py",&rMuon_py,&b_rMuon_py);
    readingtree->SetBranchAddress("Muon_pz",&rMuon_pz,&b_rMuon_pz);
    readingtree->SetBranchAddress("Muon_jetdr",&rMuon_jetdr,&b_rMuon_jetdr);
    readingtree->SetBranchAddress("Muon_gen_pt",&rMuon_gen_pt,&b_rMuon_gen_pt);
    readingtree->SetBranchAddress("Muon_gen_eta",&rMuon_gen_eta,&b_rMuon_gen_eta);
    readingtree->SetBranchAddress("Muon_gen_phi",&rMuon_gen_phi,&b_rMuon_gen_phi);
    readingtree->SetBranchAddress("Muon_gen_en",&rMuon_gen_en,&b_rMuon_gen_en);
    readingtree->SetBranchAddress("Muon_gen_pdgId",&rMuon_gen_pdgId,&b_rMuon_gen_pdgId);
    readingtree->SetBranchAddress("Muon_genMother_pt",&rMuon_genMother_pt,&b_rMuon_genMother_pt);
    readingtree->SetBranchAddress("Muon_genMother_eta",&rMuon_genMother_eta,&b_rMuon_genMother_eta);
    readingtree->SetBranchAddress("Muon_genMother_phi",&rMuon_genMother_phi,&b_rMuon_genMother_phi);
    readingtree->SetBranchAddress("Muon_genMother_en",&rMuon_genMother_en,&b_rMuon_genMother_en);
    readingtree->SetBranchAddress("Muon_genMother_pdgId",&rMuon_genMother_pdgId,&b_rMuon_genMother_pdgId);
    readingtree->SetBranchAddress("Muon_genGrandMother_pt",&rMuon_genGrandMother_pt,&b_rMuon_genGrandMother_pt);
    readingtree->SetBranchAddress("Muon_genGrandMother_eta",&rMuon_genGrandMother_eta,&b_rMuon_genGrandMother_eta);
    readingtree->SetBranchAddress("Muon_genGrandMother_phi",&rMuon_genGrandMother_phi,&b_rMuon_genGrandMother_phi);
    readingtree->SetBranchAddress("Muon_genGrandMother_en",&rMuon_genGrandMother_en,&b_rMuon_genGrandMother_en);
    readingtree->SetBranchAddress("Muon_genGrandMother_pdgId",&rMuon_genGrandMother_pdgId,&b_rMuon_genGrandMother_pdgId);
    readingtree->SetBranchAddress("Muon_gen_isPromptFinalState",&rMuon_gen_isPromptFinalState,&b_rMuon_gen_isPromptFinalState);
    readingtree->SetBranchAddress("Muon_gen_isDirectPromptTauDecayProductFinalState",&rMuon_gen_isDirectPromptTauDecayProductFinalState,&b_rMuon_gen_isDirectPromptTauDecayProductFinalState);
    //Electron
    readingtree->SetBranchAddress("patElectron_pt",&rpatElectron_pt,&b_rpatElectron_pt);
    readingtree->SetBranchAddress("patElectron_eta",&rpatElectron_eta,&b_rpatElectron_eta);
    readingtree->SetBranchAddress("patElectron_phi",&rpatElectron_phi,&b_rpatElectron_phi);
    readingtree->SetBranchAddress("patElectron_energy",&rpatElectron_energy,&b_rpatElectron_energy);
    readingtree->SetBranchAddress("patElectron_IP3Dsig",&rpatElectron_IP3Dsig,&b_rpatElectron_IP3Dsig);
    readingtree->SetBranchAddress("patElectron_miniIsoRel",&rpatElectron_miniIsoRel,&b_rpatElectron_miniIsoRel);
    readingtree->SetBranchAddress("patElectron_charge",&rpatElectron_charge,&b_rpatElectron_charge);
    readingtree->SetBranchAddress("patElectron_pdgId",&rpatElectron_pdgId,&b_rpatElectron_pdgId);
    readingtree->SetBranchAddress("patElectron_gsfTrack_dxy_pv",&rpatElectron_gsfTrack_dxy_pv,&b_rpatElectron_gsfTrack_dxy_pv);
    readingtree->SetBranchAddress("patElectron_gsfTrack_dz_pv",&rpatElectron_gsfTrack_dz_pv,&b_rpatElectron_gsfTrack_dz_pv);
    readingtree->SetBranchAddress("patElectron_jetptratio",&rpatElectron_jetptratio,&b_rpatElectron_jetptratio);
    readingtree->SetBranchAddress("patElectron_jetcsv",&rpatElectron_jetcsv,&b_rpatElectron_jetcsv);
    readingtree->SetBranchAddress("patElectron_jetpt",&rpatElectron_jetpt,&b_rpatElectron_jetpt);
    readingtree->SetBranchAddress("patElectron_lepjetchtrks",&rpatElectron_lepjetchtrks,&b_rpatElectron_lepjetchtrks);
    readingtree->SetBranchAddress("patElectron_miniIsoCh",&rpatElectron_miniIsoCh,&b_rpatElectron_miniIsoCh);
    readingtree->SetBranchAddress("patElectron_miniIsoPUsub",&rpatElectron_miniIsoPUsub,&b_rpatElectron_miniIsoPUsub);
    readingtree->SetBranchAddress("patElectron_ptrel",&rpatElectron_ptrel,&b_rpatElectron_ptrel);
    readingtree->SetBranchAddress("patElectron_px",&rpatElectron_px,&b_rpatElectron_px);
    readingtree->SetBranchAddress("patElectron_py",&rpatElectron_py,&b_rpatElectron_py);
    readingtree->SetBranchAddress("patElectron_pz",&rpatElectron_pz,&b_rpatElectron_pz);
    readingtree->SetBranchAddress("patElectron_jetdr",&rpatElectron_jetdr,&b_rpatElectron_jetdr);
    readingtree->SetBranchAddress("patElectron_gen_pt",&rpatElectron_gen_pt,&b_rpatElectron_gen_pt);
    readingtree->SetBranchAddress("patElectron_gen_eta",&rpatElectron_gen_eta,&b_rpatElectron_gen_eta);
    readingtree->SetBranchAddress("patElectron_gen_phi",&rpatElectron_gen_phi,&b_rpatElectron_gen_phi);
    readingtree->SetBranchAddress("patElectron_gen_en",&rpatElectron_gen_en,&b_rpatElectron_gen_en);
    readingtree->SetBranchAddress("patElectron_gen_pdgId",&rpatElectron_gen_pdgId,&b_rpatElectron_gen_pdgId);
    readingtree->SetBranchAddress("patElectron_genMother_pt",&rpatElectron_genMother_pt,&b_rpatElectron_genMother_pt);
    readingtree->SetBranchAddress("patElectron_genMother_eta",&rpatElectron_genMother_eta,&b_rpatElectron_genMother_eta);
    readingtree->SetBranchAddress("patElectron_genMother_phi",&rpatElectron_genMother_phi,&b_rpatElectron_genMother_phi);
    readingtree->SetBranchAddress("patElectron_genMother_en",&rpatElectron_genMother_en,&b_rpatElectron_genMother_en);
    readingtree->SetBranchAddress("patElectron_genMother_pdgId",&rpatElectron_genMother_pdgId,&b_rpatElectron_genMother_pdgId);
    readingtree->SetBranchAddress("patElectron_genGrandMother_pt",&rpatElectron_genGrandMother_pt,&b_rpatElectron_genGrandMother_pt);
    readingtree->SetBranchAddress("patElectron_genGrandMother_eta",&rpatElectron_genGrandMother_eta,&b_rpatElectron_genGrandMother_eta);
    readingtree->SetBranchAddress("patElectron_genGrandMother_phi",&rpatElectron_genGrandMother_phi,&b_rpatElectron_genGrandMother_phi);
    readingtree->SetBranchAddress("patElectron_genGrandMother_en",&rpatElectron_genGrandMother_en,&b_rpatElectron_genGrandMother_en);
    readingtree->SetBranchAddress("patElectron_genGrandMother_pdgId",&rpatElectron_genGrandMother_pdgId,&b_rpatElectron_genGrandMother_pdgId);
    readingtree->SetBranchAddress("patElectron_gen_isPromptFinalState",&rpatElectron_gen_isPromptFinalState,&b_rpatElectron_gen_isPromptFinalState);
    readingtree->SetBranchAddress("patElectron_gen_isDirectPromptTauDecayProductFinalState",&rpatElectron_gen_isDirectPromptTauDecayProductFinalState,&b_rpatElectron_gen_isDirectPromptTauDecayProductFinalState);
    readingtree->SetBranchAddress("patElectron_SCeta",&rpatElectron_SCeta,&b_rpatElectron_SCeta);
    readingtree->SetBranchAddress("patElectron_mvaValue_HZZ",&rpatElectron_mvaValue_HZZ,&b_rpatElectron_mvaValue_HZZ);
    readingtree->SetBranchAddress("patElectron_expectedMissingInnerHits",&rpatElectron_expectedMissingInnerHits,&b_rpatElectron_expectedMissingInnerHits);
    readingtree->SetBranchAddress("patElectron_full5x5_sigmaIetaIeta",&rpatElectron_full5x5_sigmaIetaIeta,&b_rpatElectron_full5x5_sigmaIetaIeta);
    readingtree->SetBranchAddress("patElectron_hOverE",&rpatElectron_hOverE,&b_rpatElectron_hOverE);
    readingtree->SetBranchAddress("patElectron_dEtaIn",&rpatElectron_dEtaIn,&b_rpatElectron_dEtaIn);
    readingtree->SetBranchAddress("patElectron_dPhiIn",&rpatElectron_dPhiIn,&b_rpatElectron_dPhiIn);
    readingtree->SetBranchAddress("patElectron_ooEmooP",&rpatElectron_ooEmooP,&b_rpatElectron_ooEmooP);
    readingtree->SetBranchAddress("patElectron_isGsfCtfScPixChargeConsistent",&rpatElectron_isGsfCtfScPixChargeConsistent,&b_rpatElectron_isGsfCtfScPixChargeConsistent);
    readingtree->SetBranchAddress("patElectron_isGsfScPixChargeConsistent",&rpatElectron_isGsfScPixChargeConsistent,&b_rpatElectron_isGsfScPixChargeConsistent);
    readingtree->SetBranchAddress("patElectron_passConversionVeto",&rpatElectron_passConversionVeto,&b_rpatElectron_passConversionVeto);
    //Tau
    readingtree->SetBranchAddress("Tau_pt",&rTau_pt,&b_rTau_pt);
    readingtree->SetBranchAddress("Tau_eta",&rTau_eta,&b_rTau_eta);
    readingtree->SetBranchAddress("Tau_phi",&rTau_phi,&b_rTau_phi);
    readingtree->SetBranchAddress("Tau_energy",&rTau_energy,&b_rTau_energy);
    readingtree->SetBranchAddress("Tau_charge",&rTau_charge,&b_rTau_charge);
    readingtree->SetBranchAddress("Tau_packedLeadTauCand_dz",&rTau_packedLeadTauCand_dz,&b_rTau_packedLeadTauCand_dz);
    readingtree->SetBranchAddress("Tau_packedLeadTauCand_dxy",&rTau_packedLeadTauCand_dxy,&b_rTau_packedLeadTauCand_dxy);
    readingtree->SetBranchAddress("Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT,&b_rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    readingtree->SetBranchAddress("Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT,&b_rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    readingtree->SetBranchAddress("Tau_decayModeFinding",&rTau_decayModeFinding,&b_rTau_decayModeFinding);
    // Jets
    readingtree->SetBranchAddress("Jet_pt",&rJet_pt,&b_rJet_pt);
    readingtree->SetBranchAddress("Jet_eta",&rJet_eta,&b_rJet_eta);
    readingtree->SetBranchAddress("Jet_phi",&rJet_phi,&b_rJet_phi);
    readingtree->SetBranchAddress("Jet_energy",&rJet_energy,&b_rJet_energy);
    readingtree->SetBranchAddress("Jet_genMother_pt",&rJet_genMother_pt,&b_rJet_genMother_pt);
    readingtree->SetBranchAddress("Jet_genMother_eta",&rJet_genMother_eta,&b_rJet_genMother_eta);
    readingtree->SetBranchAddress("Jet_genMother_phi",&rJet_genMother_phi,&b_rJet_genMother_phi);
    readingtree->SetBranchAddress("Jet_genMother_en",&rJet_genMother_en,&b_rJet_genMother_en);
    readingtree->SetBranchAddress("Jet_genMother_pdgId",&rJet_genMother_pdgId,&b_rJet_genMother_pdgId);
    readingtree->SetBranchAddress("Jet_genGrandMother_pt",&rJet_genGrandMother_pt,&b_rJet_genGrandMother_pt);
    readingtree->SetBranchAddress("Jet_genGrandMother_eta",&rJet_genGrandMother_eta,&b_rJet_genGrandMother_eta);
    readingtree->SetBranchAddress("Jet_genGrandMother_phi",&rJet_genGrandMother_phi,&b_rJet_genGrandMother_phi);
    readingtree->SetBranchAddress("Jet_genGrandMother_en",&rJet_genGrandMother_en,&b_rJet_genGrandMother_en);
    readingtree->SetBranchAddress("Jet_genGrandMother_pdgId",&rJet_genGrandMother_pdgId,&b_rJet_genGrandMother_pdgId);
    readingtree->SetBranchAddress("Jet_Uncorr_pt",&rJet_Uncorr_pt,&b_rJet_Uncorr_pt);
    readingtree->SetBranchAddress("Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags",&rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags,&b_rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
    readingtree->SetBranchAddress("Jet_pfCombinedMVAV2BJetTags",&rJet_pfCombinedMVAV2BJetTags,&b_rJet_pfCombinedMVAV2BJetTags);
    readingtree->SetBranchAddress("Jet_qg",&rJet_qg,&b_rJet_qg);
    readingtree->SetBranchAddress("Jet_axis2",&rJet_axis2,&b_rJet_axis2);
    readingtree->SetBranchAddress("Jet_ptD",&rJet_ptD,&b_rJet_ptD);
    readingtree->SetBranchAddress("Jet_mult",&rJet_mult,&b_rJet_mult);
    readingtree->SetBranchAddress("Jet_partonFlavour",&rJet_partonFlavour,&b_rJet_partonFlavour);
    readingtree->SetBranchAddress("Jet_hadronFlavour",&rJet_hadronFlavour,&b_rJet_hadronFlavour);
    readingtree->SetBranchAddress("Jet_genpt",&rJet_genpt,&b_rJet_genpt);
    readingtree->SetBranchAddress("Jet_geneta",&rJet_geneta,&b_rJet_geneta);
    readingtree->SetBranchAddress("Jet_genphi",&rJet_genphi,&b_rJet_genphi);
    readingtree->SetBranchAddress("Jet_genenergy",&rJet_genenergy,&b_rJet_genenergy);
    readingtree->SetBranchAddress("Jet_JesSF",&rJet_JesSF,&b_rJet_JesSF);
    readingtree->SetBranchAddress("Jet_JesSFup",&rJet_JesSFup,&b_rJet_JesSFup);
    readingtree->SetBranchAddress("Jet_JesSFdown",&rJet_JesSFdown,&b_rJet_JesSFdown);
    readingtree->SetBranchAddress("Jet_JerSF",&rJet_JerSF,&b_rJet_JerSF);
    readingtree->SetBranchAddress("Jet_JerSFup",&rJet_JerSFup,&b_rJet_JerSFup);
    readingtree->SetBranchAddress("Jet_JerSFdown",&rJet_JerSFdown,&b_rJet_JerSFdown);
    readingtree->SetBranchAddress("Jet_neutralHadEnergyFraction",&rJet_neutralHadEnergyFraction,&b_rJet_neutralHadEnergyFraction);
    readingtree->SetBranchAddress("Jet_neutralEmEnergyFraction",&rJet_neutralEmEnergyFraction,&b_rJet_neutralEmEnergyFraction);
    readingtree->SetBranchAddress("Jet_chargedMultiplicity",&rJet_chargedMultiplicity,&b_rJet_chargedMultiplicity);
    readingtree->SetBranchAddress("Jet_numberOfConstituents",&rJet_numberOfConstituents,&b_rJet_numberOfConstituents);
    readingtree->SetBranchAddress("Jet_chargedHadronEnergyFraction",&rJet_chargedHadronEnergyFraction,&b_rJet_chargedHadronEnergyFraction);
    readingtree->SetBranchAddress("Jet_chargedEmEnergyFraction",&rJet_chargedEmEnergyFraction,&b_rJet_chargedEmEnergyFraction);
    readingtree->SetBranchAddress("Jet_btag_sf",&rJet_btag_sf,&b_rJet_btag_sf);
    readingtree->SetBranchAddress("Jet_btag_jesup",&rJet_btag_jesup,&b_rJet_btag_jesup);
    readingtree->SetBranchAddress("Jet_btag_jesdown",&rJet_btag_jesdown,&b_rJet_btag_jesdown);
    readingtree->SetBranchAddress("Jet_btag_hfup",&rJet_btag_hfup,&b_rJet_btag_hfup);
    readingtree->SetBranchAddress("Jet_btag_hfdown",&rJet_btag_hfdown,&b_rJet_btag_hfdown);
    readingtree->SetBranchAddress("Jet_btag_hfstat1up",&rJet_btag_hfstat1up,&b_rJet_btag_hfstat1up);
    readingtree->SetBranchAddress("Jet_btag_hfstat1down",&rJet_btag_hfstat1down,&b_rJet_btag_hfstat1down);
    readingtree->SetBranchAddress("Jet_btag_hfstat2up",&rJet_btag_hfstat2up,&b_rJet_btag_hfstat2up);
    readingtree->SetBranchAddress("Jet_btag_hfstat2down",&rJet_btag_hfstat2down,&b_rJet_btag_hfstat2down);
    readingtree->SetBranchAddress("Jet_btag_lfup",&rJet_btag_lfup,&b_rJet_btag_lfup);
    readingtree->SetBranchAddress("Jet_btag_lfdown",&rJet_btag_lfdown,&b_rJet_btag_lfdown);
    readingtree->SetBranchAddress("Jet_btag_lfstat1up",&rJet_btag_lfstat1up,&b_rJet_btag_lfstat1up);
    readingtree->SetBranchAddress("Jet_btag_lfstat1down",&rJet_btag_lfstat1down,&b_rJet_btag_lfstat1down);
    readingtree->SetBranchAddress("Jet_btag_lfstat2up",&rJet_btag_lfstat2up,&b_rJet_btag_lfstat2up);
    readingtree->SetBranchAddress("Jet_btag_lfstat2down",&rJet_btag_lfstat2down,&b_rJet_btag_lfstat2down);
    readingtree->SetBranchAddress("Jet_btag_cerr1up",&rJet_btag_cerr1up,&b_rJet_btag_cerr1up);
    readingtree->SetBranchAddress("Jet_btag_cerr1down",&rJet_btag_cerr1down,&b_rJet_btag_cerr1down);
    readingtree->SetBranchAddress("Jet_btag_cerr2up",&rJet_btag_cerr2up,&b_rJet_btag_cerr2up);
    readingtree->SetBranchAddress("Jet_btag_cerr2down",&rJet_btag_cerr2down,&b_rJet_btag_cerr2down);
    //BoostJet
    readingtree->SetBranchAddress("BoostedJet_pt",&rBoostedJet_pt,&b_rBoostedJet_pt);
    readingtree->SetBranchAddress("BoostedJet_eta",&rBoostedJet_eta,&b_rBoostedJet_eta);
    readingtree->SetBranchAddress("BoostedJet_phi",&rBoostedJet_phi,&b_rBoostedJet_phi);
    readingtree->SetBranchAddress("BoostedJet_energy",&rBoostedJet_energy,&b_rBoostedJet_energy);
    readingtree->SetBranchAddress("BoostedJet_Uncorr_pt",&rBoostedJet_Uncorr_pt,&b_rBoostedJet_Uncorr_pt);
    readingtree->SetBranchAddress("BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags",&rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags,&b_rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
    readingtree->SetBranchAddress("BoostedJet_pfCombinedMVAV2BJetTags",&rBoostedJet_pfCombinedMVAV2BJetTags,&b_rBoostedJet_pfCombinedMVAV2BJetTags);
    readingtree->SetBranchAddress("BoostedJet_JesSF",&rBoostedJet_JesSF,&b_rBoostedJet_JesSF);
    readingtree->SetBranchAddress("BoostedJet_JesSFup",&rBoostedJet_JesSFup,&b_rBoostedJet_JesSFup);
    readingtree->SetBranchAddress("BoostedJet_JesSFdown",&rBoostedJet_JesSFdown,&b_rBoostedJet_JesSFdown);
    readingtree->SetBranchAddress("BoostedJet_JerSF",&rBoostedJet_JerSF,&b_rBoostedJet_JerSF);
    readingtree->SetBranchAddress("BoostedJet_JerSFup",&rBoostedJet_JerSFup,&b_rBoostedJet_JerSFup);
    readingtree->SetBranchAddress("BoostedJet_JerSFdown",&rBoostedJet_JerSFdown,&b_rBoostedJet_JerSFdown);
    readingtree->SetBranchAddress("BoostedJet_tau1",&rBoostedJet_tau1,&b_rBoostedJet_tau1);
    readingtree->SetBranchAddress("BoostedJet_tau2",&rBoostedJet_tau2,&b_rBoostedJet_tau2);
    readingtree->SetBranchAddress("BoostedJet_tau3",&rBoostedJet_tau3,&b_rBoostedJet_tau3);
    readingtree->SetBranchAddress("BoostedJet_softdrop_mass",&rBoostedJet_softdrop_mass,&b_rBoostedJet_softdrop_mass);
    readingtree->SetBranchAddress("BoostedJet_pruned_mass",&rBoostedJet_pruned_mass,&b_rBoostedJet_pruned_mass);
    //Gen
    readingtree->SetBranchAddress("Gen_pdg_id",&rGen_pdg_id,&b_rGen_pdg_id);
    readingtree->SetBranchAddress("Gen_pt",&rGen_pt,&b_rGen_pt);
    readingtree->SetBranchAddress("Gen_eta",&rGen_eta,&b_rGen_eta);
    readingtree->SetBranchAddress("Gen_phi",&rGen_phi,&b_rGen_phi);
    readingtree->SetBranchAddress("Gen_energy",&rGen_energy,&b_rGen_energy);
    readingtree->SetBranchAddress("Gen_motherpdg_id",&rGen_motherpdg_id,&b_rGen_motherpdg_id);
    readingtree->SetBranchAddress("Gen_BmotherIndex",&rGen_BmotherIndex,&b_rGen_BmotherIndex);
    readingtree->SetBranchAddress("Gen_numMother",&rGen_numMother,&b_rGen_numMother);
    readingtree->SetBranchAddress("Gen_status",&rGen_status,&b_rGen_status);
};


void wSetBranchAddress(TTree* newtree, string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //Write setbranchaddress
    newtree->Branch("EVENT_event",&EVENT_event);
    newtree->Branch("EVENT_genWeight",&EVENT_genWeight);
    newtree->Branch("HiggsDecay",&HiggsDecay);
    newtree->Branch("PUWeight",&PUWeight);
    newtree->Branch("Gen_type1PF_Met",&Gen_type1PF_Met);
    newtree->Branch("Gen_type1PF_Metpx",&Gen_type1PF_Metpx);
    newtree->Branch("Gen_type1PF_Metpy",&Gen_type1PF_Metpy);
    newtree->Branch("Gen_type1PF_Metpz",&Gen_type1PF_Metpz);
    newtree->Branch("Gen_type1PF_Meteta",&Gen_type1PF_Meteta);
    newtree->Branch("Gen_type1PF_Metphi",&Gen_type1PF_Metphi);
    newtree->Branch("Gen_type1PF_Meten",&Gen_type1PF_Meten);
    newtree->Branch("Met_type1PF_pt",&Met_type1PF_pt);
    newtree->Branch("Met_type1PF_px",&Met_type1PF_px);
    newtree->Branch("Met_type1PF_py",&Met_type1PF_py);
    newtree->Branch("Met_type1PF_pz",&Met_type1PF_pz);
    newtree->Branch("Met_type1PF_phi",&Met_type1PF_phi);
    newtree->Branch("Met_type1PF_sumEt",&Met_type1PF_sumEt);
    newtree->Branch("Met_type1PF_shiftedPtUp",&Met_type1PF_shiftedPtUp);
    newtree->Branch("Met_type1PF_shiftedPtDown",&Met_type1PF_shiftedPtDown);
    newtree->Branch("Muon_numLoose",&Muon_numLoose);
    newtree->Branch("Muon_numFake",&Muon_numFake);
    newtree->Branch("Muon_numTight",&Muon_numTight);
    newtree->Branch("patElectron_numLoose",&patElectron_numLoose);
    newtree->Branch("patElectron_numFake",&patElectron_numFake);
    newtree->Branch("patElectron_numTight",&patElectron_numTight);
    newtree->Branch("Tau_numLoose",&Tau_numLoose);
    newtree->Branch("Tau_numMedium",&Tau_numMedium);
    newtree->Branch("Jet_numLoose",&Jet_numLoose);
    newtree->Branch("Jet_numbLoose",&Jet_numbLoose);
    newtree->Branch("Jet_numbMedium",&Jet_numbMedium);
    newtree->Branch("Jet_numbTight",&Jet_numbTight);
    newtree->Branch("BWeight",&BWeight);
    newtree->Branch("BWeightLFup",&BWeightLFup);
    newtree->Branch("BWeightLFdown",&BWeightLFdown);
    newtree->Branch("BWeightHFup",&BWeightHFup);
    newtree->Branch("BWeightHFdown",&BWeightHFdown);
    newtree->Branch("BWeightJESup",&BWeightJESup);
    newtree->Branch("BWeightJESdown",&BWeightJESdown);
    newtree->Branch("BWeightLFStats1up",&BWeightLFStats1up);
    newtree->Branch("BWeightLFStats1down",&BWeightLFStats1down);
    newtree->Branch("BWeightLFStats2up",&BWeightLFStats2up);
    newtree->Branch("BWeightLFStats2down",&BWeightLFStats2down);
    newtree->Branch("BWeightHFStats1up",&BWeightHFStats1up);
    newtree->Branch("BWeightHFStats1down",&BWeightHFStats1down);
    newtree->Branch("BWeightHFStats2up",&BWeightHFStats2up);
    newtree->Branch("BWeightHFStats2down",&BWeightHFStats2down);
    newtree->Branch("BWeightCErr1up",&BWeightCErr1up);
    newtree->Branch("BWeightCErr1down",&BWeightCErr1down);
    newtree->Branch("BWeightCErr2up",&BWeightCErr2up);
    newtree->Branch("BWeightCErr2down",&BWeightCErr2down);
    newtree->Branch("Top_numSoft",&Top_numSoft);
    newtree->Branch("Top_numLoose",&Top_numLoose);
    newtree->Branch("Top_numMedium",&Top_numMedium);
    newtree->Branch("Top_numTight",&Top_numTight);
    newtree->Branch("W_numLoose",&W_numLoose);
    newtree->Branch("W_numTight",&W_numTight);
    newtree->Branch("TTHLep_2Mu",&TTHLep_2Mu);
    newtree->Branch("TTHLep_2Ele",&TTHLep_2Ele);
    newtree->Branch("TTHLep_MuEle",&TTHLep_MuEle);
    newtree->Branch("TTHLep_3L4L",&TTHLep_3L4L);
    newtree->Branch("metLD",&metLD);
    newtree->Branch("mhtT_met",&mhtT_met);
    newtree->Branch("mht_met",&mht_met);
    newtree->Branch("mhtT",&mhtT);
    newtree->Branch("mht",&mht);
    //Lepton
    newtree->Branch("Lepton_pt",&Lepton_pt);
    newtree->Branch("Lepton_eta",&Lepton_eta);
    newtree->Branch("Lepton_phi",&Lepton_phi);
    newtree->Branch("Lepton_energy",&Lepton_energy);
    newtree->Branch("Lepton_dxy_pv",&Lepton_dxy_pv);
    newtree->Branch("Lepton_dz_pv",&Lepton_dz_pv);
    newtree->Branch("Lepton_IP3Dsig",&Lepton_IP3Dsig);
    newtree->Branch("Lepton_loose",&Lepton_loose);
    newtree->Branch("Lepton_miniIsoRel",&Lepton_miniIsoRel);
    newtree->Branch("Lepton_charge",&Lepton_charge);
    newtree->Branch("Lepton_pdgId",&Lepton_pdgId);
    newtree->Branch("Lepton_isGlobal",&Lepton_isGlobal);
    newtree->Branch("Lepton_chi2",&Lepton_chi2);
    newtree->Branch("Lepton_chi2LocalPosition",&Lepton_chi2LocalPosition);
    newtree->Branch("Lepton_trkKink",&Lepton_trkKink);
    newtree->Branch("Lepton_validFraction",&Lepton_validFraction);
    newtree->Branch("Lepton_segmentCompatibility",&Lepton_segmentCompatibility);
    newtree->Branch("Lepton_jetptratio",&Lepton_jetptratio);
    newtree->Branch("Lepton_jetcsv",&Lepton_jetcsv);
    newtree->Branch("Lepton_jetpt",&Lepton_jetpt);
    newtree->Branch("Lepton_lepjetchtrks",&Lepton_lepjetchtrks);
    newtree->Branch("Lepton_miniIsoCh",&Lepton_miniIsoCh);
    newtree->Branch("Lepton_miniIsoPUsub",&Lepton_miniIsoPUsub);
    newtree->Branch("Lepton_ptrel",&Lepton_ptrel);
    newtree->Branch("Lepton_pTErrOVpT_it",&Lepton_pTErrOVpT_it);
    newtree->Branch("Lepton_px",&Lepton_px);
    newtree->Branch("Lepton_py",&Lepton_py);
    newtree->Branch("Lepton_pz",&Lepton_pz);
    newtree->Branch("Lepton_jetdr",&Lepton_jetdr);
    newtree->Branch("Lepton_SCeta",&Lepton_SCeta);
    newtree->Branch("Lepton_mvaValue_HZZ",&Lepton_mvaValue_HZZ);
    newtree->Branch("Lepton_expectedMissingInnerHits",&Lepton_expectedMissingInnerHits);
    newtree->Branch("Lepton_full5x5_sigmaIetaIeta",&Lepton_full5x5_sigmaIetaIeta);
    newtree->Branch("Lepton_hOverE",&Lepton_hOverE);
    newtree->Branch("Lepton_dEtaIn",&Lepton_dEtaIn);
    newtree->Branch("Lepton_dPhiIn",&Lepton_dPhiIn);
    newtree->Branch("Lepton_ooEmooP",&Lepton_ooEmooP);
    // gen variables
    newtree->Branch("Lepton_gen_pt",&Lepton_gen_pt);
    newtree->Branch("Lepton_gen_eta",&Lepton_gen_eta);
    newtree->Branch("Lepton_gen_phi",&Lepton_gen_phi);
    newtree->Branch("Lepton_gen_en",&Lepton_gen_en);
    newtree->Branch("Lepton_gen_pdgId",&Lepton_gen_pdgId);
    newtree->Branch("Lepton_genMother_pt",&Lepton_genMother_pt);
    newtree->Branch("Lepton_genMother_eta",&Lepton_genMother_eta);
    newtree->Branch("Lepton_genMother_phi",&Lepton_genMother_phi);
    newtree->Branch("Lepton_genMother_en",&Lepton_genMother_en);
    newtree->Branch("Lepton_genMother_pdgId",&Lepton_genMother_pdgId);
    newtree->Branch("Lepton_genGrandMother_pt",&Lepton_genGrandMother_pt);
    newtree->Branch("Lepton_genGrandMother_eta",&Lepton_genGrandMother_eta);
    newtree->Branch("Lepton_genGrandMother_phi",&Lepton_genGrandMother_phi);
    newtree->Branch("Lepton_genGrandMother_en",&Lepton_genGrandMother_en);
    newtree->Branch("Lepton_genGrandMother_pdgId",&Lepton_genGrandMother_pdgId);
    newtree->Branch("Lepton_gen_isPromptFinalState",&Lepton_gen_isPromptFinalState);
    newtree->Branch("Lepton_gen_isDirectPromptTauDecayProductFinalState",&Lepton_gen_isDirectPromptTauDecayProductFinalState);
    // new variables
    newtree->Branch("Lepton_cut",&Lepton_cut);
    newtree->Branch("Lepton_BDT",&Lepton_BDT);
    newtree->Branch("Lepton_corrpt",&Lepton_corrpt);
    newtree->Branch("Lepton_FR",&Lepton_FR);
    newtree->Branch("Lepton_CF",&Lepton_CF);
    newtree->Branch("Lepton_passConversion",&Lepton_passConversion);
    newtree->Branch("Lepton_passMuTightCharge",&Lepton_passMuTightCharge);
    newtree->Branch("Lepton_passEleTightCharge",&Lepton_passEleTightCharge);
    newtree->Branch("Lepton_passMissHit",&Lepton_passMissHit);
    newtree->Branch("Lepton_isMatchRightCharge",&Lepton_isMatchRightCharge);
    //Tau
    newtree->Branch("Tau_pt",&Tau_pt);
    newtree->Branch("Tau_eta",&Tau_eta);
    newtree->Branch("Tau_phi",&Tau_phi);
    newtree->Branch("Tau_energy",&Tau_energy);
    newtree->Branch("Tau_charge",&Tau_charge);
    newtree->Branch("Tau_packedLeadTauCand_dz",&Tau_packedLeadTauCand_dz);
    newtree->Branch("Tau_packedLeadTauCand_dxy",&Tau_packedLeadTauCand_dxy);
    newtree->Branch("Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    newtree->Branch("Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    newtree->Branch("Tau_decayModeFinding",&Tau_decayModeFinding);
    newtree->Branch("Tau_cut",&Tau_cut);
    // Jets
    newtree->Branch("Jet_pt",&Jet_pt);
    newtree->Branch("Jet_eta",&Jet_eta);
    newtree->Branch("Jet_phi",&Jet_phi);
    newtree->Branch("Jet_energy",&Jet_energy);
    newtree->Branch("Jet_genMother_pt",&Jet_genMother_pt);
    newtree->Branch("Jet_genMother_eta",&Jet_genMother_eta);
    newtree->Branch("Jet_genMother_phi",&Jet_genMother_phi);
    newtree->Branch("Jet_genMother_en",&Jet_genMother_en);
    newtree->Branch("Jet_genMother_pdgId",&Jet_genMother_pdgId);
    newtree->Branch("Jet_genGrandMother_pt",&Jet_genGrandMother_pt);
    newtree->Branch("Jet_genGrandMother_eta",&Jet_genGrandMother_eta);
    newtree->Branch("Jet_genGrandMother_phi",&Jet_genGrandMother_phi);
    newtree->Branch("Jet_genGrandMother_en",&Jet_genGrandMother_en);
    newtree->Branch("Jet_genGrandMother_pdgId",&Jet_genGrandMother_pdgId);
    newtree->Branch("Jet_Uncorr_pt",&Jet_Uncorr_pt);
    newtree->Branch("Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags",&Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
    newtree->Branch("Jet_pfCombinedMVAV2BJetTags",&Jet_pfCombinedMVAV2BJetTags);
    newtree->Branch("Jet_qg",&Jet_qg);
    newtree->Branch("Jet_axis2",&Jet_axis2);
    newtree->Branch("Jet_ptD",&Jet_ptD);
    newtree->Branch("Jet_mult",&Jet_mult);
    newtree->Branch("Jet_partonFlavour",&Jet_partonFlavour);
    newtree->Branch("Jet_hadronFlavour",&Jet_hadronFlavour);
    newtree->Branch("Jet_genpt",&Jet_genpt);
    newtree->Branch("Jet_geneta",&Jet_geneta);
    newtree->Branch("Jet_genphi",&Jet_genphi);
    newtree->Branch("Jet_genenergy",&Jet_genenergy);
    // new variables
    newtree->Branch("Jet_cut",&Jet_cut);
    newtree->Branch("Jet_bcut",&Jet_bcut);
    //BoostJet
    newtree->Branch("BoostedJet_pt",&BoostedJet_pt);
    newtree->Branch("BoostedJet_eta",&BoostedJet_eta);
    newtree->Branch("BoostedJet_phi",&BoostedJet_phi);
    newtree->Branch("BoostedJet_energy",&BoostedJet_energy);
    newtree->Branch("BoostedJet_Uncorr_pt",&BoostedJet_Uncorr_pt);
    newtree->Branch("BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags",&BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags);
    newtree->Branch("BoostedJet_pfCombinedMVAV2BJetTags",&BoostedJet_pfCombinedMVAV2BJetTags);
    newtree->Branch("BoostedJet_JesSF",&BoostedJet_JesSF);
    newtree->Branch("BoostedJet_JesSFup",&BoostedJet_JesSFup);
    newtree->Branch("BoostedJet_JesSFdown",&BoostedJet_JesSFdown);
    newtree->Branch("BoostedJet_JerSF",&BoostedJet_JerSF);
    newtree->Branch("BoostedJet_JerSFup",&BoostedJet_JerSFup);
    newtree->Branch("BoostedJet_JerSFdown",&BoostedJet_JerSFdown);
    newtree->Branch("BoostedJet_tau1",&BoostedJet_tau1);
    newtree->Branch("BoostedJet_tau2",&BoostedJet_tau2);
    newtree->Branch("BoostedJet_tau3",&BoostedJet_tau3);
    newtree->Branch("BoostedJet_softdrop_mass",&BoostedJet_softdrop_mass);
    newtree->Branch("BoostedJet_pruned_mass",&BoostedJet_pruned_mass);
    newtree->Branch("BoostedJet_tau21",&BoostedJet_tau21);
    newtree->Branch("BoostedJet_tau32",&BoostedJet_tau32);
    newtree->Branch("BoostedJet_wCut",&BoostedJet_wCut);
    newtree->Branch("BoostedJet_topCut",&BoostedJet_topCut);
    //Gen
    newtree->Branch("Gen_pdg_id",&Gen_pdg_id);
    newtree->Branch("Gen_pt",&Gen_pt);
    newtree->Branch("Gen_eta",&Gen_eta);
    newtree->Branch("Gen_phi",&Gen_phi);
    newtree->Branch("Gen_energy",&Gen_energy);
    newtree->Branch("Gen_motherpdg_id",&Gen_motherpdg_id);
    newtree->Branch("Gen_BmotherIndex",&Gen_BmotherIndex);
    newtree->Branch("Gen_numMother",&Gen_numMother);
    newtree->Branch("Gen_status",&Gen_status);
};


void wClearInitialization(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //Initialize Number
    EVENT_event= -999;
    EVENT_genWeight= -999;
    HiggsDecay= -999;
    PUWeight= -999;
    Gen_type1PF_Met= -999;
    Gen_type1PF_Metpx= -999;
    Gen_type1PF_Metpy= -999;
    Gen_type1PF_Metpz= -999;
    Gen_type1PF_Meteta= -999;
    Gen_type1PF_Metphi= -999;
    Gen_type1PF_Meten= -999;
    Met_type1PF_pt= -999;
    Met_type1PF_px= -999;
    Met_type1PF_py= -999;
    Met_type1PF_pz= -999;
    Met_type1PF_phi= -999;
    Met_type1PF_sumEt= -999;
    Met_type1PF_shiftedPtUp= -999;
    Met_type1PF_shiftedPtDown= -999;
    Muon_numLoose= -999;
    Muon_numFake= -999;
    Muon_numTight= -999;
    patElectron_numLoose= -999;
    patElectron_numFake= -999;
    patElectron_numTight= -999;
    Tau_numLoose= -999;
    Tau_numMedium= -999;
    Jet_numLoose= -999;
    Jet_numbLoose= -999;
    Jet_numbMedium= -999;
    Jet_numbTight= -999;
    BWeight= -999;
    BWeightLFup= -999;
    BWeightLFdown= -999;
    BWeightHFup= -999;
    BWeightHFdown= -999;
    BWeightJESup= -999;
    BWeightJESdown= -999;
    BWeightLFStats1up= -999;
    BWeightLFStats1down= -999;
    BWeightLFStats2up= -999;
    BWeightLFStats2down= -999;
    BWeightHFStats1up= -999;
    BWeightHFStats1down= -999;
    BWeightHFStats2up= -999;
    BWeightHFStats2down= -999;
    BWeightCErr1up= -999;
    BWeightCErr1down= -999;
    BWeightCErr2up= -999;
    BWeightCErr2down= -999;
    Top_numSoft= -999;
    Top_numLoose= -999;
    Top_numMedium= -999;
    Top_numTight= -999;
    W_numLoose= -999;
    W_numTight= -999;
    TTHLep_2Mu= -999;
    TTHLep_2Ele= -999;
    TTHLep_MuEle= -999;
    TTHLep_3L4L= -999;
    metLD= -999;
    mhtT_met= -999;
    mht_met= -999;
    mhtT= -999;
    mht= -999;
    // Lepton
    leptons->clear();
    Lepton_pt->clear();
    Lepton_eta->clear();
    Lepton_phi->clear();
    Lepton_energy->clear();
    Lepton_dxy_pv->clear();
    Lepton_dz_pv->clear();
    Lepton_IP3Dsig->clear();
    Lepton_loose->clear();
    Lepton_miniIsoRel->clear();
    Lepton_charge->clear();
    Lepton_pdgId->clear();
    Lepton_isGlobal->clear();
    Lepton_chi2->clear();
    Lepton_chi2LocalPosition->clear();
    Lepton_trkKink->clear();
    Lepton_validFraction->clear();
    Lepton_segmentCompatibility->clear();
    Lepton_jetptratio->clear();
    Lepton_jetcsv->clear();
    Lepton_jetpt->clear();
    Lepton_lepjetchtrks->clear();
    Lepton_miniIsoCh->clear();
    Lepton_miniIsoPUsub->clear();
    Lepton_ptrel->clear();
    Lepton_pTErrOVpT_it->clear();
    Lepton_px->clear();
    Lepton_py->clear();
    Lepton_pz->clear();
    Lepton_jetdr->clear();
    Lepton_SCeta->clear();
    Lepton_mvaValue_HZZ->clear();
    Lepton_expectedMissingInnerHits->clear();
    Lepton_full5x5_sigmaIetaIeta->clear();
    Lepton_hOverE->clear();
    Lepton_dEtaIn->clear();
    Lepton_dPhiIn->clear();
    Lepton_ooEmooP->clear();
    // gen
    Lepton_gen_pt->clear();
    Lepton_gen_eta->clear();
    Lepton_gen_phi->clear();
    Lepton_gen_en->clear();
    Lepton_gen_pdgId->clear();
    Lepton_genMother_pt->clear();
    Lepton_genMother_eta->clear();
    Lepton_genMother_phi->clear();
    Lepton_genMother_en->clear();
    Lepton_genMother_pdgId->clear();
    Lepton_genGrandMother_pt->clear();
    Lepton_genGrandMother_eta->clear();
    Lepton_genGrandMother_phi->clear();
    Lepton_genGrandMother_en->clear();
    Lepton_genGrandMother_pdgId->clear();
    Lepton_gen_isPromptFinalState->clear();
    Lepton_gen_isDirectPromptTauDecayProductFinalState->clear();
    // new variables
    Lepton_cut->clear();
    Lepton_BDT->clear();
    Lepton_corrpt->clear();
    Lepton_FR->clear();
    Lepton_CF->clear();
    Lepton_passConversion->clear();
    Lepton_passMuTightCharge->clear();
    Lepton_passEleTightCharge->clear();
    Lepton_passMissHit->clear();
    Lepton_isMatchRightCharge->clear();
    //Tau
    taus->clear();
    Tau_pt->clear();
    Tau_eta->clear();
    Tau_phi->clear();
    Tau_energy->clear();
    Tau_charge->clear();
    Tau_packedLeadTauCand_dz->clear();
    Tau_packedLeadTauCand_dxy->clear();
    Tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->clear();
    Tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->clear();
    Tau_decayModeFinding->clear();
    Tau_cut->clear();
    //Jet
    jets->clear();
    Jet_pt->clear();
    Jet_eta->clear();
    Jet_phi->clear();
    Jet_energy->clear();
    Jet_genMother_pt->clear();
    Jet_genMother_eta->clear();
    Jet_genMother_phi->clear();
    Jet_genMother_en->clear();
    Jet_genMother_pdgId->clear();
    Jet_genGrandMother_pt->clear();
    Jet_genGrandMother_eta->clear();
    Jet_genGrandMother_phi->clear();
    Jet_genGrandMother_en->clear();
    Jet_genGrandMother_pdgId->clear();
    Jet_Uncorr_pt->clear();
    Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags->clear();
    Jet_pfCombinedMVAV2BJetTags->clear();
    Jet_qg->clear();
    Jet_axis2->clear();
    Jet_ptD->clear();
    Jet_mult->clear();
    Jet_partonFlavour->clear();
    Jet_hadronFlavour->clear();
    Jet_genpt->clear();
    Jet_geneta->clear();
    Jet_genphi->clear();
    Jet_genenergy->clear();
    Jet_cut->clear();
    Jet_bcut->clear();
    // BoostJet
    boostjets->clear();
    BoostedJet_pt->clear();
    BoostedJet_eta->clear();
    BoostedJet_phi->clear();
    BoostedJet_energy->clear();
    BoostedJet_Uncorr_pt->clear();
    BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->clear();
    BoostedJet_pfCombinedMVAV2BJetTags->clear();
    BoostedJet_JesSF->clear();
    BoostedJet_JesSFup->clear();
    BoostedJet_JesSFdown->clear();
    BoostedJet_JerSF->clear();
    BoostedJet_JerSFup->clear();
    BoostedJet_JerSFdown->clear();
    BoostedJet_tau1->clear();
    BoostedJet_tau2->clear();
    BoostedJet_tau3->clear();
    BoostedJet_softdrop_mass->clear();
    BoostedJet_pruned_mass->clear();
    BoostedJet_tau21->clear();
    BoostedJet_tau32->clear();
    BoostedJet_wCut->clear();
    BoostedJet_topCut->clear();
    // Gen
    Gen_pdg_id->clear();
    Gen_pt->clear();
    Gen_eta->clear();
    Gen_phi->clear();
    Gen_energy->clear();
    Gen_motherpdg_id->clear();
    Gen_BmotherIndex->clear();
    Gen_numMother->clear();
    Gen_status->clear();
};


void rGetEntry(Long64_t tentry, string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //GetEntry
    b_rEVENT_event->GetEntry(tentry);
    b_rEVENT_genWeight->GetEntry(tentry);
    b_rHiggsDecay->GetEntry(tentry);
    b_rPUWeight->GetEntry(tentry);
    b_rGen_type1PF_Met->GetEntry(tentry);
    b_rGen_type1PF_Metpx->GetEntry(tentry);
    b_rGen_type1PF_Metpy->GetEntry(tentry);
    b_rGen_type1PF_Metpz->GetEntry(tentry);
    b_rGen_type1PF_Meteta->GetEntry(tentry);
    b_rGen_type1PF_Metphi->GetEntry(tentry);
    b_rGen_type1PF_Meten->GetEntry(tentry);
    b_rMet_type1PF_pt->GetEntry(tentry);
    b_rMet_type1PF_px->GetEntry(tentry);
    b_rMet_type1PF_py->GetEntry(tentry);
    b_rMet_type1PF_pz->GetEntry(tentry);
    b_rMet_type1PF_phi->GetEntry(tentry);
    b_rMet_type1PF_sumEt->GetEntry(tentry);
    b_rMet_type1PF_shiftedPtUp->GetEntry(tentry);
    b_rMet_type1PF_shiftedPtDown->GetEntry(tentry);
    b_rHLT_DiMu9_Ele9_CaloIdL_TrackIdL->GetEntry(tentry);
    b_rHLT_Mu8_DiEle12_CaloIdL_TrackIdL->GetEntry(tentry);
    b_rHLT_TripleMu_12_10_5->GetEntry(tentry);
    b_rHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL->GetEntry(tentry);
    b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL->GetEntry(tentry);
    b_rHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ->GetEntry(tentry);
    b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL->GetEntry(tentry);
    b_rHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ->GetEntry(tentry);
    b_rHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->GetEntry(tentry);
    b_rHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ->GetEntry(tentry);
    b_rHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ->GetEntry(tentry);
    b_rHLT_IsoMu22->GetEntry(tentry);
    b_rHLT_IsoTkMu22->GetEntry(tentry);
    b_rHLT_IsoMu22_eta2p1->GetEntry(tentry);
    b_rHLT_IsoTkMu22_eta2p1->GetEntry(tentry);
    b_rHLT_IsoMu24->GetEntry(tentry);
    b_rHLT_IsoTkMu24->GetEntry(tentry);
    b_rHLT_Ele27_WPTight_Gsf->GetEntry(tentry);
    b_rHLT_Ele25_eta2p1_WPTight_Gsf->GetEntry(tentry);
    b_rHLT_Ele27_eta2p1_WPLoose_Gsf->GetEntry(tentry);
    // Muon
    b_rMuon_pt->GetEntry(tentry);
    b_rMuon_eta->GetEntry(tentry);
    b_rMuon_phi->GetEntry(tentry);
    b_rMuon_energy->GetEntry(tentry);
    b_rMuon_dxy_pv->GetEntry(tentry);
    b_rMuon_dz_pv->GetEntry(tentry);
    b_rMuon_IP3Dsig_it->GetEntry(tentry);
    b_rMuon_loose->GetEntry(tentry);
    b_rMuon_miniIsoRel->GetEntry(tentry);
    b_rMuon_charge->GetEntry(tentry);
    b_rMuon_pdgId->GetEntry(tentry);
    b_rMuon_isGlobal->GetEntry(tentry);
    b_rMuon_chi2->GetEntry(tentry);
    b_rMuon_chi2LocalPosition->GetEntry(tentry);
    b_rMuon_trkKink->GetEntry(tentry);
    b_rMuon_validFraction->GetEntry(tentry);
    b_rMuon_segmentCompatibility->GetEntry(tentry);
    b_rMuon_jetptratio->GetEntry(tentry);
    b_rMuon_jetcsv->GetEntry(tentry);
    b_rMuon_jetpt->GetEntry(tentry);
    b_rMuon_lepjetchtrks->GetEntry(tentry);
    b_rMuon_miniIsoCh->GetEntry(tentry);
    b_rMuon_miniIsoPUsub->GetEntry(tentry);
    b_rMuon_ptrel->GetEntry(tentry);
    b_rMuon_pTErrOVpT_it->GetEntry(tentry);
    b_rMuon_px->GetEntry(tentry);
    b_rMuon_py->GetEntry(tentry);
    b_rMuon_pz->GetEntry(tentry);
    b_rMuon_jetdr->GetEntry(tentry);
    b_rMuon_gen_pt->GetEntry(tentry);
    b_rMuon_gen_eta->GetEntry(tentry);
    b_rMuon_gen_phi->GetEntry(tentry);
    b_rMuon_gen_en->GetEntry(tentry);
    b_rMuon_gen_pdgId->GetEntry(tentry);
    b_rMuon_genMother_pt->GetEntry(tentry);
    b_rMuon_genMother_eta->GetEntry(tentry);
    b_rMuon_genMother_phi->GetEntry(tentry);
    b_rMuon_genMother_en->GetEntry(tentry);
    b_rMuon_genMother_pdgId->GetEntry(tentry);
    b_rMuon_genGrandMother_pt->GetEntry(tentry);
    b_rMuon_genGrandMother_eta->GetEntry(tentry);
    b_rMuon_genGrandMother_phi->GetEntry(tentry);
    b_rMuon_genGrandMother_en->GetEntry(tentry);
    b_rMuon_genGrandMother_pdgId->GetEntry(tentry);
    b_rMuon_gen_isPromptFinalState->GetEntry(tentry);
    b_rMuon_gen_isDirectPromptTauDecayProductFinalState->GetEntry(tentry);
    //Electron
    b_rpatElectron_pt->GetEntry(tentry);
    b_rpatElectron_eta->GetEntry(tentry);
    b_rpatElectron_phi->GetEntry(tentry);
    b_rpatElectron_energy->GetEntry(tentry);
    b_rpatElectron_IP3Dsig->GetEntry(tentry);
    b_rpatElectron_miniIsoRel->GetEntry(tentry);
    b_rpatElectron_charge->GetEntry(tentry);
    b_rpatElectron_pdgId->GetEntry(tentry);
    b_rpatElectron_gsfTrack_dxy_pv->GetEntry(tentry);
    b_rpatElectron_gsfTrack_dz_pv->GetEntry(tentry);
    b_rpatElectron_jetptratio->GetEntry(tentry);
    b_rpatElectron_jetcsv->GetEntry(tentry);
    b_rpatElectron_jetpt->GetEntry(tentry);
    b_rpatElectron_lepjetchtrks->GetEntry(tentry);
    b_rpatElectron_miniIsoCh->GetEntry(tentry);
    b_rpatElectron_miniIsoPUsub->GetEntry(tentry);
    b_rpatElectron_ptrel->GetEntry(tentry);
    b_rpatElectron_px->GetEntry(tentry);
    b_rpatElectron_py->GetEntry(tentry);
    b_rpatElectron_pz->GetEntry(tentry);
    b_rpatElectron_jetdr->GetEntry(tentry);
    b_rpatElectron_gen_pt->GetEntry(tentry);
    b_rpatElectron_gen_eta->GetEntry(tentry);
    b_rpatElectron_gen_phi->GetEntry(tentry);
    b_rpatElectron_gen_en->GetEntry(tentry);
    b_rpatElectron_gen_pdgId->GetEntry(tentry);
    b_rpatElectron_genMother_pt->GetEntry(tentry);
    b_rpatElectron_genMother_eta->GetEntry(tentry);
    b_rpatElectron_genMother_phi->GetEntry(tentry);
    b_rpatElectron_genMother_en->GetEntry(tentry);
    b_rpatElectron_genMother_pdgId->GetEntry(tentry);
    b_rpatElectron_genGrandMother_pt->GetEntry(tentry);
    b_rpatElectron_genGrandMother_eta->GetEntry(tentry);
    b_rpatElectron_genGrandMother_phi->GetEntry(tentry);
    b_rpatElectron_genGrandMother_en->GetEntry(tentry);
    b_rpatElectron_genGrandMother_pdgId->GetEntry(tentry);
    b_rpatElectron_gen_isPromptFinalState->GetEntry(tentry);
    b_rpatElectron_gen_isDirectPromptTauDecayProductFinalState->GetEntry(tentry);
    b_rpatElectron_SCeta->GetEntry(tentry);
    b_rpatElectron_expectedMissingInnerHits->GetEntry(tentry);
    b_rpatElectron_full5x5_sigmaIetaIeta->GetEntry(tentry);
    b_rpatElectron_hOverE->GetEntry(tentry);
    b_rpatElectron_dEtaIn->GetEntry(tentry);
    b_rpatElectron_dPhiIn->GetEntry(tentry);
    b_rpatElectron_ooEmooP->GetEntry(tentry);
    b_rpatElectron_mvaValue_HZZ->GetEntry(tentry);
    b_rpatElectron_isGsfCtfScPixChargeConsistent->GetEntry(tentry);
    b_rpatElectron_isGsfScPixChargeConsistent->GetEntry(tentry);
    b_rpatElectron_passConversionVeto->GetEntry(tentry);
    //Tau
    b_rTau_pt->GetEntry(tentry);
    b_rTau_eta->GetEntry(tentry);
    b_rTau_phi->GetEntry(tentry);
    b_rTau_energy->GetEntry(tentry);
    b_rTau_charge->GetEntry(tentry);
    b_rTau_packedLeadTauCand_dz->GetEntry(tentry);
    b_rTau_packedLeadTauCand_dxy->GetEntry(tentry);
    b_rTau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->GetEntry(tentry);
    b_rTau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->GetEntry(tentry);
    b_rTau_decayModeFinding->GetEntry(tentry);
    //Jet
    b_rJet_pt->GetEntry(tentry);
    b_rJet_eta->GetEntry(tentry);
    b_rJet_phi->GetEntry(tentry);
    b_rJet_energy->GetEntry(tentry);
    b_rJet_genMother_pt->GetEntry(tentry);
    b_rJet_genMother_eta->GetEntry(tentry);
    b_rJet_genMother_phi->GetEntry(tentry);
    b_rJet_genMother_en->GetEntry(tentry);
    b_rJet_genMother_pdgId->GetEntry(tentry);
    b_rJet_genGrandMother_pt->GetEntry(tentry);
    b_rJet_genGrandMother_eta->GetEntry(tentry);
    b_rJet_genGrandMother_phi->GetEntry(tentry);
    b_rJet_genGrandMother_en->GetEntry(tentry);
    b_rJet_genGrandMother_pdgId->GetEntry(tentry);
    b_rJet_Uncorr_pt->GetEntry(tentry);
    b_rJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->GetEntry(tentry);
    b_rJet_pfCombinedMVAV2BJetTags->GetEntry(tentry);
    b_rJet_qg->GetEntry(tentry);
    b_rJet_axis2->GetEntry(tentry);
    b_rJet_ptD->GetEntry(tentry);
    b_rJet_mult->GetEntry(tentry);
    b_rJet_partonFlavour->GetEntry(tentry);
    b_rJet_hadronFlavour->GetEntry(tentry);
    b_rJet_genpt->GetEntry(tentry);
    b_rJet_geneta->GetEntry(tentry);
    b_rJet_genphi->GetEntry(tentry);
    b_rJet_genenergy->GetEntry(tentry);
    b_rJet_JesSF->GetEntry(tentry);
    b_rJet_JesSFup->GetEntry(tentry);
    b_rJet_JesSFdown->GetEntry(tentry);
    b_rJet_JerSF->GetEntry(tentry);
    b_rJet_JerSFup->GetEntry(tentry);
    b_rJet_JerSFdown->GetEntry(tentry);
    b_rJet_neutralHadEnergyFraction->GetEntry(tentry);
    b_rJet_neutralEmEnergyFraction->GetEntry(tentry);
    b_rJet_chargedMultiplicity->GetEntry(tentry);
    b_rJet_numberOfConstituents->GetEntry(tentry);
    b_rJet_chargedHadronEnergyFraction->GetEntry(tentry);
    b_rJet_chargedEmEnergyFraction->GetEntry(tentry);
    b_rJet_btag_sf->GetEntry(tentry);
    b_rJet_btag_jesup->GetEntry(tentry);
    b_rJet_btag_jesdown->GetEntry(tentry);
    b_rJet_btag_hfup->GetEntry(tentry);
    b_rJet_btag_hfdown->GetEntry(tentry);
    b_rJet_btag_hfstat1up->GetEntry(tentry);
    b_rJet_btag_hfstat1down->GetEntry(tentry);
    b_rJet_btag_hfstat2up->GetEntry(tentry);
    b_rJet_btag_hfstat2down->GetEntry(tentry);
    b_rJet_btag_lfup->GetEntry(tentry);
    b_rJet_btag_lfdown->GetEntry(tentry);
    b_rJet_btag_lfstat1up->GetEntry(tentry);
    b_rJet_btag_lfstat1down->GetEntry(tentry);
    b_rJet_btag_lfstat2up->GetEntry(tentry);
    b_rJet_btag_lfstat2down->GetEntry(tentry);
    b_rJet_btag_cerr1up->GetEntry(tentry);
    b_rJet_btag_cerr1down->GetEntry(tentry);
    b_rJet_btag_cerr2up->GetEntry(tentry);
    b_rJet_btag_cerr2down->GetEntry(tentry);
    //BoostJet
    b_rBoostedJet_pt->GetEntry(tentry);
    b_rBoostedJet_eta->GetEntry(tentry);
    b_rBoostedJet_phi->GetEntry(tentry);
    b_rBoostedJet_energy->GetEntry(tentry);
    b_rBoostedJet_Uncorr_pt->GetEntry(tentry);
    b_rBoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags->GetEntry(tentry);
    b_rBoostedJet_pfCombinedMVAV2BJetTags->GetEntry(tentry);
    b_rBoostedJet_JesSF->GetEntry(tentry);
    b_rBoostedJet_JesSFup->GetEntry(tentry);
    b_rBoostedJet_JesSFdown->GetEntry(tentry);
    b_rBoostedJet_JerSF->GetEntry(tentry);
    b_rBoostedJet_JerSFup->GetEntry(tentry);
    b_rBoostedJet_JerSFdown->GetEntry(tentry);
    b_rBoostedJet_tau1->GetEntry(tentry);
    b_rBoostedJet_tau2->GetEntry(tentry);
    b_rBoostedJet_tau3->GetEntry(tentry);
    b_rBoostedJet_softdrop_mass->GetEntry(tentry);
    b_rBoostedJet_pruned_mass->GetEntry(tentry);
    //Gen
    b_rGen_pdg_id->GetEntry(tentry);
    b_rGen_pt->GetEntry(tentry);
    b_rGen_eta->GetEntry(tentry);
    b_rGen_phi->GetEntry(tentry);
    b_rGen_energy->GetEntry(tentry);
    b_rGen_motherpdg_id->GetEntry(tentry);
    b_rGen_BmotherIndex->GetEntry(tentry);
    b_rGen_numMother->GetEntry(tentry);
    b_rGen_status->GetEntry(tentry);
};
