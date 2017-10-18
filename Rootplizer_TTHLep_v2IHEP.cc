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
        Muon.CF = 1.;
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
       
        // calculate new variables 
        // pass ismedium boolean, always true for electron
        patElectron.set_Wp_tthlep(true, ele_numLoose, ele_numFake, ele_numTight ); 
    }
    patElectron_numLoose = ele_numLoose;
    patElectron_numFake =  ele_numFake;
    patElectron_numTight = ele_numTight;
}


void Lep_sel(){
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
    Met_type1PF_pt = rMet_type1PF_pt;
    Met_type1PF_px = rMet_type1PF_px;
    Met_type1PF_py = rMet_type1PF_py;
    Met_type1PF_pz = rMet_type1PF_pz;
    Met_type1PF_phi = rMet_type1PF_phi;
    Met_type1PF_sumEt = rMet_type1PF_sumEt;
    Met_type1PF_shiftedPtUp = rMet_type1PF_shiftedPtUp;
    Met_type1PF_shiftedPtDown = rMet_type1PF_shiftedPtDown;
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
    readingtree->SetBranchAddress("Met_type1PF_pt",&rMet_type1PF_pt,&b_rMet_type1PF_pt);
    readingtree->SetBranchAddress("Met_type1PF_px",&rMet_type1PF_px,&b_rMet_type1PF_px);
    readingtree->SetBranchAddress("Met_type1PF_py",&rMet_type1PF_py,&b_rMet_type1PF_py);
    readingtree->SetBranchAddress("Met_type1PF_pz",&rMet_type1PF_pz,&b_rMet_type1PF_pz);
    readingtree->SetBranchAddress("Met_type1PF_phi",&rMet_type1PF_phi,&b_rMet_type1PF_phi);
    readingtree->SetBranchAddress("Met_type1PF_sumEt",&rMet_type1PF_sumEt,&b_rMet_type1PF_sumEt);
    readingtree->SetBranchAddress("Met_type1PF_shiftedPtUp",&rMet_type1PF_shiftedPtUp,&b_rMet_type1PF_shiftedPtUp);
    readingtree->SetBranchAddress("Met_type1PF_shiftedPtDown",&rMet_type1PF_shiftedPtDown,&b_rMet_type1PF_shiftedPtDown);
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
};


void wSetBranchAddress(TTree* newtree, string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //Write setbranchaddress
    newtree->Branch("EVENT_event",&EVENT_event);
    newtree->Branch("EVENT_genWeight",&EVENT_genWeight);
    newtree->Branch("HiggsDecay",&HiggsDecay);
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
};


void wClearInitialization(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //Initialize Number
    EVENT_event= -999;
    EVENT_genWeight= -999;
    HiggsDecay= -999;
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
};


void rGetEntry(Long64_t tentry, string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    //GetEntry
    b_rEVENT_event->GetEntry(tentry);
    b_rEVENT_genWeight->GetEntry(tentry);
    b_rHiggsDecay->GetEntry(tentry);
    b_rMet_type1PF_pt->GetEntry(tentry);
    b_rMet_type1PF_px->GetEntry(tentry);
    b_rMet_type1PF_py->GetEntry(tentry);
    b_rMet_type1PF_pz->GetEntry(tentry);
    b_rMet_type1PF_phi->GetEntry(tentry);
    b_rMet_type1PF_sumEt->GetEntry(tentry);
    b_rMet_type1PF_shiftedPtUp->GetEntry(tentry);
    b_rMet_type1PF_shiftedPtDown->GetEntry(tentry);
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
};
