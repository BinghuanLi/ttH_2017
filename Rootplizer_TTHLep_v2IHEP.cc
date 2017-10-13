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
    //Fill new branches
    int nen = nentries; if(nentries==-1) nen = readingtree->GetEntries();
    for(Int_t en=0; en<nen; en++){
        //Initialization
        wClearInitialization(sample);
        Long64_t tentry = readingtree->LoadTree(en); 
        rGetEntry(tentry, sample);
        //Muon
        Muon_sel(sample);
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
//Object Selection
////
//Muon
void Muon_sel(string sample){
    if(sample == "data"){
        cout<<" Hi, I'm data!"<<endl;
    }
    for(uint mu_en = 0; mu_en<rMuon_pt->size(); mu_en++){
        Lepton Muon;
        // initialize needed variables for mu_isLoose_tthlep()
        Muon.pt= rMuon_pt->at(mu_en);
        Muon.eta= rMuon_eta->at(mu_en);
        Muon.phi= rMuon_phi->at(mu_en);
        Muon.energy= rMuon_energy->at(mu_en);
        Muon.dxy_pv= rMuon_dxy_pv->at(mu_en);
        Muon.dz_pv= rMuon_dz_pv->at(mu_en);
        Muon.IP3Dsig_it= rMuon_IP3Dsig_it->at(mu_en);
        Muon.loose= rMuon_loose->at(mu_en);
        Muon.miniIsoRel= rMuon_miniIsoRel->at(mu_en);
        // check whether Muon pass tth loose lepton selection
        if(!(Muon.mu_isLoose_tthlep()))continue;
        Muon.charge= rMuon_charge->at(mu_en);
        Muon.pdgId= rMuon_pdgId->at(mu_en);
        Muon.set_Wp_tthlep();
        leptons->push_back(Muon);    
    }
}


void Lep_sel(){
    for(uint lep_en=0; lep_en < leptons->size(); lep_en++){
        Lepton_pt->push_back(leptons->at(lep_en).pt);
        Lepton_eta->push_back(leptons->at(lep_en).eta);
        Lepton_phi->push_back(leptons->at(lep_en).phi);
        Lepton_energy->push_back(leptons->at(lep_en).energy);
        Lepton_dxy_pv->push_back(leptons->at(lep_en).dxy_pv);
        Lepton_dz_pv->push_back(leptons->at(lep_en).dz_pv);
        Lepton_IP3Dsig_it->push_back(leptons->at(lep_en).IP3Dsig_it);
        Lepton_loose->push_back(leptons->at(lep_en).loose);
        Lepton_miniIsoRel->push_back(leptons->at(lep_en).loose);
        Lepton_charge->push_back(leptons->at(lep_en).charge);
        Lepton_pdgId->push_back(leptons->at(lep_en).pdgId);
        Lepton_cut->push_back(leptons->at(lep_en).cut);
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
    //Lepton
    newtree->Branch("Lepton_pt",&Lepton_pt);
    newtree->Branch("Lepton_eta",&Lepton_eta);
    newtree->Branch("Lepton_phi",&Lepton_phi);
    newtree->Branch("Lepton_energy",&Lepton_energy);
    newtree->Branch("Lepton_dxy_pv",&Lepton_dxy_pv);
    newtree->Branch("Lepton_dz_pv",&Lepton_dz_pv);
    newtree->Branch("Lepton_IP3Dsig_it",&Lepton_IP3Dsig_it);
    newtree->Branch("Lepton_loose",&Lepton_loose);
    newtree->Branch("Lepton_miniIsoRel",&Lepton_miniIsoRel);
    newtree->Branch("Lepton_charge",&Lepton_charge);
    newtree->Branch("Lepton_pdgId",&Lepton_pdgId);
    // new variables
    newtree->Branch("Lepton_cut",&Lepton_cut);
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
    // Lepton
    leptons->clear();
    Lepton_pt->clear();
    Lepton_eta->clear();
    Lepton_phi->clear();
    Lepton_energy->clear();
    Lepton_dxy_pv->clear();
    Lepton_dz_pv->clear();
    Lepton_IP3Dsig_it->clear();
    Lepton_loose->clear();
    Lepton_miniIsoRel->clear();
    Lepton_charge->clear();
    Lepton_pdgId->clear();
    // new variables
    Lepton_cut->clear();
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
    //Muon
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
};
