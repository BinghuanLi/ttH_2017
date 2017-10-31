import os
import re
import sys
import commands
import subprocess

#########################
############### Information provided by user
#################################
# Path you run this script
workpath = "/publicfs/cms/user/libh/Test/Rootplizer/analyzer"
#Specify needed variables
#  Type   Reading  Saving  Example
#  CaseA   Yes      Yes     Jet_pt: you want to save it and need also to read it (using SetBranchStatus+CopyEntries would require more processing time) 
#  CaseB   Yes      No      Jet_muonEnergyFraction: you want to read it to define jet ID, but may not want to save it  
#  CaseC   No       Yes     Jet_LooseID: you can not read it from ntuplas if there is not, but you want to save in rootplas 

# Variable Case
#Case = "CaseA"
#Case = "CaseB"
Case = "CaseC"
# Variable Definition

VariableType  = "double"
VariableNames = [
# CaseA Variables
#"Muon_pt","Muon_eta","Muon_phi","Muon_energy","Muon_px",
#"Muon_py","Muon_pz","Muon_p","Muon_dB","Muon_pt_it",
#"Muon_ptErr_it","Muon_pt_bt","Muon_pTErrOVpT_it","Muon_ptErr_bt","Muon_pTErrOVpT_bt",
#"Muon_pt_tunePbt","Muon_charge","Muon_isoR04Charged","Muon_isoR04NeutralHadron","Muon_isoR04Photon",
#"Muon_isoR04PU","Muon_relIsoDeltaBetaR04","Muon_isoR04CharParPt","Muon_isoR03Charged","Muon_isoR03NeutralHadron",
#"Muon_isoR03Photon","Muon_isoR03PU","Muon_relIsoDeltaBetaR03","Muon_isoR03CharParPt","Muon_trackIso",
#"Muon_TrackerIso","Muon_ecalIso","Muon_hcalIso","Muon_caloIso","Muon_isoSum",
#"Muon_pfEcalEnergy","Muon_chi2","Muon_chi2LocalPosition","Muon_matchedStat","Muon_validHits",
#"Muon_validHitsInner","Muon_TLayers","Muon_ndof","Muon_validFraction","Muon_pixelLayersWithMeasurement",
#"Muon_qualityhighPurity","Muon_trkKink","Muon_segmentCompatibility","Muon_dz_pv","Muon_dxy_pv",
#"Muon_dz_bt","Muon_dxy_bt","Muon_dz_bs","Muon_dxy_bs","Muon_dzError",
#"Muon_dxyError","Muon_vtx","Muon_vty","Muon_vtz","Muon_track_PCAx_bs",
#"Muon_track_PCAy_bs","Muon_track_PCAz_bs","Muon_track_PCAx_pv","Muon_track_PCAy_pv","Muon_track_PCAz_pv",
#"Muon_trackFitErrorMatrix_00","Muon_trackFitErrorMatrix_01","Muon_trackFitErrorMatrix_02","Muon_trackFitErrorMatrix_11","Muon_trackFitErrorMatrix_12",
#"Muon_trackFitErrorMatrix_22","Muon_miniIsoRel","Muon_miniIsoCh","Muon_miniIsoNeu","Muon_miniIsoPUsub",
#"Muon_jetdr","Muon_jetpt","Muon_jetptratio","Muon_jetcsv","Muon_ptrel",
#"Muon_IP3Dsig_it","Muon_pvass","Muon_etarel","Muon_ptOVen","Muon_mujet_pfJetProbabilityBJetTag",
#"Muon_mujet_pfCombinedMVABJetTags","Muon_mujet_qgl","Muon_mumass","Muon_mujet_mass","Muon_mujet_Wmass",
#"Muon_mujet_Topmass","Muon_mujet_WTopmass","Muon_IP3D_val","Muon_IP3D_err","Muon_IP3D_sig",
#"Muon_IP2D_val","Muon_IP2D_err","Muon_IP2D_sig","Muon_sIP3D_val","Muon_sIP3D_err",
#"Muon_sIP3D_sig","Muon_sIP2D_val","Muon_sIP2D_err","Muon_sIP2D_sig","Muon_IP1D_val",
#"Muon_IP1D_err","Muon_IP1D_sig","Muon_sIP1D_val","Muon_sIP1D_err","Muon_sIP1D_sig",
#"Muon_lepjetMaxIP3D_val","Muon_lepjetMaxIP3D_sig","Muon_lepjetMaxsIP3D_val","Muon_lepjetMaxsIP3D_sig","Muon_lepjetMaxIP2D_val",
#"Muon_lepjetMaxIP2D_sig","Muon_lepjetMaxsIP2D_val","Muon_lepjetMaxsIP2D_sig","Muon_lepjetMaxIP1D_val","Muon_lepjetMaxIP1D_sig",
#"Muon_lepjetMaxsIP1D_val","Muon_lepjetMaxsIP1D_sig","Muon_lepjetAvIP3D_val","Muon_lepjetAvIP3D_sig","Muon_lepjetAvsIP3D_val",
#"Muon_lepjetAvsIP3D_sig","Muon_lepjetAvIP2D_val","Muon_lepjetAvIP2D_sig","Muon_lepjetAvsIP2D_val","Muon_lepjetAvsIP2D_sig",
#"Muon_lepjetAvIP1D_val","Muon_lepjetAvIP1D_sig","Muon_lepjetAvsIP1D_val","Muon_lepjetAvsIP1D_sig","Muon_lepjetchtrks",
#"Muon_lepjetpvchtrks","Muon_lepjetnonpvchtrks","Muon_lepjetndaus","Muon_lepjetpvchi2","Muon_lepjetnumno2trk",
#"Muon_gen_pt","Muon_gen_eta","Muon_gen_phi","Muon_gen_en","Gen_pdg_id",
#"Gen_motherpdg_id","Gen_numMother"
#"Leph_Gen_pt","Leph_Gen_eta","Leph_Gen_phi","Leph_Gen_energy",
#"Leptop_Gen_pt","Leptop_Gen_eta","Leptop_Gen_phi","Leptop_Gen_energy",
#"MotherLeph_Gen_pt","MotherLeph_Gen_eta","MotherLeph_Gen_phi","MotherLeph_Gen_energy",
#"Leptop_Gen_Id","Leph_Gen_Id"
#"Mt_LephRecoMET",
#"Gen_pt","Gen_eta","Gen_phi","Gen_energy",
# CaseB Variables

# CaseC Variables
#"Lvh_Gen_pt","Lvh_Gen_eta","Lvh_Gen_phi","Lvh_Gen_energy","Lvh_Gen_Id",
#"Lvtop_Gen_pt","Lvtop_Gen_eta","Lvtop_Gen_phi","Lvtop_Gen_energy",
#"Lvtop_Gen_Id",
#"Lv_gmother_pdgID",
#"Lv_gmother_Index",
#"LvTag_gen_pt",
#"LvTag_Genpt",
#"LvTag_GMpdgID",
#"LvTag_GMIndex",
#"LvTag_Index",
#"Lv_Index",
#"Gen_pt","Gen_eta","Gen_phi","Gen_energy"

#"bJettopl_Gen_pt","bJettopl_Gen_eta","bJettopl_Gen_phi","bJettopl_Gen_energy","bJettopl_Gen_id",
#"Jet_origin","Jet_dr_htopb","Jet_dr_htopW","Jet_dr_ltopb","Jet_dr_higgs",
#"Jet_newpfCombinedInclusiveSecondaryVertexV2BJetTags","Jet_newpfCombinedMVAV2BJetTags",
#"Lep_BDT","FakeLep_BDT"
#"Jet_genFlavour",
#"Jet_source","Jet_ptratio_htopb","Jet_ptratio_htopW","Jet_ptratio_ltopb","Jet_ptratio_higgs","Jet_isHJetV3",

#"Lep_chRelIso","Lep_neuRelIso","Lep_jetPtRel","Lep_jetPtRatio","Lep_JetBTagCSV",
#"Lep_JetNDauCharged","Lep_sig3d","Lep_dxy","Lep_dz","Lep_mvaId","Lep_SegCompat",
#"FakeLep_chRelIso","FakeLep_neuRelIso","FakeLep_jetPtRel","FakeLep_jetPtRatio","FakeLep_JetBTagCSV",
#"FakeLep_JetNDauCharged","FakeLep_sig3d","FakeLep_dxy","FakeLep_dz","FakeLep_mvaId","FakeLep_SegCompat",

#"patElectron_genMother_pt","patElectron_genMother_eta","patElectron_genMother_phi","patElectron_genMother_en","patElectron_genMother_pdgId",
#"patElectron_genGrandMother_pt","patElectron_genGrandMother_eta","patElectron_genGrandMother_phi","patElectron_genGrandMother_en","patElectron_genGrandMother_pdgId",
#"Muon_genMother_pt","Muon_genMother_eta","Muon_genMother_phi","Muon_genMother_en","Muon_genMother_pdgId",
#"Muon_genGrandMother_pt","Muon_genGrandMother_eta","Muon_genGrandMother_phi","Muon_genGrandMother_en","Muon_genGrandMother_pdgId",
#"Tau_genMother_pt","Tau_genMother_eta","Tau_genMother_phi","Tau_genMother_en","Tau_genMother_pdgId",
#"Tau_genGrandMother_pt","Tau_genGrandMother_eta","Tau_genGrandMother_phi","Tau_genGrandMother_en","Tau_genGrandMother_pdgId",
#"Jet_genMother_pt","Jet_genMother_eta","Jet_genMother_phi","Jet_genMother_en","Jet_genMother_pdgId",
#"Jet_genGrandMother_pt","Jet_genGrandMother_eta","Jet_genGrandMother_phi","Jet_genGrandMother_en","Jet_genGrandMother_pdgId",
#"Jet_btag_sf","Jet_btag_jesup","Jet_btag_hfup","Jet_btag_hfstat1up","Jet_btag_hfstat2up","Jet_btag_lfup","Jet_btag_lfstat1up","Jet_btag_lfstat2up","Jet_btag_cerr1up","Jet_btag_cerr2up",
#"Jet_btag_jesdown","Jet_btag_hfdown","Jet_btag_hfstat1down","Jet_btag_hfstat2down","Jet_btag_lfdown","Jet_btag_lfstat1down","Jet_btag_lfstat2down","Jet_btag_cerr1down","Jet_btag_cerr2down",


#"Lep_genMother_pt","Lep_genMother_eta","Lep_genMother_phi","Lep_genMother_en","Lep_genMother_pdgId",
#"Lep_genGrandMother_pt","Lep_genGrandMother_eta","Lep_genGrandMother_phi","Lep_genGrandMother_en","Lep_genGrandMother_pdgId",
#"FakeLep_genMother_pt","FakeLep_genMother_eta","FakeLep_genMother_phi","FakeLep_genMother_en","FakeLep_genMother_pdgId",
#"FakeLep_genGrandMother_pt","FakeLep_genGrandMother_eta","FakeLep_genGrandMother_phi","FakeLep_genGrandMother_en","FakeLep_genGrandMother_pdgId",
#"Lep_genMother_pdgId","Lep_genGrandMother_pdgId",
#"FakeLep_genMother_pdgId","FakeLep_genGrandMother_pdgId",
#"Muon_mcPromptGamma",
#"patElectron_mcPromptGamma",
#"Lep_mcPromptGamma",
#"FakeLep_mcPromptGamma",

#"pvertex_ndof","pvertex_z","pvertex_Rho"
#"Muon_mcMatchId","patElectron_mcMatchId","Lep_mcMatchId","FakeLep_mcMatchId",
#"hadW_Gen_Mother_pdgId","hadW_Gen_drjj","hadW_Gen_mass",

#"Gen_numDaught","Gen_BmotherIndices"

"FakeLep_loosejetcsv","FakeLep_loosejetdr"
]

#VariableType  = "int"
#VariableNames = [
# CaseA Variables
#"Muon_soft","Muon_loose","Muon_medium","Muon_tight","Muon_isHighPt",
#"Muon_POGisGood","Muon_pdgId","Muon_pf","Muon_isGlobal","Muon_isTrackerMuon",
#"Muon_tunePBestTrackType","Muon_gen_pdgId","Muon_gen_isPromptFinalState","Muon_gen_isDirectPromptTauDecayProductFinalState"

# CaseB Variables

# CaseC Variables
#]

# Create the variable file
Vectorname = VariableType+Case+"_Vector.cc"
vector     = file (Vectorname,"w")

#ReadTreeptr & WriteTreeptr
RTreeptr = "readingtree"
WTreeptr = "newtree"

#Name of Current Entry
ParEntry = "tentry"
#Name of index in Push_back
ParSel = "gen_en"

###################
### Script itself
#################
if Case == "CaseA":
 print >> vector, "//This is CaseA"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* r"+Variable+"; TBranch* b_r"+Variable+" =0;"
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* "+Variable+" = new std::vector<" + VariableType+">;"
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Clear vector"
 for Variable in VariableNames:
     print >> vector, "    "+Variable+"->clear();"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, "    b_r"+Variable+"->GetEntry("+ParEntry+");"
 
 print >> vector, "//Push_back_variables"
 for Variable in VariableNames:
     print >> vector, "        "+Variable+"->push_back(r"+Variable+"->at("+ParSel+"));"

elif Case == "CaseB":
 print >> vector, "//This is CaseB"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* r"+Variable+"; TBranch* b_r"+Variable+" =0;"
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, " b_r"+Variable+"->GetEntry("+ParEntry+");"
 
elif Case == "CaseC":
 print >> vector, "//This is CaseC"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* "+Variable+" = new std::vector<" + VariableType+">;"
 
 print >> vector, "//source file definition"
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Clear vector"
 for Variable in VariableNames:
     print >> vector, "    "+Variable+"->clear();"

 print >> vector, "//Push_back_fist_variable"
 for Variable in VariableNames:
     print >> vector, "    "+Variable+"->push_back(r"+Variable+"->at(0));"

 print >> vector, "//sort iterator"
 for Variable in VariableNames:
     print >> vector, "     auto it_"+Variable+" = "+Variable+"->begin();"

 print >> vector, "//sort insert"
 for Variable in VariableNames:
     print >> vector, "      "+Variable+"->insert(it_"+Variable+"+index,"+Variable+"->at("+ParSel+"));"

 print >> vector, "//Push_back_variables"
 for Variable in VariableNames:
     print >> vector, "     "+Variable+"->push_back(r"+Variable+"->at("+ParSel+"));"


else:
 print >> vector, "//!!!!Please specify a correct Case!!!!!!//"
