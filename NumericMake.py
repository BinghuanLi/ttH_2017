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
#### CaseA Varialbes
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
#"Muon_gen_pt","Muon_gen_eta","Muon_gen_phi","Muon_gen_en"
# CaseB Variables

# CaseC Variable
#"Top_fGen_pt","Top_fGen_eta","Top_fGen_phi","Top_fGen_energy",
#"Tbar_fGen_pt","Tbar_fGen_eta","Tbar_fGen_phi","Tbar_fGen_energy",
#"H0_fGen_pt","H0_fGen_eta","H0_fGen_phi","H0_fGen_energy",
#"Xguds_fGen_pt","Xguds_fGen_eta","Xguds_fGen_phi","Xguds_fGen_energy",
#"num_Gen_Lvh","num_Gen_Lvtop",
#"EVENT_genWeight"
#"num_hJet","num_htopb","num_ltopb","num_htopW","num_otherj"
#"Gen_type1PF_Metpx","Gen_type1PF_Metpy","Gen_type1PF_Metpz","Gen_type1PF_Meteta","Gen_type1PF_Metphi",
#"Gen_type1PF_Meten",
#"lepControl_dilepttbar",
#"lepControl_semilepttbar",
#"lepControl_zjets",
#"PUWeight"

#"2lSubCat","3lSubCat"
#"leadingCorrpt","subleadingCorrpt","thirdCorrpt",
#"leadingJetdr","subleadingJetdr","thirdJetdr",
#"massll","minMass_OSAF","Sum2lCharge","Sum3lCharge",

#"BWeight","BWeightLFup","BWeightLFdown","BWeightHFup","BWeightHFdown",
#"BWeightHFStats1up","BWeightHFStats1down","BWeightLFStats1up","BWeightLFStats1down","BWeightHFStats2up",
#"BWeightHFStats2down","BWeightLFStats2up","BWeightLFStats2down","BWeightCErr1up","BWeightCErr1down",
#"BWeightCErr2up","BWeightCErr2down","BWeightJESup","BWeightJESdown",

#"isGoodVtx"
#"EVENT_filterBadGlobalMuonTagger","EVENT_filtercloneGlobalMuonTagger"
#"DiLepSSJetCR","DiLepMVAJetCR","DiLepOSJetCR",
#"2lep_bestMVA","2lep_worseMVA","2lep_flavor","2lep_htllv","2lep_mtWmin",
#"2lep_nTight","leadJetCSV","secondJetCSV","thirdJet_CSV","fourthJet_CSV",
#"HighestJetCSV","BadMuonEvent","HtJet","nLepFO","nLepTight",
#"minMllAFAS","minMllAFOS","minMllSFOS",
#"leadLep_pt","leadLep_eta","leadLep_charge","leadLep_ptratio","leadLep_BDT",
#"leadLep_segment","leadLep_dxy","leadLep_dz","leadLep_mvaId","leadLep_jetcsv",
#"leadLep_ptrel","leadLep_minIso","leadLep_minIsoCh","leadLep_minIsoNeu","leadLep_pdgId",
#"secondLep_pt","secondLep_eta","secondLep_charge","secondLep_ptratio","secondLep_BDT",
#"secondLep_segment","secondLep_dxy","secondLep_dz","secondLep_mvaId","secondLep_jetcsv",
#"secondLep_ptrel","secondLep_minIso","secondLep_minIsoCh","secondLep_minIsoNeu","secondLep_pdgId",
#"thirdLep_pt","thirdLep_eta","thirdLep_charge","thirdLep_ptratio","thirdLep_BDT",
#"thirdLep_segment","thirdLep_dxy","thirdLep_dz","thirdLep_mvaId","thirdLep_jetcsv",
#"thirdLep_ptrel","thirdLep_minIso","thirdLep_minIsoCh","thirdLep_minIsoNeu","thirdLep_pdgId",
#"leadLep_phi","secondLep_phi","thirdLep_phi",

#"bjet_hadTop_pt","wjet1_hadTop_pt","wjet2_hadTop_pt","bjet_hadTop_csv","bjet_lepTop_csv","reco_hadTop_mass","reco_hadTop_pt","reco_WhadTop_mass","PtRatio_leptOverleph","Dr_lept_bfromlTop","Dr_lept_bfromhTop","Dr_leph_bfromlTop",

#"minljj_pt","minljj_eta","minljj_phi","minljj_energy","minljj_mass",
#"jj_pt","jj_eta","jj_phi","jj_energy","jj_mass",
#"Mtrue","minlep_dphiMet",

#"W_numMedium","W_numSoft",
#"MediumW_numMatch","SoftW_numMatch","LooseW_numMatch","TightW_numMatch",
#"MediumW_numMatch_fromH","SoftW_numMatch_fromH","LooseW_numMatch_fromH","TightW_numMatch_fromH",
#"MediumW_numMatch_fromNonH","SoftW_numMatch_fromNonH","LooseW_numMatch_fromNonH","TightW_numMatch_fromNonH",
#"hadW_numGen_cone8"

# ttH evt weight
#"tthWeight_SR","tthWeight_OS","tthWeight_DiLepMVA","tthWeight_TriLepMVA"
"isDiMuMVAAR","isDiEleMVAAR","isEleMuMVAAR","isTriLepMVAAR","isQuaLepMVAAR",
"isDiMuOSAR","isDiEleOSAR","isEleMuOSAR",
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
Vectorname = VariableType+Case+"_Numeric.cc"
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
 print >> vector, "//This is Case A"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector,VariableType+" r"+Variable+"; TBranch* b_r"+Variable+" =0;" 
 
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector,VariableType+" "+Variable+";" 
 
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Initialize Number"
 for Variable in VariableNames:
     print >> vector, " "+Variable+"= -999;"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, " b_r"+Variable+"->GetEntry("+ParEntry+");"
 
 print >> vector, "//Write_variables"
 for Variable in VariableNames:
     print >> vector, "  "+Variable+" = r"+Variable+"->at(0);"

elif Case == "CaseB":
 print >> vector, "//This is Case B"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector,VariableType+" r"+Variable+"; TBranch* b_r"+Variable+" =0;" 
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, " b_r"+Variable+"->GetEntry("+ParEntry+");"

elif Case == "CaseC":
 print >> vector, "//This is Case C"
 print >> vector, "//Head file declaration"
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector,VariableType+" "+Variable+";" 
 
 
 print >> vector, "//source file definition"
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Initialize Number"
 for Variable in VariableNames:
     print >> vector, "    "+Variable+"= -999;"
 
 print >> vector, "//Write_variables"
 for Variable in VariableNames:
     print >> vector, "    "+Variable+" = r"+Variable+";"
else :
 print >> vector, "//!!!!!!Please specify a correct Case!!!!!!//"
