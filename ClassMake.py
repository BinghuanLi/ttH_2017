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
Case = "CaseB"
#Case = "CaseC"
# Variable Definition

rObject = "Jet"
wObject = "Jet"

VariableType  = "double"
VariableNames = [
# CaseA Variables
"pt","eta","phi","energy",
"pt","eta","phi","energy",
#"charge",
#"IP3Dsig","miniIsoRel","pdgId","gsfTrack_dxy_pv","gsfTrack_dz_pv",


#"jetptratio","jetcsv","lepjetchtrks","miniIsoCh","miniIsoPUsub","ptrel",
#"px","py","pz","jetdr",
"gen_pt","gen_eta","gen_phi","gen_en","gen_pdgId",
"genMother_pt","genMother_eta","genMother_phi","genMother_en","genMother_pdgId",
"genGrandMother_pt","genGrandMother_eta","genGrandMother_phi","genGrandMother_en","genGrandMother_pdgId",
"gen_isPromptFinalState","gen_isDirectPromptTauDecayProductFinalState",

# Electron only
#"SCeta","expectedMissingInnerHits","full5x5_sigmaIetaIeta","hOverE","dEtaIn","dPhiIn","ooEmooP", 
#"isGsfCtfScPixChargeConsistent","isGsfScPixChargeConsistent",

#"isGlobal","chi2","chi2LocalPosition","trkKink","validFraction","segmentCompatibility","pTErrOVpT_it", # Muon only

#Taus
#"packedLeadTauCand_dz","packedLeadTauCand_dxy","byLooseIsolationMVArun2v1DBdR03oldDMwLT","decayModeFinding"
#"byMediumIsolationMVArun2v1DBdR03oldDMwLT",

#Jets
"Uncorr_pt",
"pfCombinedInclusiveSecondaryVertexV2BJetTags","pfCombinedMVAV2BJetTags",
#"px","py","pz","mass",
"qg","axis2","ptD","mult",
"partonFlavour","hadronFlavour","genpt","geneta","genphi","genenergy",

#"JesSF","JesSFup","JesSFdown","JerSF","JerSFup","JerSFdown",
#"neutralHadEnergyFraction","neutralEmEnergyFraction","chargedMultiplicity","numberOfConstituents","chargedHadronEnergyFraction", "chargedEmEnergyFraction",
#"btag_sf",
#"btag_jesup","btag_jesdown",
#"btag_hfup","btag_hfdown",
#"btag_hfstat1up","btag_hfstat1down",
#"btag_hfstat2up","btag_hfstat2down",
#"btag_lfup","btag_lfdown",
#"btag_lfstat1up","btag_lfstat1down",
#"btag_lfstat2up","btag_lfstat2down",
#"btag_cerr1up","btag_cerr1down",
#"btag_cerr2up","btag_cerr2down",

#Case B

# Case C

#"cut",

#"BDT","isMedium_ST","corrpt","FR","CF",
#"passConversion","passMuTightCharge","passEleTightCharge","passMissHit","isMatchRightCharge",


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
Vectorname = VariableType+Case+"_Class.cc"
vector     = file (Vectorname,"w")

#ReadTreeptr & WriteTreeptr
RTreeptr = "readingtree"
WTreeptr = "newtree"

#Name of Current Entry
ParEntry = "tentry"
#Name of index in Push_back
ParSel = "tau_en"
ParWrite = "tau_en"

###################
### Script itself
#################
if Case == "CaseA":
 print >> vector, "//This is CaseA"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* r"+rObject+"_"+Variable+"; TBranch* b_r"+rObject+"_"+Variable+" =0;"
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* "+wObject+"_"+Variable+" = new std::vector<" + VariableType+">;"
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+RTreeptr+'->SetBranchAddress("'+rObject+"_"+Variable+'",&r'+rObject+"_"+Variable+",&b_r"+rObject+"_"+Variable+");"
 
 print >> vector, "   //Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+WTreeptr+'->Branch("'+wObject+"_"+Variable+'",&'+wObject+"_"+Variable+");"
 
 print >> vector, "   //Clear vector"
 for Variable in VariableNames:
     print >> vector, "    "+wObject+"_"+Variable+"->clear();"
 
 print >> vector, "   //GetEntry"
 for Variable in VariableNames:
     print >> vector, "    b_r"+rObject+"_"+Variable+"->GetEntry("+ParEntry+");"
 
 print >> vector, "   //class member"
 for Variable in VariableNames:
     print >> vector, "        "+VariableType+" "+Variable+" = -999.;"
 
 print >> vector, "   //Intialize variables"
 for Variable in VariableNames:
     print >> vector, "        jet."+Variable+"= r"+rObject+"_"+Variable+"->at("+ParSel+");"

 print >> vector, "   //Write variables"
 for Variable in VariableNames:
     print >> vector, "        "+wObject+"_"+Variable+"->push_back(jets->at("+ParWrite+")."+Variable+");"

elif Case == "CaseB":
 print >> vector, "//This is CaseB"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* r"+rObject+"_"+Variable+"; TBranch* b_r"+rObject+"_"+Variable+" =0;"
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+RTreeptr+'->SetBranchAddress("'+rObject+"_"+Variable+'",&r'+rObject+"_"+Variable+",&b_r"+rObject+"_"+Variable+");"
 
 print >> vector, "   //GetEntry"
 for Variable in VariableNames:
     print >> vector, "    b_r"+rObject+"_"+Variable+"->GetEntry("+ParEntry+");"
 
 print >> vector, "   //class member"
 for Variable in VariableNames:
     print >> vector, "        "+VariableType+" "+Variable+" = -999;"
 
 print >> vector, "   //Intialize variables"
 for Variable in VariableNames:
     print >> vector, "        "+rObject+"."+Variable+"= r"+rObject+"_"+Variable+"->at("+ParSel+");"

elif Case == "CaseC":
 print >> vector, "//This is CaseC"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector, "vector<"+VariableType+">* "+wObject+"_"+Variable+" = new std::vector<" + VariableType+">;"
 
 print >> vector, "//source file definition"
 
 print >> vector, "   //Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, "    "+WTreeptr+'->Branch("'+wObject+"_"+Variable+'",&'+wObject+"_"+Variable+");"
 
 print >> vector, "   //Clear vector"
 for Variable in VariableNames:
     print >> vector, "    "+wObject+"_"+Variable+"->clear();"
 
 print >> vector, "   //class member"
 for Variable in VariableNames:
     print >> vector, "        "+VariableType+" "+Variable+" = -999;"
 
 print >> vector, "   //Intialize variables"
 for Variable in VariableNames:
     print >> vector, "        "+rObject+"."+Variable+"= r"+rObject+"_"+Variable+"->at("+ParSel+");"

 print >> vector, "   //Write variables"
 for Variable in VariableNames:
     print >> vector, "        "+wObject+"_"+Variable+"->push_back(jets->at("+ParWrite+")."+Variable+");"
