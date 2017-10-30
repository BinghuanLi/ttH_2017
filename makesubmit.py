import re
import sys
import os
import glob
import string
#####
##   Parameters to be specified by the user
#####
#analysis and task
analysis = "TTHLep"
taskname = "v2IHEP"
#for the queue
#workpath    = "/publicfs/cms/user/libh/Submit_Condor/"+analysis+"/Sample80X"
workpath    = "/publicfs/cms/user/libh/Submit_Condor/"+analysis+"/BoostJet"
jobDir      = workpath+"/"+"Jobs"
AnalyzerDir = workpath+"/"+"Analyzer"
task        = analysis+"_"+taskname
rootplizer  = "Rootplizer_"+task+".cc"
headplizer  = "Rootplizer_"+task+".h"
#Directory of input files
sample={
#Data
#"SEleB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockB1/170531_111453/0000/',#1758
##"SEleB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockB1/170531_111453/0001/',#1758
#"SEleC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockC1/170531_112614/0000/',#580
#"SEleD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockD1/170531_112857/0000/',#972
#"SEleE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockE1/170531_113137/0000/',#826
#"SEleF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockF1/170531_113418/0000/',
#"SEleF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockF2/170531_113728/0000/',
#"SEleG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockG1/170531_114004/0000/',
##"SEleG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockG1/170531_114004/0001/',
#"SEleH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockH1/170531_114244/0000/',
##"SEleH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockH1/170531_114244/0001/',
#"SEleH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockH2/170530_164040/0000/',#41
#
#"SMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockB1/170531_114536/0000/',
##"SMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockB1/170531_114536/0001/',
#"SMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockC1/170531_114815/0000/',
#"SMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockD1/170531_115112/0000/',
#"SMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockE1/170531_115354/0000/',
#"SMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockF1/170531_115635/0000/',
#"SMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockF2/170531_115917/0000/',
#"SMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockG1/170531_120200/0000/',
##"SMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockG1/170531_120200/0001/',
#"SMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockH1/170531_120434/0000/',
##"SMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockH1/170531_120434/0001/',
#"SMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockH2/170531_120710/0000/',
#
#"DEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockB1/170531_120945/0000/',
##"DEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockB1/170531_120945/0001/',
#"DEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockC1/170531_121222/0000/',
#"DEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockD1/170531_121509/0000/',
#"DEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockE1/170531_121746/0000/',
#"DEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockF1/170531_122029/0000/',
#"DEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockF2/170531_122307/0000/',
#"DEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockG1/170531_122559/0000/',
##"DEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockG1/170531_122559/0001/',
#"DEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockH1/170531_122837/0000/',
##"DEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockH1/170531_122837/0001/',
#"DEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockH2/170531_123117/0000/',
#
################ There is a mis assign in multicrab_DATA between DoubleMuon and MuonEG #############
#"DMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockB1/170531_125358/0000/',
##"DMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockB1/170531_125358/0001/',
#"DMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockC1/170531_125627/0000/',
#"DMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockD1/170531_125857/0000/',
#"DMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockE1/170531_130126/0000/',
#"DMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockF1/170531_130344/0000/',
#"DMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockF2/170531_130606/0000/',
#"DMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockG1/170531_130824/0000/',
##"DMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockG1/170531_130824/0001/',
#"DMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockH1/170531_131041/0000/',
##"DMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockH1/170531_131041/0001/',
#"DMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockH2/170531_131302/0000/',
#
#"MuEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockB1/170531_123341/0000/',
##"MuEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockB1/170531_123341/0001/',
#"MuEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockC1/170531_123555/0000/',
#"MuEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockD1/170531_123810/0000/',
#"MuEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockE1/170531_124023/0000/',
#"MuEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockF1/170531_124239/0000/',
#"MuEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockF2/170531_124449/0000/',
#"MuEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockG1/170531_124700/0000/',
##"MuEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockG1/170531_124700/0001/',
#"MuEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockH1/170531_124911/0000/',
##"MuEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockH1/170531_124911/0001/',
#"MuEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockH2/170531_125134/0000/',
#
#### Data Recover ####
#"RSEleB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockB1_recover/170607_125832/0000/', #290
#"RSEleC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockC1_recover/170607_130129/0000/', #120
#"RSEleD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockD1_recover/170607_130415/0000/', #222
#"RSEleE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockE1_recover/170607_130651/0000/', #195
#"RSEleF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockF1_recover/170607_130939/0000/', #160
#"RSEleF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockF2_recover/170607_131218/0000/', #9
#"RSEleG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockG1_recover/170607_131452/0000/', #383
#"RSEleH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleElectron/crab_FullMorV2_SEleBlockH1_recover/170607_131714/0000/', #245
#
#"RSMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockB1_recover/170607_131929/0000/', 
#"RSMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockC1_recover/170607_132148/0000/',
#"RSMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockD1_recover/170607_132411/0000/',
#"RSMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockE1_recover/170607_132644/0000/',
#"RSMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockF1_recover/170607_132915/0000/',
#"RSMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockF2_recover/170607_133134/0000/',
#"RSMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockG1_recover/170607_133354/0000/',
#"RSMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockH1_recover/170607_133623/0000/',
#"RSMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockH2_recover/170607_133844/0000/',
#
#"RDEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockB1_recover/170607_134111/0000/',
#"RDEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockC1_recover/170607_134332/0000/',
#"RDEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockD1_recover/170607_134624/0000/',
#"RDEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockE1_recover/170607_134958/0000/',
#"RDEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockF1_recover/170607_182838/0000/',
#"RDEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockF2_recover/170607_140245/0000/',
#"RDEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockG1_recover/170607_140855/0000/',
#"RDEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockH1_recover/170607_141422/0000/',
#"RDEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleEG/crab_FullMorV2_DblEGBlockH2_recover/170607_142119/0000/',
#
#"RDMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockB1_recover/170607_145709/0000/',
#"RDMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockC1_recover/170607_145932/0000/',
#"RDMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockD1_recover/170607_150245/0000/',
#"RDMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockE1_recover/170607_150602/0000/',
#"RDMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockF1_recover/170607_150850/0000/',
#"RDMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockF2_recover/170607_151148/0000/',
#"RDMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockG1_recover/170607_151446/0000/',
#"RDMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockH1_recover/170607_151806/0000/',
#"RDMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/DoubleMuon/crab_FullMorV2_MuEGBlockH2_recover/170607_152113/0000/',
#
#"RMuEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockB1_recover/170607_142447/0000/',
#"RMuEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockC1_recover/170607_143737/0000/',
#"RMuEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockD1_recover/170607_144007/0000/',
#"RMuEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockE1_recover/170607_144236/0000/',
#"RMuEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockF1_recover/170607_144503/0000/',
#"RMuEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockF2_recover/170607_144730/0000/',
#"RMuEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockG1_recover/170607_144953/0000/',
#"RMuEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockH1_recover/170607_145205/0000/',
#"RMuEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/MuonEG/crab_FullMorV2_DblMuBlockH2_recover/170607_145442/0000/',

#Mc
##samples not used in final ttHleptons
#"TT_Tune":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/FullMorV2_TT/170531_195914/0000/',
##samples used in final ttHleptons
##Sig
"TTHnobb":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/FullMorV2_ttHnobb/170530_161519/0000/',#'FullMorV1_ttHnobb', 108
###Bkg
"TTWToLNuext2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/FullMorV2_amcTTWJetsToLNuext2/170531_182459/0000/',
"TTWToLNuext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/FullMorV2_amcTTWJetsToLNuext1/170531_182712/0000/',
"TTZToLLNuNu":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_amcTTZToLLNuNu_M-10_ext1/170531_182920/0000/',
"TTZToLL_M1to10":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTZToLL_M1to10/170531_183129/0000/',

"TTJets_sinLepTbar_v1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_sinLepTbar/170531_190717/0000/',
"TTJets_sinLepTbar_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_sinLepTbar_ext1/170531_190925/0000/',
"TTJets_sinLepT_v1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_sinLepT/170531_191441/0000/',
"TTJets_sinLepT_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_sinLepT_ext1/170531_191648/0000/',
"TTJets_diLep_v1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_diLep/170531_191854/0000/',
"TTJets_diLep_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_TTJets_diLep_ext1/170531_192103/0000/',
#
#"TTGJets_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/FullMorV2_TTGJets_ext1/170621_150600/0000/',#239 files
#"WGToLNuG_ext2":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_WGToLNuG_ext2/170621_150118/0000/', #124 files
#"TGJets_v1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/FullMorV2_TGJets/170531_184103/0000/',
#"WGToLNuG":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_WGToLNuG_ext1/170614_105550/0000/', #36G
#"ZGTo2LG":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_ZGTo2LG_ext1/170614_105812/0000/', #124G
#"TGJets_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/FullMorV2_TGJets_ext1/170614_110340/0000/',#20G
#
#
#"WpWpJJ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/FullMorV2_WpWpJJ_EWK-QCD/170614_110838/0000/',#2.3G
#"WWTo2L2Nu_DS":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/FullMorV2_WWTo2L2Nu_DS/170614_111109/0000/',#8.5G
#"WWW_4F":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WWW_4F/170614_111337/0000/',#3.0G
#"WWZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WWZ/170614_111630/0000/',#3.3G
#"WZZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WZZ/170614_111852/0000/',#2.5G
#"ZZZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_ZZZ/170614_112130/0000/',#3.1G
#"tZq":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/tZq_ll_4f_13TeV-amcatnlo-pythia8/FullMorV2_tZq_ext1/170614_112447/0000/',#195G
#"TTTT":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_TTTT/170614_112737/0000/',#5.2G
#"tWll":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/FullMorV2_tWll/170614_112953/0000/',#1G
#"amcWJets":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_amcWJets/170614_113235/0000/',#179G
#"WZTo3LNu":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/FullMorV2_WZTo3LNu/170614_113533/0000/',#21G
#"WWTo2L2Nu":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWTo2L2Nu_13TeV-powheg/FullMorV2_WWTo2L2Nu/170614_113820/0000/',#21G
#"ZZTo4L":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZZTo4L_13TeV_powheg_pythia8/FullMorV2_ZZTo4L/170614_114133/0000/',#69G
#
#
#
#"powhegST":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/FullMorV2_powhegST/170531_192315/0000/',
#"powhegSaT":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/FullMorV2_powhegSaT/170531_192520/0000/',
#"powhegSTt":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/FullMorV2_powhegSTt/170531_192724/0000/',
#"powhegSTat":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/FullMorV2_powhegSTat/170531_193925/0000/',
#"amcSTs":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/FullMorV2_amcSTs/170531_194132/0000/',
#"DY_M10to50":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_DY_M10to50/170531_194340/0000/',
#"DY_M50_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_DY_M50_ext1-v2/170707_103414/0000/',
#"madWJets":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/FullMorV2_WJetsToLNu/170707_103639/0000/',
#
############################# OW STAT ##############
#"WGToLNuG":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_WGToLNuG_ext1/170531_183654/0000/',
#"ZGTo2LG":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_ZGTo2LG_ext1/170531_183859/0000/',
#"TGJets_v1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/FullMorV2_TGJets/170531_184103/0000/',
#"TGJets_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/FullMorV2_TGJets_ext1/170531_184315/0000/',
#"TTGJets":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/FullMorV2_TTGJets/170531_184600/0000/',
#"WpWpJJ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/FullMorV2_WpWpJJ_EWK-QCD/170531_184917/0000/',
#"WWTo2L2Nu_DS":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/FullMorV2_WWTo2L2Nu_DS/170531_185141/0000/',
#"WWW_4F":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WWW_4F/170531_185359/0000/',
#"WWZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WWZ/170531_185615/0000/',
#"WZZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_WZZ/170531_182020/0000/',
#"ZZZ":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_ZZZ/170531_190035/0000/',
#"tZq":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/tZq_ll_4f_13TeV-amcatnlo-pythia8/FullMorV2_tZq_ext1/170531_190247/0000/',
#"TTTT":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/FullMorV2_TTTT/170531_190457/0000/',
#"amcWJets":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/FullMorV2_amcWJets/170531_194800/0000/',
#"WZTo3LNu":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/FullMorV2_WZTo3LNu/170531_195006/0000/',
#"WWTo2L2Nu":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/WWTo2L2Nu_13TeV-powheg/FullMorV2_WWTo2L2Nu/170531_195212/0000/',

#"ZZTo4L":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ZZTo4L_13TeV_powheg_pythia8/FullMorV2_ZZTo4L/170531_195459/0000/',

#"TTGJets":'/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/FullMorV2_TTGJets/170614_110604/0000/',#42G
       }
sampleout={
##data
"SEleB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleB1',
"SEleC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleC1',
"SEleD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleD1',
"SEleE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleE1',
"SEleF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleF1',
"SEleF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleF2',
"SEleG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleG1',
"SEleH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleH1',
"SEleH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SEleH2',
"SMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuB1',
"SMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuC1',
"SMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuD1',
"SMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuE1',
"SMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuF1',
"SMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuF2',
"SMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuG1',
"SMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuH1',
"SMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/SMuH2',
"DEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmB1',
"DEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmC1',
"DEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmD1',
"DEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmE1',
"DEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmF1',
"DEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmF2',
"DEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmG1',
"DEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmH1',
"DEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DEleGmH2',
"MuEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmB1',
"MuEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmC1',
"MuEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmD1',
"MuEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmE1',
"MuEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmF1',
"MuEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmF2',
"MuEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmG1',
"MuEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmH1',
"MuEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/MuEleGmH2',
"DMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuB1',
"DMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuC1',
"DMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuD1',
"DMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuE1',
"DMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuF1',
"DMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuF2',
"DMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuG1',
"DMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuH1',
"DMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DMuH2',

##data recover
"RSEleB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleB1', 
"RSEleC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleC1',
"RSEleD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleD1',
"RSEleE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleE1',
"RSEleF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleF1',
"RSEleF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleF2',
"RSEleG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleG1',
"RSEleH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleH1',
"RSEleH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSEleH2',
"RSMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuB1',
"RSMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuC1',
"RSMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuD1',
"RSMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuE1',
"RSMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuF1',
"RSMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuF2',
"RSMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuG1',
"RSMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuH1',
"RSMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RSMuH2',
"RDEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmB1',
"RDEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmC1',
"RDEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmD1',
"RDEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmE1',
"RDEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmF1',
"RDEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmF2',
"RDEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmG1',
"RDEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmH1',
"RDEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDEleGmH2',
"RMuEleGmB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmB1',
"RMuEleGmC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmC1',
"RMuEleGmD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmD1',
"RMuEleGmE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmE1',
"RMuEleGmF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmF1',
"RMuEleGmF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmF2',
"RMuEleGmG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmG1',
"RMuEleGmH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmH1',
"RMuEleGmH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RMuEleGmH2',
"RDMuB1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuB1',
"RDMuC1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuC1',
"RDMuD1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuD1',
"RDMuE1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuE1',
"RDMuF1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuF1',
"RDMuF2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuF2',
"RDMuG1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuG1',
"RDMuH1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuH1',
"RDMuH2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/RDMuH2',

#Sig
"TT_Tune":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TT_Tune',

"TTHnobb":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTHnobb',

"TTWToLNuext2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTWToLNuext2',
"TTWToLNuext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTWToLNuext1',
"TTZToLLNuNu":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTZToLLNuNu',
"TTZToLL_M1to10":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTZToLL_M1to10',
"WGToLNuG":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WGToLNuG',
"WGToLNuG_ext2":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WGToLNuG_ext2',
"ZGTo2LG":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/ZGTo2LG',
"TGJets_v1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TGJets_v1',
"TGJets_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TGJets_ext1',
"TTGJets":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTGJets',
"TTGJets_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTGJets_ext1',
"WpWpJJ":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WpWpJJ',
"WWTo2L2Nu_DS":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WWTo2L2Nu_DS',
"WWW_4F":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WWW_4F',
"WWZ":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WWZ',
"WZZ":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WZZ',
"ZZZ":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/ZZZ',
"tZq":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/tZq',
"TTTT":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTTT',
"tWll":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/tWll',
"TTJets_sinLepTbar_v1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_sinLepTbar_v1',
"TTJets_sinLepTbar_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_sinLepTbar_ext1',
"TTJets_sinLepT_v1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_sinLepT_v1',
"TTJets_sinLepT_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_sinLepT_ext1',
"TTJets_diLep_v1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_diLep_v1',
"TTJets_diLep_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/TTJets_diLep_ext1',
"powhegST":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/powhegST',
"powhegSaT":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/powhegSaT',
"powhegSTt":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/powhegSTt',
"powhegSTat":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/powhegSTat',
"amcSTs":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/amcSTs',
"DY_M10to50":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DY_M10to50',
"DY_M50_ext1":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/DY_M50_ext1',
"amcWJets":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/amcWJets',
"madWJets":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/madWJets',
"WZTo3LNu":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WZTo3LNu',
"WWTo2L2Nu":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/WWTo2L2Nu',
"ZZTo4L":'/publicfs/cms/data/TopQuark/cms13TeV/Binghuan/'+analysis+'/FullMorV2/ZZTo4L',


          }
#####
##   The script itsself
#####
cshFilePath = jobDir+"/"+"sh"
logFilePath = jobDir+"/"+"log"
if os.path.exists(jobDir):
	os.popen('rm -fr '+jobDir)
if os.path.exists(AnalyzerDir):
        os.popen('rm -fr '+AnalyzerDir)
os.popen('mkdir -p '+cshFilePath)
os.popen('mkdir -p '+logFilePath)

allJobFileName = "all.sh"
allJobFile      = file(allJobFileName,"w")
print >> allJobFile, "#!/bin/bash"
print >> allJobFile, "cd ",cshFilePath
#print >> allJobFile, "chmod +x "+cshFilePath+"/*.sh"

MergeFileName = "merge.sh"
MergeFile      = file(MergeFileName,"w")
MergeSourceFile = " "
def prepareSubmitJob(submitFileName,cshFileName):
	cshFile      = file(submitFileName,"w")
	print >> cshFile, "Universe     = vanilla"
	print >> cshFile, "getenv       = true"
	print >> cshFile, "accounting_group    = cms"
	print >> cshFile, "Executable   = ",cshFileName
	#print >> cshFile, "Output       = ",outPutFileName
	#print >> cshFile, "Error        = ",errorFileName
	print >> cshFile, "Queue"

def prepareCshJob(input,output,submitFileName,analyzerpath,stdOutFile,stdErrFile):
        subFile =file(submitFileName,"w")
        os.popen('chmod +x '+submitFileName)
        print >> subFile, "#!/bin/bash"
        print >> subFile, "/bin/hostname"
        #print >> subFile, "cd /cvmfs/cms.cern.ch"
        #print >> subFile, "source cmsset_default.csh"
        #print >> subFile, "cd slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_3"
        #print >> subFile, "cmsenv"
        #print >> subFile, "cd -"
        print >> subFile, "gcc -v"
        #print >> subFile, "root -b -v -q"
        print >> subFile, "pwd"
	#print >> subFile, "cd /publicfs/cms/data/TopQuark/cms13TeV/software/root/bin/"
	#print >> subFile, "source thisroot.csh"
	#print >> subFile, "cd /publicfs/cms/user/libh/CMSSW_5_3_9/src/ttH_13Tev"
	#print >> subFile, "setenv SCRAM_ARCH slc5_amd64_gcc462"
	#print >> subFile, "source /cvmfs/cms.cern.ch/cmsset_default.csh"
	#print >> subFile, "source  /afs/ihep.ac.cn/soft/CMS/64bit/root/profile/rootenv-entry 5.34.18"
        #print >> subFile, "eval \`scramv1 runtime -sh\`"
        print >> subFile, "cd "+analyzerpath
	#print >> subFile, "cp ${jobDir}/getAhist.C ."
        print >> subFile, "root -b -q -l "+rootplizer+"+'(\""+input+"\",\""+output+"\")' 2>"+ stdErrFile+">" + stdOutFile
        #print >> subFile, "root -b -q -l "+rootplizer+"'(\""+input+"\",\""+output+"\")'"
 
#for iroot in range(nroot):
N_totjobs =0
for k in sample:
	print k
	print sample[k]
	sampleName = k
	rootDirectory   = sample[k] # +"/"+k    ## you can use your naming convension to set up the inputDirectory
	outputDirectory = sampleout[k] # you can use your naming convension to set up the outputDirectory
	AnalyzerSampleDir = AnalyzerDir + "/" + sampleName
        os.popen('mkdir -p ' +AnalyzerSampleDir)
        os.popen('mkdir -p '+outputDirectory)
        os.chdir(rootDirectory)
	roots = glob.glob('*.root')
	roots.sort()
	nroot = len(roots)
	N_totjobs += nroot
	print str(nroot)
	for iroot in range(nroot):
	#for iroot in range(10):
		#print rootDirectory+"/"+roots[iroot]
		input  = rootDirectory+"/"+roots[iroot]
		output = outputDirectory+"/"+roots[iroot].replace(".root","_rootplas.root")
	        analyzerpath = AnalyzerSampleDir+"/"+roots[iroot]
                os.popen('mkdir -p '+analyzerpath)
                #os.system('pwd ')
                command_cp_cc = 'cp '+workpath+"/"+rootplizer+" "+analyzerpath
                command_cp_h = 'cp '+workpath+"/"+headplizer+" "+analyzerpath
                os.system(command_cp_cc)
                #os.system(command_cp_h)
                MergeSourceFile += " "+roots[iroot].replace(".root","_rootplas.root")	
                #make the name of submitFile and InputFile matched
                name = roots[iroot].replace('.root','')
                
                submitFileName = sampleName+"_"+name+".submit"

		
		result = re.search('OutTree_(.*).root', roots[iroot])
                S_result = result.group(1)
                index = int(S_result)-1
	
		#cshFileName = cshFilePath+"/"+sampleName+".sh."+str(index)
		cshFileName = cshFilePath+"/"+sampleName+".sh."+str(iroot)
		#cshFileName = cshFilePath+"/"+sampleName+".sh."+str(S_result)
		logFileName = logFilePath+"/"+sampleName+"_"+name+".log"
		errFileName = logFilePath+"/"+sampleName+"_"+name+".err"
	
		#cshFileName = cshFilePath+"/"+sampleName+"_"+name+".sh"
		#logFileName = logFilePath+"/"+sampleName+"_"+name+".log"
		#errFileName = logFilePath+"/"+sampleName+"_"+name+".err"
	
		#prepareSubmitJob(cshFilePath+"/"+submitFileName, cshFileName, logFileName, errFileName)
		#prepareCshJob(input,output,cshFileName,analyzerpath)
		#prepareSubmitJob(cshFilePath+"/"+submitFileName, cshFileName)
		prepareCshJob(input,output,cshFileName,analyzerpath, logFileName, errFileName)
		
		#print >> allJobFile, "condor_submit "+ submitFileName + " -group cms -name job@schedd01.ihep.ac.cn"	
		#print >> allJobFile, "hep_sub "+ cshFileName + " -g cms"
	print >> allJobFile, "hep_sub "+ sampleName + ".sh.%{ProcId} -n "+str(nroot)+" -g cms"

print "number of total jobs are ",N_totjobs
print >> MergeFile, "cd",outputDirectory
print >> MergeFile, "hadd Merged_rootplas.root",MergeSourceFile
print >> allJobFile, "cd -"

