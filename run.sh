
if [ "$1" == "v1" ]; then
root -b -q -l Rootplizer_TTHLep_IHEP.cc+'("/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/FullMorV2_ttHnobb/170530_161519/0000//OutTree_2.root","OutTree_2_rootplas.root")'

elif [ "$1" == "data" ]; then
root -b -q -l Rootplizer_TTHLep_v2IHEP.cc+'("/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/data/SingleMuon/crab_FullMorV2_SMuBlockF1_recover/170607_132915/0000/OutTree_90.root","RSMu_90.root")'

else
root -b -q -l Rootplizer_TTHLep_v2IHEP.cc+'("/publicfs/cms/data/TopQuark/cms13TeV/FullMorV2/mc/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/FullMorV2_ttHnobb/170530_161519/0000//OutTree_2.root","OutTree_2_v2rootplas.root")'
fi
