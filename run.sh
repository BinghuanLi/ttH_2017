
if [ "$1" == "v1" ]; then
root -l -b -q Rootplizer_TTHLep_IHEP.cc+'("/home/binghuan/Work/RootTestFiles/TTHLep_2017/data/NTuple/TTHnobb_Ntuple_99.root","TTHnobb_IHEP_99.root")'

else
root -l -b -q Rootplizer_TTHLep_v2IHEP.cc+'("/home/binghuan/Work/RootTestFiles/TTHLep_2017/data/NTuple/TTHnobb_Ntuple_99.root","TTHnobb_rootpla_99.root")'

fi
