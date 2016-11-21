# Start root 6
#source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh;
#source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh;
#source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh;
#source /opt/root6/bin/thisroot.sh
root6
# Start PoD
source /opt/PoD/PoD_env.sh
# Start PAF
#source /nfs/fanae/PAF_releases/head/PAF_setup.sh
source /opt/PAF/PAF_setup.sh

resetpaf

root -l -b -q 'RunStopAnalysis.C("TTbar_Powheg_ext", 50, true)'
root -l -b -q 'RunStopAnalysis.C("WJetsToLNu_aMCatNLO"    ,  8, true)'
root -l -b -q 'RunStopAnalysis.C("ZZ"   , 8, true)'
root -l -b -q 'RunStopAnalysis.C("WW"   , 8, true)'
root -l -b -q 'RunStopAnalysis.C("WZ"   , 8, true)'
root -l -b -q 'RunStopAnalysis.C("TW"   , 30, true)'
root -l -b -q 'RunStopAnalysis.C("TbarW", 30, true)'
#
root -l -b -q 'RunStopAnalysis.C("TTWToLNu"    , 8, true)'
root -l -b -q 'RunStopAnalysis.C("TTZToQQ"   , 8, true)'
root -l -b -q 'RunStopAnalysis.C("TTZToLLNuNu"   , 8, true)'
root -l -b -q 'RunStopAnalysis.C("TTWToQQ"   , 8, true)'

#root -l -b -q 'RunStopAnalysis.C("TTGJets"   , 5, true)'

root -l -b -q 'RunStopAnalysis.C("MuonEG"    , 20, true)'
root -l -b -q 'RunStopAnalysis.C("DoubleMuon", 20, true)'
root -l -b -q 'RunStopAnalysis.C("DoubleEG"  , 20, true)'
root -l -b -q 'RunStopAnalysis.C("SingleElectron"  , 20, true)'
root -l -b -q 'RunStopAnalysis.C("SingleMuon"  , 20, true)'


root -l -b -q 'RunStopAnalysis.C("DYJetsToLL_M50_aMCatNLO"  , 30, true)'
root -l -b -q 'RunStopAnalysis.C("DYJetsToLL_M10to50_aMCatNLO",  30, true)'
#root -l -b -q 'RunStopAnalysis.C("DYJetsToLL_M50_aMCatNLO_ext"  , 30, true)'
#root -l -b -q 'RunStopAnalysis.C("DYJetsToLL_M10to50_aMCatNLO_ext", 30, true)'
source runStop_signal.sh


#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M1",     1, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M10",    1, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M50",    4, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M100",   4, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M200",   4, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M500",   4, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M600",   4, true)'
#root -l -b -q 'RunStopAnalysis.C("TTDMJets_EFT_M1000",  4, true)'
