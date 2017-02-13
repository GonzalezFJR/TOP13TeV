# Start root 6
#source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh; 
#source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh;
#source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh; 
#source /opt/root6/bin/thisroot.sh
root6

# Start PoD
source /opt/PoD/PoD_env.sh

# Start PAF
source /opt/PAF/PAF_setup.sh

resetpaf

#root -l -b -q 'RunTree_ReReco.C("TestHeppy", 1, 0, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg", 50, 0, true)'
root -l -b -q 'RunTree_ReReco.C("WJetsToLNu_aMCatNLO"    ,  15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("ZZ_ext"   , 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("WW_ext"   , 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("WZ_ext"   , 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TW_ext"   , 30, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TbarW_ext", 30, 0, true)'
#
resetpaf -a
root -l -b -q 'RunTree_ReReco.C("TTWToLNu_ext2"    , 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TTZToQQ"   , 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TTZToLLNuNu_ext", 15, 0, true)'
root -l -b -q 'RunTree_ReReco.C("TTWToQQ"   , 15, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("TTGJets"   , 5, 0, true)'

resetpaf -a
root -l -b -q 'RunTree_ReReco.C("MuonEG"    , 40, 0, true)'
root -l -b -q 'RunTree_ReReco.C("DoubleMuon", 40, 0, true)'
root -l -b -q 'RunTree_ReReco.C("DoubleEG"  , 40, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("SingleElectron"  , 40, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("SingleMuon"  , 40, 0, true)'

resetpaf -a
root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_MLM_ext"  , 50, 0, true)'
root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M5to50_MLM",  30, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("Tree_TTJets_MLM", 30, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("Tree_TTJets_MLM_FastSim", 30, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_aMCatNLO_ext"  , 30, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M10to50_aMCatNLO_ext", 30, 0, true)'

