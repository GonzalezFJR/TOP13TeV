# Start root 6
source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh; 
source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh;
source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh; 
source /opt/root6/bin/thisroot.sh
# Start PoD
source /opt/PoD/PoD_env.sh
# Start PAF
source /nfs/fanae/PAF_releases/head/PAF_setup.sh

resetpaf

#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg", 1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg", 1, 1, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("WW"          , 1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("WZ"          , 1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("TW"          , 1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("TbarW"       , 1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("DYJetsToLL_M50_aMCatNLO",  4, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("DYJetsToLL_M10to50_aMCatNLO",  1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("WJetsToLNu_aMCatNLO"    ,  1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("Data_SingleMu"    ,  1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("Data_SingleElec"    ,  1, 0, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_ScaleDown" ,  1, 0, true)'    
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_ScaleUp"   ,  1, 0, true)'    
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_Herwig", 1, 0, true)'



# FIDUCIAL
root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg"            , 4, 4, true)'
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_ScaleDown"  , 4, 5, true)'    
#root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_ScaleUp"    , 4, 5, true)'    
root -l -b -q 'RunTree_ReReco5TeV.C("TTbar_Powheg_Herwig"     , 4, 4, true)'
