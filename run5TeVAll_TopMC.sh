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

#root -l -b -q 'RunTree_ReReco.C("Test", 1, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg", 5, true)'
#root -l -b -q 'RunTree_ReReco.C("WW"   , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("WZ"   , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("TW"   , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("TbarW", 5, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_aMCatNLO",  5, true)'
#root -l -b -q 'RunTree_ReReco.C("WJetsToLNu_aMCatNLO"	  ,  5, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_scaleDown" , 5, true)'    
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_scaleUp"   , 5, true)'    
