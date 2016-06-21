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

#root -l -b -q 'RunTree_ReReco.C("TestHeppy", 1, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg", 15, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg", 15, 1, true)' # Semileptonic selection
#root -l -b -q 'RunTree_ReReco.C("ZZ"   , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("WW"   , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("WZ"   , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TW"   , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TbarW", 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_aMCatNLO"	 , 15, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M10to50_aMCatNLO_ext",  5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("WJetsToLNu_aMCatNLO"		 ,  5, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_scaleDown_ext", 15, 0, true)'    
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_scaleUp"  , 15, 0, true)'    
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1665" , 15, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1695" , 15, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1715" , 15, 0, true)'	
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1735" , 15, 0, true)'	
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1755" , 15, 0, true)'	
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1785" , 15, 0, true)'	
#root -l -b -q 'RunTree_ReReco.C("TTJets_aMCatNLO"	 , 15, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_Herwig"   , 15, 0, true)'	      
#root -l -b -q 'RunTree_ReReco.C("TTWToLNu"	  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTZToQQ"	  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTGJets"	  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TTWToQQ"	  , 5, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("MuonEG"    , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("DoubleMuon", 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("DoubleEG"  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("SingleEle" , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("SingleMu"  , 5, 0, true)'

#root -l -b -q 'RunTree_ReReco.C("TbarW_mtop1695"  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TbarW_mtop1755"  , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("TW_mtop1695"	   , 5, 0, true)'
#root -l -b -q 'RunTree_ReReco.C("Tree_TW_mtop1755", 5, 0, true)'

####root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_Herwig", 5, 0, true)'           
