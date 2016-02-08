# Start root 6
source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh; 
source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh;
source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh; 
source /opt/root6/bin/thisroot.sh
# Start PoD
source /opt/PoD/PoD_env.sh
# Start PAF
source /nfs/fanae/PAF_releases/head/PAF_setup.sh


#root -l -b -q 'RunTree_ReReco.C("TestHeppy", 1, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg", 30, true)'
root -l -b -q 'RunTree_ReReco.C("ZZ"	, 5, true)'
root -l -b -q 'RunTree_ReReco.C("WW"   , 5, true)'
root -l -b -q 'RunTree_ReReco.C("WZ"   , 5, true)'
root -l -b -q 'RunTree_ReReco.C("TW"   , 10, true)'
root -l -b -q 'RunTree_ReReco.C("TbarW" , 10, true)'
root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_aMCatNLO"    , 30, true)'
root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M10to50_aMCatNLO", 10, true)'
root -l -b -q 'RunTree_ReReco.C("WJetsToLNu_aMCatNLO", 10, true)'
#root -l -b -q 'RunTree_ReReco.C("TTJets_amcatnlo", 1, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg" 	 , 5, true)'    #/TT_TuneCUETP8M1_13TeV-powheg-pythia8                 /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2_ext3-v1/MINIAODSIM
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_Herwig"	, 20, true)'     #/TT_TuneEE5C_13TeV-powheg-herwigpp		       /RunIISpring15MiniAODv2-Asympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM        
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_ScaleDown", 20, true)'    #/TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8       /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM  
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_ScaleUp"  , 20, true)'    #/TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8	       /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM   
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1695" , 20, true)'    #/TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8        /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM   
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_mtop1755" , 20, true)'    #/TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8        /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM	
#root -l -b -q 'RunTree_ReReco.C("TTJets_aMCatNLO"       , 5, true)'    #/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8       /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM"
root -l -b -q 'RunTree_ReReco.C("MuonEG"	    , 10, true)'
root -l -b -q 'RunTree_ReReco.C("DoubleMuon"	    , 10, true)'
root -l -b -q 'RunTree_ReReco.C("DoubleEG"	    , 10, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M50_MLM"	     , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("DYJetsToLL_M10to50_MLM"     , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("WJetsToLNu_MLM"     , 5, true)'

#root -l -b -q 'RunTree_ReReco.C("T2tt_500_325"       , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_850_100"       , 5, true)'
#root -l -b -q 'RunTree_ReReco.C("TTWToLNu"      , 1, true)'    # /TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM
#root -l -b -q 'RunTree_ReReco.C("TTWToQQ"       , 5, true)'    # /TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8 /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM
#root -l -b -q 'RunTree_ReReco.C("TTZToQQ"       , 5, true)'    # /TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8                 /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM		
#root -l -b -q 'RunTree_ReReco.C("WWZ"           , 5, true)'    # /WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8                     /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM  		       
#root -l -b -q 'RunTree_ReReco.C("WZZ"           , 5, true)'    # /WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8                     /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM  		       
#root -l -b -q 'RunTree_ReReco.C("ZZZ"           , 5, true)'    # /ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8                     /RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM  		       

