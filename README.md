# TOP13TeV
# TOP Analyzer to run with the new PAF
# It uses the last version of DatasetManager, you should have the code in "DatasetManager/DatasetManager.h", or change this path in the RunTree_ReReco.C file
# Get the dataset manager: https://gitlab.cern.ch/IFCA-UO-CMS/Utils/tree/master/DatasetManager

# The analyzer and the other libraries should be in the folder "packages".
# Root 6 and the last version of PAF are used. Both are authomaticaly loaded in "runAll_TopMC.sh.

# Download all the files:
git clone https://github.com/GonzalezFJR/TOP13TeV
cd TOP13TeV
# Remember to copy the Dataset Manager header file

# Analyze all the samples
# source runAll_TopMC.sh

#         -- or --

# Initialise root 6 and PAF and analyze one sample:
  # Start root 6
  source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh;
  source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh;
  source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh;
  source /opt/root6/bin/thisroot.sh
  # Start PoD
  source /opt/PoD/PoD_env.sh
  # Start PAF
  source /nfs/fanae/PAF_releases/head/PAF_setup.sh

  root -l -b -q 'RunTree_ReReco.C("TTJets", 1, true, 1000)'
  
# root -l -b -q 'RunTree_ReReco.C(sample, nSlots, doSyst?, nEvents)'
# nEvets == 0 means all events in the sample
