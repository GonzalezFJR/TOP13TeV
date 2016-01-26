Download the analyzer
====

Download all the code from github.

    git clone https://github.com/GonzalezFJR/TOP13TeV

Download Dataset Manager from gitlab.

    git clone https://gitlab.cern.ch/IFCA-UO-CMS/Utils.git
    mv Utils/DatasetManager/ TOP13TeV
    cd TOP13TeV

Now you can set the enviroment and process one sample or directly run all the samples.


Set the enviroment
====

Load Root 6, start PAF and PoD.

If you run in Oviedo, to load Root 6 execute:

    source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh
    source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh
    source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh
    source /opt/root6/bin/thisroot.sh
    
Start Pod:

    source /opt/PoD/PoD_env.sh
  
Load PAF:

    source /nfs/fanae/PAF_releases/head/PAF_setup.sh
  
  
Process one sample
====

You have to choose the name of the sample, number of slots, if systematic studies are done (false by default) and the number of events (all in the sample by default).

    root -l -b -q 'RunTree_ReReco.C(sample, nSlots, doSyst?, nEvents)'

As an example:

    root -l -b -q 'RunTree_ReReco.C("TTJets", 1, true, 1000)'


Run all the samples
====



    source runAll_TopMC.sh

Root 6, PAF and PoD are loaded authomatically.
